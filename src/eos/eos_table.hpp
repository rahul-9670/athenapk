//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2026, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================
//! \file eos_table.hpp
//! \brief Tabulated equation of state for a hydrogen/helium mixture (H2 dissociation,
//!        H and He ionization, H2 rotational and vibrational modes).
//!
//! The table is generated offline by src/eos/gen_eos_table.py and read once at startup.
//! It holds, on evenly spaced log10 grids,
//!
//!     (log10 rho, log10 esp) -> P, cs^2, log10 T[K]
//!
//! where esp is the specific internal energy. The file stores everything in cgs and is
//! converted to code units on load with the run's own units, so one file works for any
//! normalization. Lookups are bilinear. The inverse e(rho, P), which the Riemann solver
//! needs, is a fixed-count bisection over the (monotone) esp axis.
//!
//! The struct only holds Kokkos View handles and a few scalars, so it is cheap to copy by
//! value into kernels.
//!
//! File layout (little endian, all row major):
//!     int64  magic, version, flags, nr, ne, nT
//!     double lr0, dlr, le0, dle, lT0, dlT, (two unused provenance slots)
//!     double P[nr][ne]      erg/cm^3
//!     double cs2[nr][ne]    cm^2/s^2
//!     double logT[nr][ne]   log10 K
//!     double esp[nr][nT]    erg/g, e(rho, T) on the (log10 rho, log10 T) grid
//! The axes are log10 rho [g/cm^3], log10 esp [erg/g] and log10 T [K]. The last block is
//! not used by the hydro solver and is skipped on load.
#ifndef EOS_EOS_TABLE_HPP_
#define EOS_EOS_TABLE_HPP_

#include <cmath>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

#include <Kokkos_Core.hpp>

#include "basic_types.hpp"
#include "parthenon_arrays.hpp"
#include "utils/error_checking.hpp"

namespace EOSTable {
using parthenon::ParArray2D;
using parthenon::Real;

//! First int64 of a table file, ASCII "EOSTABL1".
constexpr std::int64_t kMagic = 0x454F535441424C31LL;
constexpr std::int64_t kVersion = 2;
//! flags bit 0: axes and arrays are stored in cgs
constexpr std::int64_t kFlagCgs = 1;

struct EosTable {
  ParArray2D<Real> P_, cs2_, logT_;
  Real lr0_ = 0.0, dlr_ = 1.0, le0_ = 0.0, dle_ = 1.0;
  int nr_ = 0, ne_ = 0;
  //! Number of bisection steps in EspFromP. 20 steps over the ~7 dex esp axis of the
  //! production table resolve e(rho, P) to ~1e-5 relative, well below the interpolation
  //! error of the table itself. The count is fixed (no early exit) so that all threads of
  //! a warp do the same amount of work.
  int n_bisect_ = 20;
  bool loaded_ = false;

  //! Out-of-range lookups. Outside the table the lookup returns the value at the nearest
  //! edge (index and weight are both clamped). That is the safest numerical choice, but
  //! it must not go unnoticed, so every clamped lookup is counted here and reported once
  //! per cycle (see Hydro::PreStepMeshUserWorkInLoop). The counter is write-only from the
  //! kernels and never feeds back into the result.
  //! Layout: [0] rho below, [1] rho above, [2] esp below, [3] esp above.
  static constexpr int kNClampCtr = 4;
  Kokkos::View<unsigned long long *> clamp_ctr_;

  //! Bilinear interpolation of A[n1 x n2] on the axes (x0, dx) and (y0, dy), clamped to
  //! the table edges. Clamps are counted when `count` is true.
  KOKKOS_INLINE_FUNCTION
  Real Bilin(const ParArray2D<Real> &A, const Real x, const Real y,
             const bool count = true) const {
    const bool do_count = count && (clamp_ctr_.data() != nullptr);
    const Real fi = (x - lr0_) / dlr_;
    int i = static_cast<int>(fi);
    // Test the low side on fi, not on i: truncation towards zero maps fi in (-1, 0) to
    // i = 0, but that lookup is still clamped (via ti below).
    if (do_count && fi < 0.0) Kokkos::atomic_add(&clamp_ctr_(0), 1ull);
    if (i < 0) i = 0;
    if (i > nr_ - 2) {
      i = nr_ - 2;
      if (do_count) Kokkos::atomic_add(&clamp_ctr_(1), 1ull);
    }
    Real ti = fi - i;
    if (ti < 0.0) ti = 0.0;
    if (ti > 1.0) ti = 1.0;
    const Real fj = (y - le0_) / dle_;
    int j = static_cast<int>(fj);
    if (do_count && fj < 0.0) Kokkos::atomic_add(&clamp_ctr_(2), 1ull);
    if (j < 0) j = 0;
    if (j > ne_ - 2) {
      j = ne_ - 2;
      if (do_count) Kokkos::atomic_add(&clamp_ctr_(3), 1ull);
    }
    Real tj = fj - j;
    if (tj < 0.0) tj = 0.0;
    if (tj > 1.0) tj = 1.0;
    return (1.0 - ti) * (1.0 - tj) * A(i, j) + ti * (1.0 - tj) * A(i + 1, j) +
           (1.0 - ti) * tj * A(i, j + 1) + ti * tj * A(i + 1, j + 1);
  }

  //! Specific internal energy at (rho, P). P is monotone in esp at fixed rho, so bisect
  //! over the table's esp axis. The bisection never leaves the esp axis; the rho axis is
  //! checked once here, so an out-of-range rho is counted once and not n_bisect_ times.
  KOKKOS_INLINE_FUNCTION
  Real EspFromP(const Real rho, const Real P) const {
    const Real lr = std::log10(rho);
    if (clamp_ctr_.data() != nullptr) {
      const Real fi = (lr - lr0_) / dlr_;
      if (fi < 0.0) Kokkos::atomic_add(&clamp_ctr_(0), 1ull);
      if (static_cast<int>(fi) > nr_ - 2) Kokkos::atomic_add(&clamp_ctr_(1), 1ull);
    }
    Real lelo = le0_, lehi = le0_ + (ne_ - 1) * dle_;
    for (int it = 0; it < n_bisect_; ++it) {
      const Real lem = 0.5 * (lelo + lehi);
      if (Bilin(P_, lr, lem, false) < P) {
        lelo = lem;
      } else {
        lehi = lem;
      }
    }
    return std::pow(10.0, 0.5 * (lelo + lehi));
  }

  //! Gas pressure from density and internal energy DENSITY (code units).
  KOKKOS_INLINE_FUNCTION
  Real PresFromRhoEint(const Real rho, const Real eint) const {
    return Bilin(P_, std::log10(rho), std::log10(eint / rho));
  }
  //! Internal energy DENSITY from density and pressure (code units), inverse of the
  //! above.
  KOKKOS_INLINE_FUNCTION
  Real EintFromRhoPres(const Real rho, const Real pres) const {
    return EspFromP(rho, pres) * rho;
  }
  //! Adiabatic sound speed squared from density and pressure (code units).
  KOKKOS_INLINE_FUNCTION
  Real AsqFromRhoPres(const Real rho, const Real pres) const {
    // EspFromP already counted the rho axis, so only count the esp axis here.
    const Real esp = EspFromP(rho, pres);
    return Bilin(cs2_, std::log10(rho), std::log10(esp), false);
  }
  //! Gas temperature in K from density and internal energy DENSITY (code units).
  KOKKOS_INLINE_FUNCTION
  Real TemperatureK(const Real rho, const Real eint) const {
    return std::pow(10.0, Bilin(logT_, std::log10(rho), std::log10(eint / rho)));
  }

  //! Read the table written by gen_eos_table.py and convert it to code units.
  //! `rho_unit` [g/cm^3] and `v_unit` [cm/s] are the code units of this run.
  void Load(const std::string &path, const Real rho_unit, const Real v_unit) {
    std::ifstream f(path, std::ios::binary);
    PARTHENON_REQUIRE_THROWS(f.good(), "EOS table: cannot open " + path +
                                           ". Generate it with src/eos/gen_eos_table.py "
                                           "(see docs/eos.md).");
    std::int64_t iv[6];
    double dv[8];
    f.read(reinterpret_cast<char *>(iv), sizeof(iv));
    f.read(reinterpret_cast<char *>(dv), sizeof(dv));
    PARTHENON_REQUIRE_THROWS(f.good(), "EOS table: cannot read the header of " + path);
    PARTHENON_REQUIRE_THROWS(iv[0] == kMagic,
                             "EOS table: " + path +
                                 " is not an EOS table written by "
                                 "src/eos/gen_eos_table.py (bad magic).");
    PARTHENON_REQUIRE_THROWS(iv[1] == kVersion,
                             "EOS table: " + path + " has format version " +
                                 std::to_string(iv[1]) + ", this build reads version " +
                                 std::to_string(kVersion) + ".");
    PARTHENON_REQUIRE_THROWS((iv[2] & kFlagCgs) != 0,
                             "EOS table: " + path + " is not stored in cgs.");
    nr_ = static_cast<int>(iv[3]);
    ne_ = static_cast<int>(iv[4]);
    const int nT = static_cast<int>(iv[5]);
    lr0_ = dv[0];
    dlr_ = dv[1];
    le0_ = dv[2];
    dle_ = dv[3];
    PARTHENON_REQUIRE_THROWS(nr_ >= 2 && ne_ >= 2 && nT >= 2 && dlr_ > 0.0 && dle_ > 0.0,
                             "EOS table: bad grid in " + path);

    // The payload must be exactly three (nr x ne) arrays and one (nr x nT) array.
    const std::streamoff payload_start = f.tellg();
    f.seekg(0, std::ios::end);
    const std::streamoff nbytes = f.tellg() - payload_start;
    const std::streamoff expected = static_cast<std::streamoff>(sizeof(double)) *
                                    (3 * static_cast<std::streamoff>(nr_) * ne_ +
                                     static_cast<std::streamoff>(nr_) * nT);
    PARTHENON_REQUIRE_THROWS(nbytes == expected,
                             "EOS table: " + path + " has " + std::to_string(nbytes) +
                                 " payload bytes, expected " + std::to_string(expected) +
                                 " (truncated or corrupt file).");
    f.seekg(payload_start, std::ios::beg);

    // cgs to code units. The axes are log10 spaced, so a change of units only shifts the
    // axis origin; the spacing is unchanged. log10 T is unit free.
    const Real esp_unit = v_unit * v_unit;
    lr0_ -= std::log10(rho_unit);
    le0_ -= std::log10(esp_unit);

    auto load2d = [&](ParArray2D<Real> &view, const char *name, const Real scale) {
      view = ParArray2D<Real>(name, nr_, ne_);
      auto h = Kokkos::create_mirror_view(view);
      std::vector<double> buf(static_cast<size_t>(nr_) * ne_);
      f.read(reinterpret_cast<char *>(buf.data()),
             static_cast<std::streamsize>(buf.size() * sizeof(double)));
      for (int i = 0; i < nr_; ++i) {
        for (int j = 0; j < ne_; ++j) {
          h(i, j) = static_cast<Real>(buf[static_cast<size_t>(i) * ne_ + j]) * scale;
        }
      }
      Kokkos::deep_copy(view, h);
    };
    load2d(P_, "eos_P", 1.0 / (rho_unit * esp_unit));
    load2d(cs2_, "eos_cs2", 1.0 / esp_unit);
    load2d(logT_, "eos_logT", 1.0);
    PARTHENON_REQUIRE_THROWS(f.good(), "EOS table: error reading " + path);

    clamp_ctr_ = Kokkos::View<unsigned long long *>("eos_clamp_ctr", kNClampCtr);
    loaded_ = true;
  }

  //! Copy the clamp counts to the host and reset them. Host only, call at most once per
  //! cycle. Returns false if no table is loaded.
  bool ReadAndResetClampCounts(unsigned long long out[kNClampCtr]) const {
    if (clamp_ctr_.data() == nullptr) return false;
    auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), clamp_ctr_);
    for (int q = 0; q < kNClampCtr; ++q) {
      out[q] = h(q);
    }
    Kokkos::deep_copy(clamp_ctr_, 0ull);
    return true;
  }

  //! Axis ranges in code units (log10 rho and log10 esp), for messages.
  void AxisRanges(Real &lr_lo, Real &lr_hi, Real &le_lo, Real &le_hi) const {
    lr_lo = lr0_;
    lr_hi = lr0_ + (nr_ - 1) * dlr_;
    le_lo = le0_;
    le_hi = le0_ + (ne_ - 1) * dle_;
  }
};

} // namespace EOSTable

#endif // EOS_EOS_TABLE_HPP_
