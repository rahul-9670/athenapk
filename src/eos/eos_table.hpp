#ifndef EOS_EOS_TABLE_HPP_
#define EOS_EOS_TABLE_HPP_
//========================================================================================
// AthenaPK -- device-callable tabulated protostellar EOS. Holds 2D Kokkos tables generated
// offline by src/eos/gen_eos_table.py (full multi-Saha: H2 dissociation + H ionization +
// He/He+ ionization + H2 rot/vib + inert He), all in the shared FHC code units. On device
// it does O(1) bilinear interpolation on evenly-spaced log10 grids; the pressure inverse
// (needed by the Riemann solver and sound speed) is a short bisection over the monotone
// esp axis. The struct is small (holds Kokkos View handles + grid scalars) and is copied
// BY VALUE into kernels, so the handles ride along to the device correctly.
//
// Forward grid  (log10 rho_code, log10 esp_code) -> P_code, cs2_code, log10 T[K]
// RT grid       (log10 rho_code, log10 T[K])     -> esp_code
// esp = specific internal energy (code v0^2 units); e_int DENSITY = esp * rho.
//========================================================================================
#include <cmath>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Kokkos_Core.hpp"
#include "basic_types.hpp"       // parthenon::Real
#include "parthenon_arrays.hpp"  // parthenon::ParArray2D

#include "eos_table_format.hpp"  // shared table header parser (audit N14)

namespace EOSTable {
using parthenon::ParArray2D;
using parthenon::Real;

struct EosTable {
  ParArray2D<Real> P_, cs2_, logT_, espT_;
  Real lr0_ = 0.0, dlr_ = 1.0, le0_ = 0.0, dle_ = 1.0, lT0_ = 0.0, dlT_ = 1.0;
  int nr_ = 0, ne_ = 0, nT_ = 0;
  bool loaded_ = false;

  // bilinear interpolation of A[n1 x n2] over axes (x0,dx) and (y0,dy); edge-clamped.
  KOKKOS_INLINE_FUNCTION
  static Real bilin(const ParArray2D<Real> &A, int n1, int n2, Real x0, Real dx, Real y0,
                    Real dy, Real x, Real y) {
    Real fi = (x - x0) / dx;
    int i = static_cast<int>(fi);
    if (i < 0) i = 0;
    if (i > n1 - 2) i = n1 - 2;
    Real ti = fi - i;
    if (ti < 0.0) ti = 0.0;
    if (ti > 1.0) ti = 1.0;
    Real fj = (y - y0) / dy;
    int j = static_cast<int>(fj);
    if (j < 0) j = 0;
    if (j > n2 - 2) j = n2 - 2;
    Real tj = fj - j;
    if (tj < 0.0) tj = 0.0;
    if (tj > 1.0) tj = 1.0;
    return (1.0 - ti) * (1.0 - tj) * A(i, j) + ti * (1.0 - tj) * A(i + 1, j) +
           (1.0 - ti) * tj * A(i, j + 1) + ti * tj * A(i + 1, j + 1);
  }

  KOKKOS_INLINE_FUNCTION
  Real Pof(Real rho, Real esp) const {
    return bilin(P_, nr_, ne_, lr0_, dlr_, le0_, dle_, std::log10(rho), std::log10(esp));
  }
  KOKKOS_INLINE_FUNCTION
  Real Cs2of(Real rho, Real esp) const {
    return bilin(cs2_, nr_, ne_, lr0_, dlr_, le0_, dle_, std::log10(rho), std::log10(esp));
  }
  KOKKOS_INLINE_FUNCTION
  Real LogTof(Real rho, Real esp) const {
    return bilin(logT_, nr_, ne_, lr0_, dlr_, le0_, dle_, std::log10(rho), std::log10(esp));
  }
  // Invert P over log10(esp) at fixed rho (P is monotone increasing in esp).
  KOKKOS_INLINE_FUNCTION
  Real EspFromP(Real rho, Real P) const {
    const Real lr = std::log10(rho);
    Real lelo = le0_, lehi = le0_ + (ne_ - 1) * dle_;
    // 20 bisections => interval 2^-20 ~ 1e-6 of the log-esp range, still ~1e4x below the ~1%
    // table interpolation accuracy. FIXED count (not early-exit) is deliberate: on GPU a warp
    // runs until its slowest lane converges, so a fixed low count is faster + divergence-free than
    // an early-exit. 20 (was 30) trims this hot per-cell inversion (sound speed etc.) ~33%.
    for (int it = 0; it < 20; ++it) {
      const Real lem = 0.5 * (lelo + lehi);
      const Real Pm = bilin(P_, nr_, ne_, lr0_, dlr_, le0_, dle_, lr, lem);
      if (Pm < P) lelo = lem; else lehi = lem;
    }
    return std::pow(10.0, 0.5 * (lelo + lehi));
  }

  // ---- code-unit interface (mirrors the analytic HydrogenEOS names) ----
  // gas pressure from (rho, internal-energy DENSITY)
  KOKKOS_INLINE_FUNCTION
  Real PresFromRhoEint(Real rho, Real eint) const { return Pof(rho, eint / rho); }
  // internal-energy DENSITY from (rho, pressure)   [inverse of the above]
  KOKKOS_INLINE_FUNCTION
  Real EintFromRhoPres(Real rho, Real pres) const { return EspFromP(rho, pres) * rho; }
  // adiabatic sound speed^2 from (rho, pressure)
  KOKKOS_INLINE_FUNCTION
  Real AsqFromRhoPres(Real rho, Real pres) const { return Cs2of(rho, EspFromP(rho, pres)); }
  // gas temperature [K] from (rho, internal-energy DENSITY)
  KOKKOS_INLINE_FUNCTION
  Real TemperatureK(Real rho, Real eint) const {
    return std::pow(10.0, LogTof(rho, eint / rho));
  }
  // gas temperature [K] from (rho, pressure). Needed by the non-ideal MHD ionization model,
  // whose T_code = p/rho identity only holds for the ideal gas; with dissociation/ionization
  // mu varies, so T must come from the EOS. (esp inverse then the log10T table.)
  KOKKOS_INLINE_FUNCTION
  Real TemperatureKFromPres(Real rho, Real pres) const {
    const Real esp = EspFromP(rho, pres);
    return std::pow(10.0, bilin(logT_, nr_, ne_, lr0_, dlr_, le0_, dle_, std::log10(rho),
                                std::log10(esp)));
  }
  // internal-energy DENSITY from (rho, T[K])   [for RT]
  KOKKOS_INLINE_FUNCTION
  Real EintFromRhoTk(Real rho, Real Tk) const {
    return bilin(espT_, nr_, nT_, lr0_, dlr_, lT0_, dlT_, std::log10(rho), std::log10(Tk)) *
           rho;
  }

  //! True when the file was a legacy (v1) table, i.e. its axes and arrays were baked in the
  //! GENERATOR's code units rather than this run's. Lets the caller warn without re-parsing.
  bool used_legacy_code_units = false;

  // ---- host-side loader (reads the flat binary from gen_eos_table.py) ----
  //! `rho_unit` [g/cm^3] and `v_unit` [cm/s] are this run's authoritative scales. A v2 (cgs)
  //! table is converted into code units with them here, so the file is independent of any
  //! particular IC (audit N14). A legacy table is consumed verbatim -- bit-identical to the
  //! historical behaviour, and carrying the historical bias, which the caller reports.
  void Load(const std::string &path, const Real rho_unit, const Real v_unit) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("EosTable: cannot open table file " + path);
    const auto hdr = ReadEosHeader(f, path);
    nr_ = hdr.nr;
    ne_ = hdr.ne;
    nT_ = hdr.nT;
    lr0_ = hdr.lr0; dlr_ = hdr.dlr;
    le0_ = hdr.le0; dle_ = hdr.dle;
    lT0_ = hdr.lT0; dlT_ = hdr.dlT;
    if (nr_ < 2 || ne_ < 2 || nT_ < 2 || !(dlr_ > 0.0) || !(dle_ > 0.0) || !(dlT_ > 0.0))
      throw std::runtime_error("EosTable: bad grid in " + path);

    used_legacy_code_units = !hdr.in_cgs;
    // cgs -> code. The axes are log10-spaced, so this is a pure shift of the ORIGIN and the
    // SPACING is invariant; the payload arrays scale by fixed factors. log10 T[K] is
    // unit-free and is left alone. All factors are 1 / no-op for a legacy table.
    const Real esp_unit = v_unit * v_unit;         // erg/g
    const Real e_unit = rho_unit * esp_unit;       // erg/cm^3
    Real scale_P = 1.0, scale_cs2 = 1.0, scale_esp = 1.0;
    if (hdr.in_cgs) {
      lr0_ -= static_cast<Real>(std::log10(static_cast<double>(rho_unit)));
      le0_ -= static_cast<Real>(std::log10(static_cast<double>(esp_unit)));
      scale_P = 1.0 / e_unit;
      scale_cs2 = 1.0 / esp_unit;
      scale_esp = 1.0 / esp_unit;
    }

    auto load2d = [&](ParArray2D<Real> &view, const char *name, int n1, int n2,
                      const Real scale) {
      view = ParArray2D<Real>(name, n1, n2);
      auto h = Kokkos::create_mirror_view(view);
      std::vector<double> buf(static_cast<size_t>(n1) * n2);
      f.read(reinterpret_cast<char *>(buf.data()),
             static_cast<std::streamsize>(buf.size() * sizeof(double)));
      for (int i = 0; i < n1; ++i)
        for (int j = 0; j < n2; ++j)
          h(i, j) = static_cast<Real>(buf[i * n2 + j]) * scale;
      Kokkos::deep_copy(view, h);
    };
    load2d(P_, "eos_P", nr_, ne_, scale_P);
    load2d(cs2_, "eos_cs2", nr_, ne_, scale_cs2);
    load2d(logT_, "eos_logT", nr_, ne_, 1.0); // log10 T[K]: unit-free
    load2d(espT_, "eos_espT", nr_, nT_, scale_esp);
    if (!f) throw std::runtime_error("EosTable: truncated table file " + path);
    loaded_ = true;
  }
};

} // namespace EOSTable

#endif // EOS_EOS_TABLE_HPP_
