//========================================================================================
// AthenaPK - opacity table file format (shared by the two readers).
//
// AUDIT N13 (2026-08-05). The original format was headerless -- int64[3] then double[4] then
// the payload -- and gen_opacity_table.py wrote kappa ALREADY CONVERTED to code units using
// its own hardcoded, rounded copies of rho0 and l0. Neither reader rescaled, so every
// tabulated kappa carried a fixed bias (+0.135 % against the exact BE normalization) that no
// C++ guard could see, because the offending constants lived in an offline Python script.
//
// Version 2 fixes the cause rather than the symptom: kappa is stored in **cgs**, and the
// loader multiplies by the RUNNING simulation's opacity_unit. A v2 table is therefore
// independent of any particular IC -- the same file is correct for every normalization.
//
// Legacy (v1) headerless files are still readable and behave EXACTLY as before, so every
// existing deck and every queued job is unaffected. They warn.
//
// Both `OpacityTable::Load` (radiation_opacity.hpp) and `BuildGroupOpacityTableFromFile`
// (radiation_groups.hpp) parse the header THROUGH THIS ONE FUNCTION. They used to carry
// independent copies of the layout, which is a standing invitation for them to drift apart.
//========================================================================================
#ifndef RADIATION_OPACITY_TABLE_FORMAT_HPP_
#define RADIATION_OPACITY_TABLE_FORMAT_HPP_

#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>

namespace Radiation {

//! First int64 of a v2+ file. Cannot collide with a legacy file, whose first int64 is `ng`
//! (a small positive group count, <= MAX_GROUP). ASCII "OPACTBL1".
constexpr std::int64_t kOpacityTableMagic = 0x4F50414354424C31LL;
constexpr std::int64_t kOpacityTableVersion = 2;

//! flags bit 0: gray kappa arrays (kP, kR, ks) are stored in CGS [cm^2/g] and MUST be
//! multiplied by the run's opacity_unit at load. Clear => stored in code units (the v1
//! behaviour, retained only so a v2 file could in principle be written the old way).
constexpr std::int64_t kOpacityFlagKappaCgs = 1;

struct OpacityFileHeader {
  int version = 1; //!< 1 = legacy headerless, 2 = self-describing
  int ng = 1, nr = 0, nT = 0;
  double lr0 = 0.0, dlr = 1.0; //!< log10(rho_cgs) axis origin + spacing
  double lT0 = 0.0, dlT = 1.0; //!< log10(T_K)     axis origin + spacing
  bool kappa_in_cgs = false;   //!< true => loader must scale by the runtime opacity_unit
  //! Provenance only, 0 when not applicable (a cgs table assumes no units at all).
  double gen_rho_unit = 0.0, gen_len_unit = 0.0, gen_opacity_unit = 0.0;
};

//! Parse the header, leaving `f` positioned at the first byte of the payload. Detects the
//! version from the first int64 and rewinds for legacy files, so a v1 file reads exactly as
//! it always has.
inline OpacityFileHeader ReadOpacityHeader(std::ifstream &f, const std::string &path) {
  OpacityFileHeader h;
  std::int64_t first = 0;
  f.read(reinterpret_cast<char *>(&first), sizeof(first));
  if (!f) throw std::runtime_error("OpacityTable: cannot read header of " + path);

  if (first != kOpacityTableMagic) {
    // Legacy v1: the int64 just consumed was ng. Rewind and read the original layout.
    f.seekg(0, std::ios::beg);
    std::int64_t hdr[3];
    f.read(reinterpret_cast<char *>(hdr), sizeof(hdr));
    double g[4];
    f.read(reinterpret_cast<char *>(g), sizeof(g));
    if (!f) throw std::runtime_error("OpacityTable: truncated legacy header in " + path);
    h.version = 1;
    h.ng = static_cast<int>(hdr[0]);
    h.nr = static_cast<int>(hdr[1]);
    h.nT = static_cast<int>(hdr[2]);
    h.lr0 = g[0]; h.dlr = g[1]; h.lT0 = g[2]; h.dlT = g[3];
    h.kappa_in_cgs = false; // v1 always stored code units -- this is the N13 bias
    return h;
  }

  std::int64_t iv[5]; // version, flags, ng, nr, nT
  f.read(reinterpret_cast<char *>(iv), sizeof(iv));
  double dv[7]; // lr0, dlr, lT0, dlT, gen_rho_unit, gen_len_unit, gen_opacity_unit
  f.read(reinterpret_cast<char *>(dv), sizeof(dv));
  if (!f) throw std::runtime_error("OpacityTable: truncated v2 header in " + path);
  h.version = static_cast<int>(iv[0]);
  if (h.version > kOpacityTableVersion)
    throw std::runtime_error("OpacityTable: " + path + " is format version " +
                             std::to_string(h.version) + ", newer than this build supports (" +
                             std::to_string(kOpacityTableVersion) +
                             "). Rebuild, or regenerate the table.");
  h.kappa_in_cgs = (iv[1] & kOpacityFlagKappaCgs) != 0;
  h.ng = static_cast<int>(iv[2]);
  h.nr = static_cast<int>(iv[3]);
  h.nT = static_cast<int>(iv[4]);
  h.lr0 = dv[0]; h.dlr = dv[1]; h.lT0 = dv[2]; h.dlT = dv[3];
  h.gen_rho_unit = dv[4];
  h.gen_len_unit = dv[5];
  h.gen_opacity_unit = dv[6];
  return h;
}

} // namespace Radiation

#endif // RADIATION_OPACITY_TABLE_FORMAT_HPP_
