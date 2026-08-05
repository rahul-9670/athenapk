//========================================================================================
// AthenaPK - protostellar EOS table file format.
//
// AUDIT N14 (2026-08-05). The original `.bin` was headerless -- int64[3] then double[6] then
// four raw arrays -- and gen_eos_table.py built its axes and arrays in CODE units using its
// own hardcoded, rounded copies of rho0 and v0. The loader could not tell, because the file
// carried no unit field at all: `eos_table.bin` is exactly 1,238,472 bytes = 3*8 + 6*8 +
// 8*(180*220*3 + 180*200), with no spare byte where a unit could hide. The consequence was
// that a run whose authoritative PhysUnits scales differ from the generator's silently
// consulted the table at slightly wrong (rho, esp) -- +0.00315 % for the FHC IC.
//
// Version 2 stores everything in **cgs** and converts at load with the RUNNING simulation's
// units, so one file is correct for every normalization. The conversion is cheap and does not
// touch the device path: the axes are log10-spaced, so a unit change is a pure shift of the
// axis ORIGIN (the spacing is invariant), and the arrays scale by fixed factors.
//
//   log10 rho_code = log10 rho_cgs - log10(rho_unit)
//   log10 esp_code = log10 esp_cgs - log10(v_unit^2)
//   P_code   = P_cgs   / (rho_unit * v_unit^2)
//   cs2_code = cs2_cgs / v_unit^2
//   esp_code = esp_cgs / v_unit^2      (the RT grid's payload)
//   log10 T[K] is already unit-free.
//
// Legacy (v1) headerless files still load and behave EXACTLY as before, so every existing
// deck and queued job is unaffected. They warn (see hydro.cpp).
//========================================================================================
#ifndef EOS_EOS_TABLE_FORMAT_HPP_
#define EOS_EOS_TABLE_FORMAT_HPP_

#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>

namespace EOSTable {

//! First int64 of a v2+ file. Cannot collide with a legacy file, whose first int64 is `nr`
//! (a modest positive grid count). ASCII "EOSTABL1".
constexpr std::int64_t kEosTableMagic = 0x454F535441424C31LL;
constexpr std::int64_t kEosTableVersion = 2;

//! flags bit 0: axes and arrays are stored in CGS and must be converted at load.
constexpr std::int64_t kEosFlagCgs = 1;

struct EosFileHeader {
  int version = 1; //!< 1 = legacy headerless (code units), 2 = self-describing
  int nr = 0, ne = 0, nT = 0;
  double lr0 = 0.0, dlr = 1.0; //!< log10 rho axis  (cgs if in_cgs, else code)
  double le0 = 0.0, dle = 1.0; //!< log10 esp axis  (cgs if in_cgs, else code)
  double lT0 = 0.0, dlT = 1.0; //!< log10 T[K] axis (always unit-free)
  bool in_cgs = false;
  //! Provenance only, 0 when not applicable (a cgs table assumes no units at all).
  double gen_rho0 = 0.0, gen_v0 = 0.0;
};

//! Parse the header, leaving `f` positioned at the first byte of the payload. Detects the
//! version from the first int64 and rewinds for legacy files, so a v1 file reads exactly as
//! it always has.
inline EosFileHeader ReadEosHeader(std::ifstream &f, const std::string &path) {
  EosFileHeader h;
  std::int64_t first = 0;
  f.read(reinterpret_cast<char *>(&first), sizeof(first));
  if (!f) throw std::runtime_error("EosTable: cannot read header of " + path);

  if (first != kEosTableMagic) {
    f.seekg(0, std::ios::beg); // legacy: the int64 just consumed was nr
    std::int64_t hdr[3];
    f.read(reinterpret_cast<char *>(hdr), sizeof(hdr));
    double g[6];
    f.read(reinterpret_cast<char *>(g), sizeof(g));
    if (!f) throw std::runtime_error("EosTable: truncated legacy header in " + path);
    h.version = 1;
    h.nr = static_cast<int>(hdr[0]);
    h.ne = static_cast<int>(hdr[1]);
    h.nT = static_cast<int>(hdr[2]);
    h.lr0 = g[0]; h.dlr = g[1];
    h.le0 = g[2]; h.dle = g[3];
    h.lT0 = g[4]; h.dlT = g[5];
    h.in_cgs = false; // v1 always stored code units -- this is the N14 defect
    return h;
  }

  std::int64_t iv[5]; // version, flags, nr, ne, nT
  f.read(reinterpret_cast<char *>(iv), sizeof(iv));
  double dv[8]; // lr0, dlr, le0, dle, lT0, dlT, gen_rho0, gen_v0
  f.read(reinterpret_cast<char *>(dv), sizeof(dv));
  if (!f) throw std::runtime_error("EosTable: truncated v2 header in " + path);
  h.version = static_cast<int>(iv[0]);
  if (h.version > kEosTableVersion)
    throw std::runtime_error("EosTable: " + path + " is format version " +
                             std::to_string(h.version) + ", newer than this build supports (" +
                             std::to_string(kEosTableVersion) +
                             "). Rebuild, or regenerate the table.");
  h.in_cgs = (iv[1] & kEosFlagCgs) != 0;
  h.nr = static_cast<int>(iv[2]);
  h.ne = static_cast<int>(iv[3]);
  h.nT = static_cast<int>(iv[4]);
  h.lr0 = dv[0]; h.dlr = dv[1];
  h.le0 = dv[2]; h.dle = dv[3];
  h.lT0 = dv[4]; h.dlT = dv[5];
  h.gen_rho0 = dv[6];
  h.gen_v0 = dv[7];
  return h;
}

} // namespace EOSTable

#endif // EOS_EOS_TABLE_FORMAT_HPP_
