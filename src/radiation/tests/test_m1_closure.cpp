// Flagship Phase 4 (radiation): characterization gate for the SHIPPED M1 closure.
//
// Compiles the REAL src/radiation/radiation_closure.hpp host-only (pure device math; the tiny
// host_shims/ stub <Kokkos_Core.hpp>/<Kokkos_Pair.hpp>/<basic_types.hpp>) and asserts the M1
// two-moment closure identities the transport relies on. This characterizes the gray-M1
// baseline; MULTIGROUP RHD (the Phase-4 deliverable) is not yet implemented -> deferred, large
// new-numerics effort for user direction.
//
// Build + run (host g++, no Kokkos/MPI):
//   cd src/radiation/tests
//   g++ -O2 -std=c++17 -Ihost_shims test_m1_closure.cpp -o /tmp/t && /tmp/t
//
// Gates (chi := ThriceEddingtonFactor = 3 * Eddington factor chi_Edd):
//   1. Eddington-factor limits: chi_Edd(f=0)=1/3 (diffusion), chi_Edd(f=1)=1 (free-streaming).
//   2. Eddington tensor trace == E for ANY reduced-flux vector (radiation-energy consistency).
//   3. limits: f->1 beam gives P = n_i n_j (P_par=1, P_perp=0); f->0 gives P = (1/3) delta.
//   4. M1 wave speeds are causal (|lambda| <= 1 in chat units); f=0 -> +-1/sqrt3;
//      f=1,mu=1 -> +1 (beam at c).
//   5. NormalizeFlux clamps |f| <= 1 preserving direction.
#include <array>
#include <cmath>
#include <cstdio>

#include "../radiation_closure.hpp"

using namespace Radiation;

static int fails = 0;
static void check(const char *name, bool ok, double val, double ref) {
  std::printf("  %-52s %s (got %.6e, ref %.6e)\n", name, ok ? "PASS" : "FAIL", val, ref);
  if (!ok) ++fails;
}
static bool close(double a, double b, double tol = 1e-12) {
  return std::fabs(a - b) <= tol * (std::fabs(b) + 1e-300) + tol;
}

int main() {
  const double t = 1e-12;

  std::printf("--- GATE 1: Eddington-factor limits ---\n");
  const double chi0 = ThriceEddingtonFactor<Closure::M1>(0.0) / 3.0;
  const double chi1 = ThriceEddingtonFactor<Closure::M1>(1.0) / 3.0;
  check("chi_Edd(f=0) == 1/3", close(chi0, 1.0 / 3.0), chi0, 1.0 / 3.0);
  check("chi_Edd(f=1) == 1", close(chi1, 1.0), chi1, 1.0);

  std::printf("--- GATE 2: Eddington tensor trace == E (any flux vector) ---\n");
  double worst_tr = 0.0;
  const std::array<std::array<double, 3>, 5> fv = {{{0, 0, 0}, {0.3, 0, 0},
      {0.2, -0.5, 0.1}, {0.6, 0.6, 0.4}, {0.99, 0.0, 0.0}}};
  for (auto &fr : fv) {
    auto P = EddingtonTensor<Closure::M1>(fr);
    worst_tr = std::max(worst_tr, std::fabs(P[0] + P[1] + P[2] - 1.0));
  }
  check("max |trace(P/E) - 1|", worst_tr < t, worst_tr, 0.0);

  std::printf("--- GATE 3: beam (f->1) and diffusion (f->0) tensor limits ---\n");
  auto Pbeam = EddingtonTensor<Closure::M1>({1.0, 0.0, 0.0});   // beam along x
  bool beam_ok = close(Pbeam[0], 1.0, 1e-6) && close(Pbeam[1], 0.0, 1e-6) &&
                 close(Pbeam[2], 0.0, 1e-6) && close(Pbeam[5], 0.0, 1e-6);
  check("f=1 beam: P_xx=1, P_yy=P_zz=P_xy=0", beam_ok, Pbeam[0], 1.0);
  auto Piso = EddingtonTensor<Closure::M1>({0.0, 0.0, 0.0});
  bool iso_ok = close(Piso[0], 1.0 / 3.0) && close(Piso[1], 1.0 / 3.0) &&
                close(Piso[2], 1.0 / 3.0) && close(Piso[3], 0.0);
  check("f=0 diffusion: P = (1/3) delta", iso_ok, Piso[0], 1.0 / 3.0);

  std::printf("--- GATE 4: causal M1 wave speeds ---\n");
  double worst_sp = 0.0;
  for (double f = 0.0; f <= 1.0001; f += 0.1)
    for (double mu = -1.0; mu <= 1.0001; mu += 0.25) {
      auto s = WaveSpeed<Closure::M1>(mu, std::min(f, 1.0));
      worst_sp = std::max(worst_sp, std::max(std::fabs(s.first), std::fabs(s.second)));
    }
  check("max |lambda| <= 1 (causal)", worst_sp <= 1.0 + 1e-9, worst_sp, 1.0);
  auto sdiff = WaveSpeed<Closure::M1>(0.0, 0.0);
  check("f=0: |lambda| == 1/sqrt(3)", close(sdiff.second, std::sqrt(1.0 / 3.0), 1e-9),
        sdiff.second, std::sqrt(1.0 / 3.0));
  auto sbeam = WaveSpeed<Closure::M1>(1.0, 1.0);
  check("f=1,mu=1: lambda_max == 1 (beam at c)", close(sbeam.second, 1.0, 1e-6), sbeam.second,
        1.0);

  std::printf("--- GATE 5: NormalizeFlux clamps |f|<=1, preserves direction ---\n");
  auto nf = NormalizeFlux(3.0, 4.0, 0.0); // |f|=5 -> clamp to 1 along (0.6,0.8,0)
  double fmag = std::sqrt(nf[0] * nf[0] + nf[1] * nf[1] + nf[2] * nf[2]);
  bool dir_ok = close(nf[0], 0.6, 1e-9) && close(nf[1], 0.8, 1e-9);
  check("clamp |f|=5 -> 1, direction (0.6,0.8,0)", close(fmag, 1.0, 1e-9) && dir_ok, fmag, 1.0);

  std::printf("\n%s (%d failures)\n", fails == 0 ? "ALL GATES PASS" : "GATE FAILURES", fails);
  return fails ? 1 : 0;
}
