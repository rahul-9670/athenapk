//========================================================================================
// AthenaPK - sink-particle package (WS-1 of PHYSICS_COMPLETION_PLAN.md)
// Swarm plumbing modeled on src/tracers/tracers.cpp (Athena-Parthenon / phoebus).
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#include "basic_types.hpp"
#include "globals.hpp"
#include "interface/metadata.hpp"
#include "kokkos_abstraction.hpp"
#include "utils/error_checking.hpp"

#include "../main.hpp"
#include "sinks.hpp"

namespace Sinks {
using namespace parthenon::package::prelude;

std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin) {
  auto pkg = std::make_shared<StateDescriptor>("sinks");
  const bool enabled = pin->GetOrAddBoolean("sinks", "enabled", false);
  pkg->AddParam<bool>("enabled", enabled);
  // OFF-state: register nothing else (no swarm, no work function) so a run without a
  // <sinks> block is bit-identical to the pre-sinks code.
  if (!enabled) return pkg;

  const std::string swarm_name = "sinks";
  pkg->AddParam<>("swarm_name", swarm_name);

  // Swarm carries built-in {id, x, y, z} plus the fields below. Restart metadata makes
  // Parthenon checkpoint/restore the swarm (positions + values) in the .rhdf restart.
  Metadata swarm_metadata({Metadata::Provides, Metadata::None, Metadata::Restart});
  pkg->AddSwarm(swarm_name, swarm_metadata);
  Metadata real_value({Metadata::Real, Metadata::Restart});
  for (const auto &v : {"mass", "vx", "vy", "vz", "Lx", "Ly", "Lz", "t_created"}) {
    pkg->AddSwarmValue(v, swarm_name, real_value);
  }

  PARTHENON_REQUIRE_THROWS(pin->GetInteger("parthenon/mesh", "nx3") > 1,
                           "Sinks currently only supported/tested in 3D.");

  // --- Gravity of sinks on the gas (WS-1 increment 2) ------------------------
  // G in code units: G = four_pi_G / (4 pi), matching the self-gravity convention.
  const Real four_pi_G = pin->GetOrAddReal("sinks", "four_pi_G", 1.0);
  pkg->AddParam("four_pi_G", four_pi_G);
  // Plummer softening length = soft_cells * dx_min (avoids the r->0 singularity).
  const Real soft_cells = pin->GetOrAddReal("sinks", "soft_cells", 2.0);
  pkg->AddParam("soft_cells", soft_cells);
  // Whether the sinks' gravity acts on the gas (WS-1 inc 2). Off isolates the N-body (inc 3).
  pkg->AddParam("gas_gravity", pin->GetOrAddBoolean("sinks", "gas_gravity", true));
  // Whether sinks feel each other's gravity (WS-1 inc 3 two-body dynamics).
  pkg->AddParam("sink_gravity", pin->GetOrAddBoolean("sinks", "sink_gravity", true));
  // Sink-sink Plummer softening (code length; 0 = pure Newtonian, e.g. for Kepler orbits).
  pkg->AddParam("soft_sink", pin->GetOrAddReal("sinks", "soft_sink", 0.0));
  // KDK-leapfrog subcycle factor: substep dt_sub = subcycle_cfl * min_i(|v_i|/|a_i|).
  pkg->AddParam("subcycle_cfl", pin->GetOrAddReal("sinks", "subcycle_cfl", 0.1));

  // --- Sink creation (WS-1 increment 4) --------------------------------------
  pkg->AddParam("creation", pin->GetOrAddBoolean("sinks", "creation", false));
  // Density threshold: rho_sink_code >= 0 uses that fixed value; < 0 => per-cell Truelove
  // rho where the local Jeans length = n_jeans_sink * dx (cell becomes Jeans-unresolved).
  pkg->AddParam("rho_sink_code", pin->GetOrAddReal("sinks", "rho_sink_code", -1.0));
  // Truelove resolution factor for the auto rho_sink (default 8; set to match refinement/njeans).
  pkg->AddParam("n_jeans_sink", pin->GetOrAddReal("sinks", "n_jeans_sink", 8.0));
  // Accretion radius (control radius) in cells; the no-nearby-sink test uses 2*racc.
  pkg->AddParam("racc_cells", pin->GetOrAddReal("sinks", "racc_cells", 4.0));
  pkg->AddParam("gamma", pin->GetReal("hydro", "gamma"));
  // New sinks get ids starting past the seeded ones (seeds use ids 0..nseed-1).
  pkg->AddParam("next_sink_id", pin->GetOrAddInteger("sinks", "nseed", 1),
                /*is_mutable=*/true);
  // Accretion (WS-1 increment 5): remove bound gas within r_acc down to a rho_sink/3 floor.
  pkg->AddParam("accretion", pin->GetOrAddBoolean("sinks", "accretion", false));
  // Every cell feels every sink, so all ranks need the full sink list each step.
  const int max_sinks = pin->GetOrAddInteger("sinks", "max_sinks", 64);
  pkg->AddParam("max_sinks", max_sinks);
  // Device-resident {x,y,z,mass} of every sink, refilled by GatherSinks each step.
  Kokkos::View<Real *[4]> sink_data("Sinks::sink_data", max_sinks);
  pkg->AddParam("sink_data", sink_data);
  pkg->AddParam("sink_count", 0, /*is_mutable=*/true);

  pkg->UserWorkBeforeLoopMesh = SeedInitialSinks;
  return pkg;
}

// Collect {x,y,z,mass} of all active sinks across all blocks and ranks into the package's
// device sink_data View (replicated on every rank). Runs once per step before the gas source.
// Assumes packs_per_rank==1 so md holds every local block.
TaskStatus GatherSinks(MeshData<Real> *md) {
  auto pm = md->GetParentPointer();
  auto pkg = pm->packages.Get("sinks");
  const int max_sinks = pkg->Param<int>("max_sinks");

  std::vector<Real> local; // flattened {x,y,z,mass} per local sink
  for (auto &pmb : pm->block_list) {
    auto &swarm = pmb->meshblock_data.Get()->GetSwarmData()->Get("sinks");
    const int maxi = swarm->GetMaxActiveIndex();
    if (maxi < 0) continue;
    auto mask_h = swarm->GetMask().GetHostMirrorAndCopy();
    auto x_h = swarm->Get<Real>(swarm_position::x::name()).Get().GetHostMirrorAndCopy();
    auto y_h = swarm->Get<Real>(swarm_position::y::name()).Get().GetHostMirrorAndCopy();
    auto z_h = swarm->Get<Real>(swarm_position::z::name()).Get().GetHostMirrorAndCopy();
    auto m_h = swarm->Get<Real>("mass").Get().GetHostMirrorAndCopy();
    for (int n = 0; n <= maxi; ++n) {
      if (!mask_h(n)) continue;
      local.push_back(x_h(n));
      local.push_back(y_h(n));
      local.push_back(z_h(n));
      local.push_back(m_h(n));
    }
  }

  std::vector<Real> global = local;
#ifdef MPI_PARALLEL
  int nranks = 1;
  PARTHENON_MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &nranks));
  int myn = static_cast<int>(local.size());
  std::vector<int> counts(nranks), displs(nranks);
  PARTHENON_MPI_CHECK(
      MPI_Allgather(&myn, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD));
  int total = 0;
  for (int r = 0; r < nranks; ++r) {
    displs[r] = total;
    total += counts[r];
  }
  global.resize(total);
  PARTHENON_MPI_CHECK(MPI_Allgatherv(local.data(), myn, MPI_PARTHENON_REAL, global.data(),
                                     counts.data(), displs.data(), MPI_PARTHENON_REAL,
                                     MPI_COMM_WORLD));
#endif
  const int nsink = static_cast<int>(global.size() / 4);
  PARTHENON_REQUIRE(nsink <= max_sinks, "Number of sinks exceeds sinks/max_sinks.");

  auto sink_data = pkg->Param<Kokkos::View<Real *[4]>>("sink_data");
  auto sink_data_h = Kokkos::create_mirror_view(sink_data);
  for (int s = 0; s < nsink; ++s)
    for (int c = 0; c < 4; ++c) sink_data_h(s, c) = global[4 * s + c];
  Kokkos::deep_copy(sink_data, sink_data_h);
  pkg->UpdateParam("sink_count", nsink);
  return TaskStatus::complete;
}

// Add the sinks' point-mass gravity to the gas as an operator-split source on the final
// stage (mirrors SelfGravity::ApplyGravitySource). Momentum kick dm = rho*g*dt with Plummer
// softening; energy gets the exact resulting kinetic-energy change (conservative).
TaskStatus ApplySinkGravity(MeshData<Real> *md, const parthenon::SimTime &, const Real dt) {
  auto pm = md->GetParentPointer();
  auto pkg = pm->packages.Get("sinks");
  if (!pkg->Param<bool>("gas_gravity")) return TaskStatus::complete;
  const int nsink = pkg->Param<int>("sink_count");
  if (nsink == 0) return TaskStatus::complete;

  auto sink_data = pkg->Param<Kokkos::View<Real *[4]>>("sink_data");
  const Real four_pi_G = pkg->Param<Real>("four_pi_G");
  const Real G = four_pi_G / (4.0 * M_PI);
  const Real soft_cells = pkg->Param<Real>("soft_cells");

  const auto &cons_pack = md->PackVariables(std::vector<std::string>{"cons"});
  const auto &prim_pack = md->PackVariables(std::vector<std::string>{"prim"});

  IndexRange ib = md->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBoundsK(IndexDomain::interior);
  const int nblocks = md->NumBlocks();

  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "Sinks::ApplySinkGravity", parthenon::DevExecSpace(), 0,
      nblocks - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        const auto &coords = prim_pack.GetCoords(b);
        const Real xc = coords.template Xc<parthenon::X1DIR>(k, j, i);
        const Real yc = coords.template Xc<parthenon::X2DIR>(k, j, i);
        const Real zc = coords.template Xc<parthenon::X3DIR>(k, j, i);
        const Real eps = soft_cells * coords.template Dxc<parthenon::X1DIR>(k, j, i);
        const Real eps2 = eps * eps;

        Real gx = 0.0, gy = 0.0, gz = 0.0;
        for (int s = 0; s < nsink; ++s) {
          const Real dx = xc - sink_data(s, 0);
          const Real dy = yc - sink_data(s, 1);
          const Real dz = zc - sink_data(s, 2);
          const Real r2 = dx * dx + dy * dy + dz * dz + eps2;
          const Real inv_r3 = 1.0 / (r2 * std::sqrt(r2));
          const Real fac = -G * sink_data(s, 3) * inv_r3;
          gx += fac * dx;
          gy += fac * dy;
          gz += fac * dz;
        }

        auto &cons = cons_pack(b);
        const Real rho = cons(IDN, k, j, i);
        const Real mx = cons(IM1, k, j, i), my = cons(IM2, k, j, i), mz = cons(IM3, k, j, i);
        const Real ke_old = 0.5 * (mx * mx + my * my + mz * mz) / rho;
        const Real nmx = mx + rho * gx * dt;
        const Real nmy = my + rho * gy * dt;
        const Real nmz = mz + rho * gz * dt;
        cons(IM1, k, j, i) = nmx;
        cons(IM2, k, j, i) = nmy;
        cons(IM3, k, j, i) = nmz;
        const Real ke_new = 0.5 * (nmx * nmx + nmy * nmy + nmz * nmz) / rho;
        cons(IEN, k, j, i) += ke_new - ke_old;
      });
  return TaskStatus::complete;
}

void SeedInitialSinks(Mesh *pmesh, ParameterInput *pin, parthenon::SimTime &tm) {
  // On restart the swarm is restored from the restart file; do not re-seed.
  if (parthenon::Globals::is_restart) return;

  // nseed==1 uses the legacy flat keys (seed_x,...); nseed>1 uses seed{i}_x,... .
  const int nseed = pin->GetOrAddInteger("sinks", "nseed", 1);
  for (int si = 0; si < nseed; ++si) {
    const std::string p =
        (nseed == 1) ? std::string("seed_") : ("seed" + std::to_string(si) + "_");
    const Real x0 = pin->GetReal("sinks", p + "x");
    const Real y0 = pin->GetReal("sinks", p + "y");
    const Real z0 = pin->GetReal("sinks", p + "z");
    const Real vx0 = pin->GetOrAddReal("sinks", p + "vx", 0.0);
    const Real vy0 = pin->GetOrAddReal("sinks", p + "vy", 0.0);
    const Real vz0 = pin->GetOrAddReal("sinks", p + "vz", 0.0);
    const Real m0 = pin->GetOrAddReal("sinks", p + "mass", 0.0);
    const std::uint64_t sid = static_cast<std::uint64_t>(si);

    // Place seed sink si on whichever block owns (x0,y0,z0). Exactly one block on exactly
    // one rank matches (half-open cell-face intervals avoid double placement).
    for (auto &pmb : pmesh->block_list) {
      IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
      IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
      IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);
      const Real xlo = pmb->coords.Xf<1>(ib.s), xhi = pmb->coords.Xf<1>(ib.e + 1);
      const Real ylo = pmb->coords.Xf<2>(jb.s), yhi = pmb->coords.Xf<2>(jb.e + 1);
      const Real zlo = pmb->coords.Xf<3>(kb.s), zhi = pmb->coords.Xf<3>(kb.e + 1);
      if (!(x0 >= xlo && x0 < xhi && y0 >= ylo && y0 < yhi && z0 >= zlo && z0 < zhi))
        continue;

      auto &swarm = pmb->meshblock_data.Get()->GetSwarmData()->Get("sinks");
      auto new_particles_context = swarm->AddEmptyParticles(1);

      auto &x = swarm->Get<Real>(swarm_position::x::name()).Get();
      auto &y = swarm->Get<Real>(swarm_position::y::name()).Get();
      auto &z = swarm->Get<Real>(swarm_position::z::name()).Get();
      auto &id = swarm->Get<std::uint64_t>(swarm_position::id::name()).Get();
      auto &mass = swarm->Get<Real>("mass").Get();
      auto &vx = swarm->Get<Real>("vx").Get();
      auto &vy = swarm->Get<Real>("vy").Get();
      auto &vz = swarm->Get<Real>("vz").Get();
      auto &Lx = swarm->Get<Real>("Lx").Get();
      auto &Ly = swarm->Get<Real>("Ly").Get();
      auto &Lz = swarm->Get<Real>("Lz").Get();
      auto &t_created = swarm->Get<Real>("t_created").Get();
      auto swarm_d = swarm->GetDeviceContext();

      pmb->par_for(
          "Sinks::SeedInitialSinks", 0, new_particles_context.GetNewParticlesMaxIndex(),
          KOKKOS_LAMBDA(const int new_n) {
            const int n = new_particles_context.GetNewParticleIndex(new_n);
            x(n) = x0; y(n) = y0; z(n) = z0;
            id(n) = sid;
            mass(n) = m0;
            vx(n) = vx0; vy(n) = vy0; vz(n) = vz0;
            Lx(n) = 0.0; Ly(n) = 0.0; Lz(n) = 0.0;
            t_created(n) = 0.0;
            bool on_current_mesh_block = true;
            swarm_d.GetNeighborBlockIndex(n, x(n), y(n), z(n), on_current_mesh_block);
          });
    }
  }
}

// Sink creation (WS-1 increment 4). Once per step, find the single densest cell that meets
// all formation criteria and spawn one sink there. Criteria (Federrath-style, minus the full
// virial term): (1) rho > rho_sink [fixed, or per-cell Truelove where the Jeans length =
// n_jeans_sink*dx]; (2) local density maximum over the 6 face neighbours; (3) converging flow
// div(v) < 0; (5) no existing sink within 2*r_acc. Picking only the global maximum guarantees
// at most one sink forms per step; criterion (5) prevents a second sink forming on the pile-up
// around a just-formed one. No gas is removed yet (mass conservation is increment 5).
TaskStatus CreateSinks(MeshData<Real> *md, const parthenon::SimTime &tm) {
  auto pm = md->GetParentPointer();
  auto pkg = pm->packages.Get("sinks");
  if (!pkg->Param<bool>("creation")) return TaskStatus::complete;

  const Real G = pkg->Param<Real>("four_pi_G") / (4.0 * M_PI);
  const Real gam = pkg->Param<Real>("gamma");
  const Real njeans = pkg->Param<Real>("n_jeans_sink");
  const Real rho_sink_fixed = pkg->Param<Real>("rho_sink_code");
  const Real racc_cells = pkg->Param<Real>("racc_cells");
  const int nsink = pkg->Param<int>("sink_count");
  auto sink_data = pkg->Param<Kokkos::View<Real *[4]>>("sink_data");

  const auto &prim = md->PackVariables(std::vector<std::string>{"prim"});
  IndexRange ib = md->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBoundsK(IndexDomain::interior);
  const int nb = md->NumBlocks();
  const int NI = ib.e - ib.s + 1, NJ = jb.e - jb.s + 1, NK = kb.e - kb.s + 1;
  const long NCELL = static_cast<long>(NI) * NJ * NK;

  using MaxLocReal = Kokkos::MaxLoc<Real, long>;
  using mlv = MaxLocReal::value_type;
  mlv res;
  Kokkos::parallel_reduce(
      "Sinks::FindCandidate", Kokkos::RangePolicy<>(0, static_cast<long>(nb) * NCELL),
      KOKKOS_LAMBDA(const long idx, mlv &lm) {
        const int b = static_cast<int>(idx / NCELL);
        const long r = idx - static_cast<long>(b) * NCELL;
        const int k = kb.s + static_cast<int>(r / (NI * NJ));
        const int j = jb.s + static_cast<int>((r / NI) % NJ);
        const int i = ib.s + static_cast<int>(r % NI);
        const auto prm = prim(b);
        const auto &coords = prim.GetCoords(b);
        const Real rho = prm(IDN, k, j, i);

        // (1) density threshold
        const Real dx = coords.template Dxc<parthenon::X1DIR>(k, j, i);
        const Real dy = coords.template Dxc<parthenon::X2DIR>(k, j, i);
        const Real dz = coords.template Dxc<parthenon::X3DIR>(k, j, i);
        Real rho_sink = rho_sink_fixed;
        if (rho_sink_fixed < 0.0) {
          const Real cs2 = gam * prm(IPR, k, j, i) / rho;
          rho_sink = M_PI * cs2 / (G * njeans * dx * njeans * dx);
        }
        if (rho < rho_sink) return;
        // (2) local density maximum (6 face neighbours)
        if (rho < prm(IDN, k, j, i - 1) || rho < prm(IDN, k, j, i + 1) ||
            rho < prm(IDN, k, j - 1, i) || rho < prm(IDN, k, j + 1, i) ||
            rho < prm(IDN, k - 1, j, i) || rho < prm(IDN, k + 1, j, i))
          return;
        // (3) converging flow
        const Real divv = (prm(IV1, k, j, i + 1) - prm(IV1, k, j, i - 1)) / (2.0 * dx) +
                          (prm(IV2, k, j + 1, i) - prm(IV2, k, j - 1, i)) / (2.0 * dy) +
                          (prm(IV3, k + 1, j, i) - prm(IV3, k - 1, j, i)) / (2.0 * dz);
        if (divv >= 0.0) return;
        // (5) no existing sink within 2*r_acc
        const Real xc = coords.template Xc<parthenon::X1DIR>(k, j, i);
        const Real yc = coords.template Xc<parthenon::X2DIR>(k, j, i);
        const Real zc = coords.template Xc<parthenon::X3DIR>(k, j, i);
        const Real r2min = (2.0 * racc_cells * dx) * (2.0 * racc_cells * dx);
        for (int s = 0; s < nsink; ++s) {
          const Real dxs = xc - sink_data(s, 0), dys = yc - sink_data(s, 1),
                     dzs = zc - sink_data(s, 2);
          if (dxs * dxs + dys * dys + dzs * dzs < r2min) return;
        }
        if (rho > lm.val) { lm.val = rho; lm.loc = idx; }
      },
      MaxLocReal(res));
  Kokkos::fence();

  // Global owner of the densest candidate via MAXLOC over (density, rank).
  const bool have_local = (res.loc >= 0);
  int owner = parthenon::Globals::my_rank;
  Real gmax = have_local ? res.val : -std::numeric_limits<Real>::max();
#ifdef MPI_PARALLEL
  struct {
    double val;
    int rank;
  } lin, gout;
  lin.val = have_local ? static_cast<double>(res.val) : -std::numeric_limits<double>::max();
  lin.rank = parthenon::Globals::my_rank;
  PARTHENON_MPI_CHECK(MPI_Allreduce(&lin, &gout, 1, MPI_DOUBLE_INT, MPI_MAXLOC,
                                    MPI_COMM_WORLD));
  gmax = gout.val;
  owner = gout.rank;
#endif
  if (gmax <= 0.0) return TaskStatus::complete; // no candidate anywhere
  if (parthenon::Globals::my_rank != owner || !have_local) return TaskStatus::complete;

  // Decode the winning cell and extract its {x,y,z,vx,vy,vz,rho*dV} to the host.
  const int b = static_cast<int>(res.loc / NCELL);
  const long r = res.loc - static_cast<long>(b) * NCELL;
  const int k = kb.s + static_cast<int>(r / (NI * NJ));
  const int j = jb.s + static_cast<int>((r / NI) % NJ);
  const int i = ib.s + static_cast<int>(r % NI);
  Kokkos::View<Real[7]> cell("Sinks::new_sink_cell");
  Kokkos::parallel_for(
      "Sinks::ExtractCell", 1, KOKKOS_LAMBDA(const int) {
        const auto prm = prim(b);
        const auto &coords = prim.GetCoords(b);
        cell(0) = coords.template Xc<parthenon::X1DIR>(k, j, i);
        cell(1) = coords.template Xc<parthenon::X2DIR>(k, j, i);
        cell(2) = coords.template Xc<parthenon::X3DIR>(k, j, i);
        cell(3) = prm(IV1, k, j, i);
        cell(4) = prm(IV2, k, j, i);
        cell(5) = prm(IV3, k, j, i);
        cell(6) = prm(IDN, k, j, i) * coords.CellVolume(k, j, i);
      });
  auto cell_h = Kokkos::create_mirror_view(cell);
  Kokkos::deep_copy(cell_h, cell);
  const Real cx = cell_h(0), cy = cell_h(1), cz = cell_h(2);
  const Real cvx = cell_h(3), cvy = cell_h(4), cvz = cell_h(5), cmass = cell_h(6);
  const Real ct = tm.time;
  const std::uint64_t newid = static_cast<std::uint64_t>(pkg->Param<int>("next_sink_id"));

  auto *pmb = md->GetBlockData(b)->GetBlockPointer();
  auto &swarm = pmb->meshblock_data.Get()->GetSwarmData()->Get("sinks");
  auto npc = swarm->AddEmptyParticles(1);
  auto &x = swarm->Get<Real>(swarm_position::x::name()).Get();
  auto &y = swarm->Get<Real>(swarm_position::y::name()).Get();
  auto &z = swarm->Get<Real>(swarm_position::z::name()).Get();
  auto &id = swarm->Get<std::uint64_t>(swarm_position::id::name()).Get();
  auto &mass = swarm->Get<Real>("mass").Get();
  auto &vx = swarm->Get<Real>("vx").Get();
  auto &vy = swarm->Get<Real>("vy").Get();
  auto &vz = swarm->Get<Real>("vz").Get();
  auto &Lx = swarm->Get<Real>("Lx").Get();
  auto &Ly = swarm->Get<Real>("Ly").Get();
  auto &Lz = swarm->Get<Real>("Lz").Get();
  auto &t_created = swarm->Get<Real>("t_created").Get();
  auto swarm_d = swarm->GetDeviceContext();
  pmb->par_for(
      "Sinks::SpawnSink", 0, npc.GetNewParticlesMaxIndex(), KOKKOS_LAMBDA(const int nn) {
        const int n = npc.GetNewParticleIndex(nn);
        x(n) = cx; y(n) = cy; z(n) = cz;
        id(n) = newid;
        mass(n) = cmass;
        vx(n) = cvx; vy(n) = cvy; vz(n) = cvz;
        Lx(n) = 0.0; Ly(n) = 0.0; Lz(n) = 0.0;
        t_created(n) = ct;
        bool on = true;
        swarm_d.GetNeighborBlockIndex(n, x(n), y(n), z(n), on);
      });
  pkg->UpdateParam("next_sink_id", static_cast<int>(newid) + 1);
  if (parthenon::Globals::my_rank == 0)
    std::cout << "[sinks] created sink id=" << newid << " at (" << cx << "," << cy << ","
              << cz << ") rho*dV=" << cmass << " t=" << ct << std::endl;
  return TaskStatus::complete;
}

// Sink accretion (WS-1 increment 5). For gas cells within r_acc of a sink that are bound to it
// and above a rho_sink/3 floor, remove a quadratic-ramp fraction of the mass (and, proportional
// to that, momentum + energy; B is NOT accreted) and add it to the sink -- conserving total
// mass and momentum. The rho_sink/3 floor caps the density in the sink region, which is what
// lets the timestep recover past first-core formation. Self-contained gather (own {id,x,v,mass}
// Allgatherv), so it does not depend on GatherSinks. Non-overlapping r_acc (creation enforces
// 2*r_acc spacing) => each cell accretes to at most one sink.
TaskStatus AccreteSinks(MeshData<Real> *md, const parthenon::SimTime &tm, const Real dt) {
  auto pm = md->GetParentPointer();
  auto pkg = pm->packages.Get("sinks");
  if (!pkg->Param<bool>("accretion")) return TaskStatus::complete;
  const Real G = pkg->Param<Real>("four_pi_G") / (4.0 * M_PI);
  const Real gam = pkg->Param<Real>("gamma");
  const Real njeans = pkg->Param<Real>("n_jeans_sink");
  const Real rho_sink_fixed = pkg->Param<Real>("rho_sink_code");
  const Real racc_cells = pkg->Param<Real>("racc_cells");

  // --- Gather global {id, x,y,z, vx,vy,vz, mass} (same pattern as AdvanceSinksNBody) ---
  std::vector<std::uint64_t> lid;
  std::vector<Real> ld;
  for (auto &pmb : pm->block_list) {
    auto &sw = pmb->meshblock_data.Get()->GetSwarmData()->Get("sinks");
    const int maxi = sw->GetMaxActiveIndex();
    if (maxi < 0) continue;
    auto mask_h = sw->GetMask().GetHostMirrorAndCopy();
    auto xh = sw->Get<Real>(swarm_position::x::name()).Get().GetHostMirrorAndCopy();
    auto yh = sw->Get<Real>(swarm_position::y::name()).Get().GetHostMirrorAndCopy();
    auto zh = sw->Get<Real>(swarm_position::z::name()).Get().GetHostMirrorAndCopy();
    auto idh =
        sw->Get<std::uint64_t>(swarm_position::id::name()).Get().GetHostMirrorAndCopy();
    auto vxh = sw->Get<Real>("vx").Get().GetHostMirrorAndCopy();
    auto vyh = sw->Get<Real>("vy").Get().GetHostMirrorAndCopy();
    auto vzh = sw->Get<Real>("vz").Get().GetHostMirrorAndCopy();
    auto mh = sw->Get<Real>("mass").Get().GetHostMirrorAndCopy();
    for (int n = 0; n <= maxi; ++n) {
      if (!mask_h(n)) continue;
      lid.push_back(idh(n));
      ld.push_back(xh(n));  ld.push_back(yh(n));  ld.push_back(zh(n));
      ld.push_back(vxh(n)); ld.push_back(vyh(n)); ld.push_back(vzh(n));
      ld.push_back(mh(n));
    }
  }
  std::vector<std::uint64_t> gid = lid;
  std::vector<Real> gd = ld;
#ifdef MPI_PARALLEL
  int nranks = 1;
  PARTHENON_MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &nranks));
  int nloc = static_cast<int>(lid.size());
  std::vector<int> cnt(nranks), dsp(nranks), cnt7(nranks), dsp7(nranks);
  PARTHENON_MPI_CHECK(
      MPI_Allgather(&nloc, 1, MPI_INT, cnt.data(), 1, MPI_INT, MPI_COMM_WORLD));
  int tot = 0;
  for (int r = 0; r < nranks; ++r) {
    dsp[r] = tot; cnt7[r] = cnt[r] * 7; dsp7[r] = tot * 7; tot += cnt[r];
  }
  gid.resize(tot); gd.resize(tot * 7);
  PARTHENON_MPI_CHECK(MPI_Allgatherv(lid.data(), nloc, MPI_UINT64_T, gid.data(), cnt.data(),
                                     dsp.data(), MPI_UINT64_T, MPI_COMM_WORLD));
  PARTHENON_MPI_CHECK(MPI_Allgatherv(ld.data(), nloc * 7, MPI_PARTHENON_REAL, gd.data(),
                                     cnt7.data(), dsp7.data(), MPI_PARTHENON_REAL,
                                     MPI_COMM_WORLD));
#endif
  const int ns = static_cast<int>(gid.size());
  if (ns == 0) return TaskStatus::complete;

  Kokkos::View<Real *[7]> state("Sinks::acc_state", ns); // x,y,z,vx,vy,vz,mass
  auto state_h = Kokkos::create_mirror_view(state);
  for (int i = 0; i < ns; ++i)
    for (int c = 0; c < 7; ++c) state_h(i, c) = gd[7 * i + c];
  Kokkos::deep_copy(state, state_h);

  // Per-sink accumulators {dM, dpx,dpy,dpz, dLx,dLy,dLz}.
  Kokkos::View<Real **> accum("Sinks::acc", ns, 7);
  Kokkos::deep_copy(accum, 0.0);

  const auto &cons = md->PackVariables(std::vector<std::string>{"cons"});
  IndexRange ib = md->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBoundsK(IndexDomain::interior);
  const int nb = md->NumBlocks();

  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "Sinks::Accrete", parthenon::DevExecSpace(), 0, nb - 1, kb.s,
      kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int b, const int k, const int j, const int i) {
        auto &c = cons(b);
        const auto &coords = cons.GetCoords(b);
        const Real dx = coords.template Dxc<parthenon::X1DIR>(k, j, i);
        const Real dV = coords.CellVolume(k, j, i);
        const Real racc = racc_cells * dx;
        const Real xc = coords.template Xc<parthenon::X1DIR>(k, j, i);
        const Real yc = coords.template Xc<parthenon::X2DIR>(k, j, i);
        const Real zc = coords.template Xc<parthenon::X3DIR>(k, j, i);
        const Real rho = c(IDN, k, j, i);
        // per-cell floor rho_sink/3 (fixed threshold, or Truelove)
        Real rho_sink = rho_sink_fixed;
        if (rho_sink_fixed < 0.0) {
          // c(IEN) is total energy; approximate p via (gamma-1)(E - KE) for the sound speed
          const Real ke = 0.5 *
                          (c(IM1, k, j, i) * c(IM1, k, j, i) + c(IM2, k, j, i) * c(IM2, k, j, i) +
                           c(IM3, k, j, i) * c(IM3, k, j, i)) / rho;
          const Real p = (gam - 1.0) * (c(IEN, k, j, i) - ke);
          const Real cs2 = gam * p / rho;
          rho_sink = M_PI * cs2 / (G * njeans * dx * njeans * dx);
        }
        const Real floor = rho_sink / 3.0;
        if (rho <= floor) return;
        for (int s = 0; s < ns; ++s) {
          const Real dxs = xc - state(s, 0), dys = yc - state(s, 1), dzs = zc - state(s, 2);
          const Real r2 = dxs * dxs + dys * dys + dzs * dzs;
          if (r2 >= racc * racc) continue;
          const Real r = std::sqrt(r2);
          const Real vgx = c(IM1, k, j, i) / rho, vgy = c(IM2, k, j, i) / rho,
                     vgz = c(IM3, k, j, i) / rho;
          const Real dvx = vgx - state(s, 3), dvy = vgy - state(s, 4), dvz = vgz - state(s, 5);
          const Real dv2 = dvx * dvx + dvy * dvy + dvz * dvz;
          // bound: 0.5 |v_rel|^2 - G m_s / r < 0
          if (0.5 * dv2 - G * state(s, 6) / (r + 1e-30) >= 0.0) continue;
          const Real ramp = ((racc - r) / racc) * ((racc - r) / racc);
          const Real dM = ramp * (rho - floor) * dV; // total mass removed (>=0)
          const Real phi = dM / (rho * dV);           // fraction removed (<=1 since ramp<=1)
          // proportional removal of density, momentum, energy (B untouched for MHD)
          c(IDN, k, j, i) = rho * (1.0 - phi);
          const Real dpx = c(IM1, k, j, i) * phi, dpy = c(IM2, k, j, i) * phi,
                     dpz = c(IM3, k, j, i) * phi;
          c(IM1, k, j, i) -= dpx; c(IM2, k, j, i) -= dpy; c(IM3, k, j, i) -= dpz;
          c(IEN, k, j, i) *= (1.0 - phi);
          // accumulate onto the sink (dp already a total momentum: density*phi*... no -> *dV)
          Kokkos::atomic_add(&accum(s, 0), dM);
          Kokkos::atomic_add(&accum(s, 1), dpx * dV);
          Kokkos::atomic_add(&accum(s, 2), dpy * dV);
          Kokkos::atomic_add(&accum(s, 3), dpz * dV);
          // angular momentum of extracted gas about the sink: dM (r_rel x v_rel)
          Kokkos::atomic_add(&accum(s, 4), dM * (dys * dvz - dzs * dvy));
          Kokkos::atomic_add(&accum(s, 5), dM * (dzs * dvx - dxs * dvz));
          Kokkos::atomic_add(&accum(s, 6), dM * (dxs * dvy - dys * dvx));
          break; // at most one sink per cell
        }
      });
  Kokkos::fence();

  auto accum_h = Kokkos::create_mirror_view(accum);
  Kokkos::deep_copy(accum_h, accum);
  std::vector<Real> ac(ns * 7);
  for (int s = 0; s < ns; ++s)
    for (int cc = 0; cc < 7; ++cc) ac[7 * s + cc] = accum_h(s, cc);
#ifdef MPI_PARALLEL
  PARTHENON_MPI_CHECK(MPI_Allreduce(MPI_IN_PLACE, ac.data(), ns * 7, MPI_PARTHENON_REAL,
                                    MPI_SUM, MPI_COMM_WORLD));
#endif

  // Update each rank's own sinks (momentum-conserving): m += dM, v = (m v + dp)/(m+dM).
  std::unordered_map<std::uint64_t, int> id2g;
  for (int s = 0; s < ns; ++s) id2g[gid[s]] = s;
  for (auto &pmb : pm->block_list) {
    auto &sw = pmb->meshblock_data.Get()->GetSwarmData()->Get("sinks");
    const int maxi = sw->GetMaxActiveIndex();
    if (maxi < 0) continue;
    auto &md_ = sw->Get<Real>("mass").Get();
    auto &vxd = sw->Get<Real>("vx").Get();
    auto &vyd = sw->Get<Real>("vy").Get();
    auto &vzd = sw->Get<Real>("vz").Get();
    auto &Lxd = sw->Get<Real>("Lx").Get();
    auto &Lyd = sw->Get<Real>("Ly").Get();
    auto &Lzd = sw->Get<Real>("Lz").Get();
    auto mask_h = sw->GetMask().GetHostMirrorAndCopy();
    auto idh = sw->Get<std::uint64_t>(swarm_position::id::name()).Get().GetHostMirrorAndCopy();
    auto mh = Kokkos::create_mirror_view(md_), vxh = Kokkos::create_mirror_view(vxd),
         vyh = Kokkos::create_mirror_view(vyd), vzh = Kokkos::create_mirror_view(vzd),
         Lxh = Kokkos::create_mirror_view(Lxd), Lyh = Kokkos::create_mirror_view(Lyd),
         Lzh = Kokkos::create_mirror_view(Lzd);
    Kokkos::deep_copy(mh, md_); Kokkos::deep_copy(vxh, vxd); Kokkos::deep_copy(vyh, vyd);
    Kokkos::deep_copy(vzh, vzd); Kokkos::deep_copy(Lxh, Lxd); Kokkos::deep_copy(Lyh, Lyd);
    Kokkos::deep_copy(Lzh, Lzd);
    for (int n = 0; n <= maxi; ++n) {
      if (!mask_h(n)) continue;
      const int s = id2g[idh(n)];
      const Real dM = ac[7 * s + 0];
      if (dM <= 0.0) continue;
      const Real m0 = mh(n), mnew = m0 + dM;
      vxh(n) = (m0 * vxh(n) + ac[7 * s + 1]) / mnew;
      vyh(n) = (m0 * vyh(n) + ac[7 * s + 2]) / mnew;
      vzh(n) = (m0 * vzh(n) + ac[7 * s + 3]) / mnew;
      mh(n) = mnew;
      Lxh(n) += ac[7 * s + 4]; Lyh(n) += ac[7 * s + 5]; Lzh(n) += ac[7 * s + 6];
    }
    Kokkos::deep_copy(md_, mh); Kokkos::deep_copy(vxd, vxh); Kokkos::deep_copy(vyd, vyh);
    Kokkos::deep_copy(vzd, vzh); Kokkos::deep_copy(Lxd, Lxh); Kokkos::deep_copy(Lyd, Lyh);
    Kokkos::deep_copy(Lzd, Lzh);
  }
  return TaskStatus::complete;
}

// Mesh-level sink advance: gather every sink {id,x,v,mass} onto all ranks, integrate the
// whole (small-N) system with a subcycled KDK leapfrog on the host (deterministic, so every
// rank gets identical results), then write each rank's own sinks back and update their
// neighbor-block index for migration. With sink_gravity off or a single sink the acceleration
// is zero, so KDK reduces to exact constant-velocity drift (WS-1 increment 1 ballistic path).
TaskStatus AdvanceSinksNBody(MeshData<Real> *md, const Real dt) {
  auto pm = md->GetParentPointer();
  auto pkg = pm->packages.Get("sinks");
  const bool sink_gravity = pkg->Param<bool>("sink_gravity");
  const Real G = pkg->Param<Real>("four_pi_G") / (4.0 * M_PI);
  const Real eps = pkg->Param<Real>("soft_sink");
  const Real eps2 = eps * eps;
  const Real subcycle_cfl = pkg->Param<Real>("subcycle_cfl");

  // --- 1. Gather global {id, x,y,z, vx,vy,vz, mass} onto every rank -----------
  std::vector<std::uint64_t> lid;
  std::vector<Real> ld; // 7 Reals per sink
  for (auto &pmb : pm->block_list) {
    auto &swarm = pmb->meshblock_data.Get()->GetSwarmData()->Get("sinks");
    const int maxi = swarm->GetMaxActiveIndex();
    if (maxi < 0) continue;
    auto mask_h = swarm->GetMask().GetHostMirrorAndCopy();
    auto x_h = swarm->Get<Real>(swarm_position::x::name()).Get().GetHostMirrorAndCopy();
    auto y_h = swarm->Get<Real>(swarm_position::y::name()).Get().GetHostMirrorAndCopy();
    auto z_h = swarm->Get<Real>(swarm_position::z::name()).Get().GetHostMirrorAndCopy();
    auto id_h =
        swarm->Get<std::uint64_t>(swarm_position::id::name()).Get().GetHostMirrorAndCopy();
    auto vx_h = swarm->Get<Real>("vx").Get().GetHostMirrorAndCopy();
    auto vy_h = swarm->Get<Real>("vy").Get().GetHostMirrorAndCopy();
    auto vz_h = swarm->Get<Real>("vz").Get().GetHostMirrorAndCopy();
    auto m_h = swarm->Get<Real>("mass").Get().GetHostMirrorAndCopy();
    for (int n = 0; n <= maxi; ++n) {
      if (!mask_h(n)) continue;
      lid.push_back(id_h(n));
      ld.push_back(x_h(n));  ld.push_back(y_h(n));  ld.push_back(z_h(n));
      ld.push_back(vx_h(n)); ld.push_back(vy_h(n)); ld.push_back(vz_h(n));
      ld.push_back(m_h(n));
    }
  }
  std::vector<std::uint64_t> gid = lid;
  std::vector<Real> gd = ld;
#ifdef MPI_PARALLEL
  int nranks = 1;
  PARTHENON_MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &nranks));
  int nloc = static_cast<int>(lid.size());
  std::vector<int> cnt(nranks), dsp(nranks), cnt7(nranks), dsp7(nranks);
  PARTHENON_MPI_CHECK(
      MPI_Allgather(&nloc, 1, MPI_INT, cnt.data(), 1, MPI_INT, MPI_COMM_WORLD));
  int tot = 0;
  for (int r = 0; r < nranks; ++r) {
    dsp[r] = tot;   cnt7[r] = cnt[r] * 7;   dsp7[r] = tot * 7;   tot += cnt[r];
  }
  gid.resize(tot);
  gd.resize(tot * 7);
  PARTHENON_MPI_CHECK(MPI_Allgatherv(lid.data(), nloc, MPI_UINT64_T, gid.data(), cnt.data(),
                                     dsp.data(), MPI_UINT64_T, MPI_COMM_WORLD));
  PARTHENON_MPI_CHECK(MPI_Allgatherv(ld.data(), nloc * 7, MPI_PARTHENON_REAL, gd.data(),
                                     cnt7.data(), dsp7.data(), MPI_PARTHENON_REAL,
                                     MPI_COMM_WORLD));
#endif
  const int ns = static_cast<int>(gid.size());
  if (ns == 0) return TaskStatus::complete;

  // Unpack into position/velocity/mass working arrays.
  std::vector<Real> x(ns), y(ns), z(ns), vx(ns), vy(ns), vz(ns), m(ns);
  for (int i = 0; i < ns; ++i) {
    x[i] = gd[7 * i + 0]; y[i] = gd[7 * i + 1]; z[i] = gd[7 * i + 2];
    vx[i] = gd[7 * i + 3]; vy[i] = gd[7 * i + 4]; vz[i] = gd[7 * i + 5];
    m[i] = gd[7 * i + 6];
  }

  // --- 2. Subcycled KDK leapfrog (host; identical on every rank) --------------
  std::vector<Real> ax(ns), ay(ns), az(ns);
  auto accel = [&]() {
    for (int i = 0; i < ns; ++i) { ax[i] = ay[i] = az[i] = 0.0; }
    if (!sink_gravity) return;
    for (int i = 0; i < ns; ++i) {
      for (int j = 0; j < ns; ++j) {
        if (j == i) continue;
        const Real dx = x[j] - x[i], dy = y[j] - y[i], dz = z[j] - z[i];
        const Real r2 = dx * dx + dy * dy + dz * dz + eps2;
        const Real f = G * m[j] / (r2 * std::sqrt(r2));
        ax[i] += f * dx; ay[i] += f * dy; az[i] += f * dz;
      }
    }
  };
  accel();
  // Subcycle count from the shortest v/a timescale (guard a==0 -> no force -> nsub=1).
  int nsub = 1;
  if (subcycle_cfl > 0.0) {
    Real tmin = std::numeric_limits<Real>::max();
    for (int i = 0; i < ns; ++i) {
      const Real amag = std::sqrt(ax[i] * ax[i] + ay[i] * ay[i] + az[i] * az[i]);
      const Real vmag = std::sqrt(vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
      if (amag > 0.0) tmin = std::min(tmin, vmag / amag);
    }
    if (tmin < std::numeric_limits<Real>::max())
      nsub = std::max(1, static_cast<int>(std::ceil(dt / (subcycle_cfl * tmin))));
  }
  const Real dts = dt / nsub;
  for (int s = 0; s < nsub; ++s) {
    for (int i = 0; i < ns; ++i) { // kick
      vx[i] += 0.5 * dts * ax[i]; vy[i] += 0.5 * dts * ay[i]; vz[i] += 0.5 * dts * az[i];
    }
    for (int i = 0; i < ns; ++i) { // drift
      x[i] += dts * vx[i]; y[i] += dts * vy[i]; z[i] += dts * vz[i];
    }
    accel();
    for (int i = 0; i < ns; ++i) { // kick
      vx[i] += 0.5 * dts * ax[i]; vy[i] += 0.5 * dts * ay[i]; vz[i] += 0.5 * dts * az[i];
    }
  }

  // --- 3. Write each rank's own sinks back (matched by id) + migrate ----------
  std::unordered_map<std::uint64_t, int> id2g;
  for (int i = 0; i < ns; ++i) id2g[gid[i]] = i;
  for (auto &pmb : pm->block_list) {
    auto &swarm = pmb->meshblock_data.Get()->GetSwarmData()->Get("sinks");
    const int maxi = swarm->GetMaxActiveIndex();
    if (maxi < 0) continue;
    auto &xd = swarm->Get<Real>(swarm_position::x::name()).Get();
    auto &yd = swarm->Get<Real>(swarm_position::y::name()).Get();
    auto &zd = swarm->Get<Real>(swarm_position::z::name()).Get();
    auto &vxd = swarm->Get<Real>("vx").Get();
    auto &vyd = swarm->Get<Real>("vy").Get();
    auto &vzd = swarm->Get<Real>("vz").Get();
    auto &idd = swarm->Get<std::uint64_t>(swarm_position::id::name()).Get();
    auto mask_h = swarm->GetMask().GetHostMirrorAndCopy();
    auto xh = Kokkos::create_mirror_view(xd), yh = Kokkos::create_mirror_view(yd),
         zh = Kokkos::create_mirror_view(zd), vxh = Kokkos::create_mirror_view(vxd),
         vyh = Kokkos::create_mirror_view(vyd), vzh = Kokkos::create_mirror_view(vzd);
    auto idh = swarm->Get<std::uint64_t>(swarm_position::id::name()).Get().GetHostMirrorAndCopy();
    for (int n = 0; n <= maxi; ++n) {
      if (!mask_h(n)) continue;
      const int g = id2g[idh(n)];
      xh(n) = x[g]; yh(n) = y[g]; zh(n) = z[g];
      vxh(n) = vx[g]; vyh(n) = vy[g]; vzh(n) = vz[g];
    }
    Kokkos::deep_copy(xd, xh); Kokkos::deep_copy(yd, yh); Kokkos::deep_copy(zd, zh);
    Kokkos::deep_copy(vxd, vxh); Kokkos::deep_copy(vyd, vyh); Kokkos::deep_copy(vzd, vzh);
    auto swarm_d = swarm->GetDeviceContext();
    pmb->par_for(
        "Sinks::MigrateSinks", 0, maxi, KOKKOS_LAMBDA(const int n) {
          if (swarm_d.IsActive(n)) {
            bool on_current_mesh_block = true;
            swarm_d.GetNeighborBlockIndex(n, xd(n), yd(n), zd(n), on_current_mesh_block);
          }
        });
  }
  return TaskStatus::complete;
}

} // namespace Sinks
