# Equation of state

The equation of state is selected with `eos` in the `<hydro>` block.

- `adiabatic`: ideal gas with adiabatic index `gamma` (default choice in all example
  inputs).
- `hydrogen`: tabulated EOS of a hydrogen/helium mixture that includes H2 dissociation,
  H and He ionization and the rotational and vibrational modes of H2. Meant for star
  formation problems where the gas is heated beyond ~1000 K, e.g., the second collapse
  in protostellar cores.

## Tabulated hydrogen EOS

### Physics

The table is computed offline by `src/eos/gen_eos_table.py` from coupled Saha equilibria
of H2, HI, HII, HeI, HeII, HeIII and electrons, for a hydrogen mass fraction `X = 0.716`
(helium `Y = 0.284`, mean molecular weight 2.33 for molecular gas). The internal energy
contains the translational energy, the H2 rotational energy for a frozen 3:1 ortho:para
mixture, the harmonic H2 vibrational energy, and the dissociation and ionization energies
(zero point: cold molecular gas). The effective adiabatic index is therefore ~5/3 below
~100 K, ~7/5 for warm molecular gas, and drops below 4/3 during H2 dissociation
(~2000 K) and H ionization (~1e4 K). `gen_eos_table.py --check` prints it as a function
of temperature.

The file holds, on evenly spaced log10 grids, pressure, adiabatic sound speed squared and
temperature as functions of `(log10 rho, log10 esp)`, where `esp` is the specific internal
energy. Everything is stored in cgs and converted to code units at startup using the
`<units>` block, so the same file can be used with any normalization. The code
interpolates bilinearly. The inverse `e(rho, P)`, needed by the Riemann solver and the
sound speed, is a fixed number of bisection steps over the esp axis.

### Usage

```
<units>
code_length_cgs = 2.81e16             # all three are required with eos = hydrogen
code_mass_cgs   = 1.21302020147e31
code_time_cgs   = 1478947368421.0527

<hydro>
fluid          = glmmhd
eos            = hydrogen
eos_table_file = /path/to/eos_table_hires_v5_cold.bin
gamma          = 1.666666666666667    # only used by problem generators (initial energy)
riemann        = hlld
#eos_bisections = 20                  # bisection steps for e(rho, P), in [10, 60]
#eos_table_clamp_report = true        # report lookups outside the table every cycle
```

The table is not part of the repository and has to be generated first (see below). The
code stops with an error if `eos_table_file` is not set or cannot be read.

Restrictions (checked at startup):
- only `fluid = glmmhd` with `riemann = hlld` (the other Riemann solvers assume an ideal
  gas), and no first order flux correction (it uses the LLF solver)
- no `Tfloor`/`Tceil`, thermal conduction or tabular cooling, which all assume an ideal
  gas. The density and pressure floors work. When the pressure floor is applied, the
  stored pressure is the one of the internal energy that is set, so primitive and
  conserved variables stay consistent.
- the `collapse_be` problem generator sets the thermal energy with its own barotropic
  closure and therefore requires `eos = adiabatic`.

Problem generators set the initial state through the conserved variables. If they compute
the internal energy as `p / (gamma - 1)`, the actual initial pressure is the tabulated
pressure of that energy.

Outside the table the lookups are clamped to the nearest edge (no extrapolation). Every
clamped lookup is counted and, unless `eos_table_clamp_report = false`, a warning with the
number of lookups below or above each axis is printed once per cycle. Note that the
conversion to primitive variables also runs over ghost zones, so short excursions of
coarse/fine ghost cells below the table can show up here.

### Generating the table

The generator only needs numpy.

```bash
# coarse table, 180 x 220 (rho, esp) nodes, 1.2 MB, ~30 s
python3 src/eos/gen_eos_table.py --out eos_table.bin
# high resolution table used for protostellar collapse runs, 17.8 MB, ~7 min
src/eos/build_hires_v5_cold.sh eos_table_hires_v5_cold.bin
```

| table | (rho, esp) nodes | rho [g/cm^3] | T [K] | esp [erg/g] | md5 |
|---|---|---|---|---|---|
| default | 180 x 220 | 1e-20 .. 1 | 8 .. 3e5 | 4.2e8 .. 7.7e13 | `ff78c284c4f8361b27c122ece302b164` |
| `build_hires_v5_cold.sh` | 400 x 1399 | 1e-20 .. 1 | 0.046 .. 3e5 | 3.4e6 .. 7.7e13 | `f712d6a806cb82889214fec95a6c2681` |

The md5 sums are for numpy 2.4 on x86_64. The generator is deterministic, so a rebuild with
the same generator reproduces the file bit for bit.

The high resolution table is the 400 x 1000 x 920 table with 399 esp nodes and 450
temperature nodes added below the axes at the same spacing (`--esp-pad`, `--T-pad`), so
that gas cooled to a few K is still inside the table while all nodes of the unpadded table
keep their values. In the same way, `--T-pad-hi` extends the temperature axis above
`--T-max` (default 3e5 K). Changing `--nT`, `--T-max` or the density range instead moves
every node.

At high density the esp axis extends beyond the energy reached at the top of the
temperature axis (the grid is rectangular). There, pressure and temperature are flat in
esp, which the generator reports. This only matters for gas approaching 3e5 K.

### Tests

`src/eos/tests/test_eos_consistency.py [table.bin]` checks charge neutrality and nucleus
conservation of the Saha solution, the thermodynamic identity
`-rho^2 (de/drho)_T = T (dP/dT)_rho - P`, and (with a table) the tabulated pressure and
sound speed against the Saha EOS for 1e-16 to 1e-4 g/cm^3 and ~16 to 2e4 K. For the
high resolution table the median relative error is 1.6e-3 in pressure and 1.6e-3 in
`cs^2`; the largest errors (5.5e-3 in pressure, 3.3e-2 in `cs^2`) sit at the H
ionization kink. For the default table the medians are 6.8e-3 and 8.5e-3 and the maxima
6.7e-2 and 0.18.
