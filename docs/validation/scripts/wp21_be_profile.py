#!/usr/bin/env python
"""WP-21: error of the BEProfile analytic approximation in collapse_be.cpp
against a directly integrated isothermal Lane-Emden (Bonnor-Ebert) solution.

Code approximation (collapse_be.cpp:66-68, "Tomida 2011"):
    rho(r)/rho_c = (1 + r^2 / rcsq)^(-3/2),   rcsq = 26/3
truncated at rc_code = 6.45 (collapse_be.cpp:40).

Truth: isothermal Lane-Emden
    psi'' + (2/xi) psi' = exp(-psi),  psi(0)=0, psi'(0)=0
    rho/rho_c = exp(-psi)
In units where the code length unit is the isothermal scale a = c_s/sqrt(4 pi G rho_c),
the code radius r IS the dimensionless xi.
"""
import numpy as np
from scipy.integrate import solve_ivp, quad

RCSQ = 26.0 / 3.0
RC = 6.45
BEMASS_CODE = 197.561  # collapse_be.cpp:42


def rhs(xi, y):
    psi, dpsi = y
    # series-expand near origin to avoid 2/xi singularity
    if xi < 1e-8:
        return [dpsi, np.exp(-psi) / 3.0]
    return [dpsi, np.exp(-psi) - 2.0 * dpsi / xi]


# integrate well past the critical radius
sol = solve_ivp(rhs, [0.0, 20.0], [0.0, 0.0], rtol=1e-12, atol=1e-14,
                dense_output=True, method="DOP853")
assert sol.success, sol.message


def rho_exact(xi):
    xi = np.atleast_1d(xi)
    out = np.empty_like(xi, dtype=float)
    small = xi < 1e-8
    out[small] = 1.0
    if np.any(~small):
        out[~small] = np.exp(-sol.sol(xi[~small])[0])
    return out


def rho_approx(xi):
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    return (1.0 + xi * xi / RCSQ) ** -1.5


# ---------------------------------------------------------------- critical xi
# The critical (marginally stable) BE sphere is at the first maximum of
# the dimensionless mass  m(xi) = xi^2 psi'(xi) * exp(-psi/2) ... use the
# standard criterion: max of  xi^2 |psi'| / exp(psi/2)  -> xi_crit ~ 6.451
xg = np.linspace(0.5, 12.0, 400001)
psi_g = sol.sol(xg)[0]
dpsi_g = sol.sol(xg)[1]
mfunc = xg ** 2 * dpsi_g * np.exp(-0.5 * psi_g)
xi_crit = xg[np.argmax(mfunc)]

print("=" * 74)
print("WP-21  BE profile: analytic approximation vs integrated Lane-Emden")
print("=" * 74)
print(f"critical xi (integrated)      : {xi_crit:.4f}   [textbook 6.451]")
print(f"rc_code used by pgen          : {RC}")
print(f"  => truncation radius error  : {100*(RC-xi_crit)/xi_crit:+.3f} %")

# ------------------------------------------------------------- profile error
xi = np.linspace(0.0, RC, 200001)
re = rho_exact(xi)
ra = rho_approx(xi)
rel = (ra - re) / re

print()
print("--- density profile error over 0 <= xi <= 6.45 (rho_approx/rho_exact - 1)")
print(f"  max |rel err|               : {np.max(np.abs(rel))*100:.3f} %  at xi = {xi[np.argmax(np.abs(rel))]:.3f}")
print(f"  RMS rel err                 : {np.sqrt(np.mean(rel**2))*100:.3f} %")
print(f"  rel err at xi=0             : {rel[0]*100:+.4f} %")
print(f"  rel err at edge xi=6.45     : {rel[-1]*100:+.3f} %")
print(f"  rho_exact(6.45)/rho_c       : {re[-1]:.6f}")
print(f"  rho_approx(6.45)/rho_c      : {ra[-1]:.6f}")
print(f"  centre-to-edge contrast exact : {1.0/re[-1]:.3f}")
print(f"  centre-to-edge contrast approx: {1.0/ra[-1]:.3f}   [textbook 14.04]")

print()
print("  xi      rho_exact    rho_approx    rel.err[%]")
for x in [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 6.45]:
    e = rho_exact(x)[0]
    a = rho_approx(x)[0]
    print(f"  {x:5.2f}  {e:11.6f}  {a:11.6f}  {100*(a-e)/e:+9.3f}")

# ----------------------------------------------------------------- mass error
# enclosed mass in code units:  M = int_0^R 4 pi xi^2 rho dxi   (rho_c = 1)
m_exact = quad(lambda x: 4*np.pi*x*x*rho_exact(x)[0], 0, RC, limit=400)[0]
m_approx = quad(lambda x: 4*np.pi*x*x*rho_approx(x)[0], 0, RC, limit=400)[0]
# exact mass via the first integral: M = 4 pi xi^2 psi'(xi)
m_exact_fi = 4*np.pi*RC**2*sol.sol(RC)[1]

print()
print("--- enclosed mass at xi = 6.45, code units (rho_c = 1, four_pi_G = 1)")
print(f"  M exact  (quadrature)       : {m_exact:.4f}")
print(f"  M exact  (first integral)   : {m_exact_fi:.4f}   [consistency check]")
print(f"  M approx (quadrature)       : {m_approx:.4f}")
print(f"  mass error of approximation : {100*(m_approx-m_exact)/m_exact:+.3f} %")
print(f"  bemass_code in pgen         : {BEMASS_CODE}")
print(f"    ratio bemass_code/M_approx: {BEMASS_CODE/m_approx:.4f}")
print(f"    ratio bemass_code/M_exact : {BEMASS_CODE/m_exact:.4f}")

# --------------------------------------------- free-fall time implication
# t_ff ~ 1/sqrt(G rho); a fractional density error d gives d(t_ff)/t_ff = -d/2
mean_rho_exact = m_exact / (4.0/3.0*np.pi*RC**3)
mean_rho_approx = m_approx / (4.0/3.0*np.pi*RC**3)
d = (mean_rho_approx - mean_rho_exact)/mean_rho_exact
print()
print("--- consequence for the collapse clock")
print(f"  mean density exact          : {mean_rho_exact:.6f}")
print(f"  mean density approx         : {mean_rho_approx:.6f}")
print(f"  mean-density error          : {100*d:+.3f} %")
print(f"  => t_ff shift               : {-50*d:+.3f} %  (t_ff ~ rho^-1/2)")
print("=" * 74)
