#!/usr/bin/env python
"""Assemble the final artifact HTML: substitute figures (base64) and extracted numbers."""
import base64, json, numpy as np

OUT = '/home/bbg6470/.claude/jobs/873b4541/tmp'
T0YR = 1.4789e12/3.1557e7
RHO0 = 5.467e-19

D = np.load(f'{OUT}/t4_data.npz')
key = json.load(open(f'{OUT}/t4_key_numbers.json'))

# e-fold time from the last two tail samples
ts, rho = D['ts'], D['rhomax']
efold_code = (ts[-1]-ts[-3]) / np.log(rho[-1]/rho[-3])
efold_yr = efold_code * T0YR
print(f'e-fold (samples -3..-1): {efold_yr:.3f} yr')

# plateau duration: cycles 15000..67500 samples
i15 = int(np.argmin(np.abs(D['cyc']-15000))); i67 = int(np.argmin(np.abs(D['cyc']-67500)))
plateau_yr = (ts[i67]-ts[i15]) * T0YR
print(f'plateau {D["cyc"][i15]}..{D["cyc"][i67]}: {plateau_yr:.1f} yr, '
      f'rho {rho[i15]*RHO0:.2e}..{rho[i67]*RHO0:.2e}')

Tc = key['central']['T_est_K']
xec = key['central']['xe']
print(f'central T_est {Tc:.0f} K, xe {xec:.2e}')

html = open(f'{OUT}/star_formation_flow.template.html').read()

def b64(p):
    return base64.b64encode(open(p, 'rb').read()).decode()

subs = {
    '{{F1}}': b64(f'{OUT}/f1_rhomax.png'),
    '{{F2}}': b64(f'{OUT}/f2_slices.png'),
    '{{F3}}': b64(f'{OUT}/f3_profiles.png'),
    '{{F4}}': b64(f'{OUT}/f4_brho.png'),
    '{{PLATEAU_TEXT}}': f'~{plateau_yr:.0f} yr quasi-static',
    '{{PLATEAU_YEARS}}': f'{plateau_yr:.0f}',
    '{{T_CENTRAL_NOTE}}': (f'The raw central temperature readout at the stop is ~{Tc:,.0f} K, but this '
                           f'is contaminated low by the STS/radiation bug (caveat below); the exact-EOS '
                           f'isentrope from the last clean state gives ~1500–2200 K. The central '
                           f'ionization fraction sits at the 10⁻¹⁵ floor — deeply decoupled — though a '
                           f'hotter corrected core may partially re-couple via thermal (K) ionization: '
                           f'a first-order question for the fossil-flux budget, answered by the re-run.'),
    '{{EFOLD_TEXT}}': f'{efold_yr:.2f} yr (measured)',
    '{{EFOLD_NOTE}}': (f'Density e-folding time measured across the final snapshots: '
                       f'{efold_yr:.2f} yr — a genuine dynamical runaway.'),
}
for k, v in subs.items():
    assert k in html, f'missing token {k}'
    html = html.replace(k, v)
open(f'{OUT}/star_formation_flow.html', 'w').write(html)
import os
print('bytes:', os.path.getsize(f'{OUT}/star_formation_flow.html'))
print('ASSEMBLED')
