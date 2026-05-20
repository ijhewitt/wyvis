"""Python port of example2.m.

Flat bed and linear surface slope, driven by a steady distributed water input,
coupled drainage system with the ice-velocity equation enabled. Solves for one
year toward steady state.

Note: MATLAB ``pd.lambda`` is named ``pd.lambda_`` here (``lambda`` is a Python
reserved word).
"""
from types import SimpleNamespace
from pathlib import Path
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')

from wyvis_parameters import wyvis_default_parameters
from wyvis_solve import wyvis_solve
from wyvis_plot import wyvis_plot_all

HERE = Path(__file__).resolve().parent
OUT = HERE/'outputs'; OUT.mkdir(exist_ok=True)


def main():
    pd = wyvis_default_parameters()
    pd.lambda_ = 1e-9                      # conduit-sheet exchange [m^2/s/Pa]
    print('pd =', pd.__dict__)

    ad = SimpleNamespace()
    ad.x = np.linspace(0, 10e3, 40)
    b0 = 5*100
    ad.Z_b = 0*ad.x - b0                    # constant negative bed
    ad.Z_s = 0.01*(ad.x[-1] - ad.x) + (pd.rho_w/pd.rho_i - 1)*b0
    ad.M_in = lambda t: np.ones_like(ad.x) * 1e3   # distributed conduit input [m^2/s]

    oo = SimpleNamespace(include_ice_velocity=1)

    print('--- Solving 1 year to steady state (ice velocity coupled) ---')
    t0 = time.time()
    td, ad_out, ud = wyvis_solve(
        np.arange(0, 366) * 24*60*60, ad, pd, oo=oo)
    print(f'wall time: {time.time()-t0:.2f}s')

    wyvis_plot_all(1, ad_out, ud[-1], savepath=OUT/'python_ex2_fig_all.png')

    # Save numerical results for comparison with MATLAB reference
    def stack(ud, name):
        return np.column_stack([getattr(u, name).ravel() for u in ud])

    fields = ['S','phi_c','h','phi','u','Q','N_c','M','q','N','m','kappa']
    out = {
        'td': td,
        'x': ad_out.x,
        'Z_b': ad_out.Z_b, 'Z_s': ad_out.Z_s,
        'phi_b': ad_out.phi_b, 'phi_s': ad_out.phi_s,
        'W': ad_out.W,
    }
    for f in fields:
        out[f] = stack(ud, f)
    np.savez(OUT/'python_example2_results.npz', **out)
    print(f'saved {OUT/"python_example2_results.npz"}')


if __name__ == '__main__':
    main()
