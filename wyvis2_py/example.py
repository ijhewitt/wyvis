"""Python port of example.m.

Runs the wyvis subglacial-hydrology model for one year toward steady state,
then 10 days with a diurnal water-input forcing. Saves PNG plots and the
numerical results (.npz) next to this script.
"""
from types import SimpleNamespace
from pathlib import Path
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')

from wyvis_parameters import wyvis_default_parameters
from wyvis_solve import wyvis_solve
from wyvis_plot import wyvis_plot_all_c, wyvis_plot_phi_ct

HERE = Path(__file__).resolve().parent
OUT = HERE/'outputs'; OUT.mkdir(exist_ok=True)


def main():
    pd = wyvis_default_parameters()
    print('pd =', pd.__dict__)

    ad = SimpleNamespace()
    ad.x = np.linspace(0, 10000, 40)
    ad.Z_b = 0*ad.x
    ad.Z_s = 0.1*(ad.x[-1] - ad.x)
    ad.M_in = lambda t: np.ones_like(ad.x) * 1e-3

    print('--- Phase 1: 1 year to steady state ---')
    t0 = time.time()
    td1, ad1, ud1 = wyvis_solve(
        np.arange(0, 366) * 24*60*60, ad, pd)
    print(f'phase 1 wall time: {time.time()-t0:.2f}s')
    wyvis_plot_all_c(1, ad1, ud1[-1], savepath=OUT/'python_fig_steady.png')

    print('--- Phase 2: 10 days with diurnal forcing ---')
    ad1.M_in = lambda t: np.ones_like(ad1.x) * 1e-3 * (1 + 0.5*np.sin(2*np.pi*t/24/60/60))
    t0 = time.time()
    td2, ad2, ud2 = wyvis_solve(
        np.arange(0, 241)/24 * 24*60*60, ad1, pd, ud=[ud1[-1]])
    print(f'phase 2 wall time: {time.time()-t0:.2f}s')
    wyvis_plot_all_c(1, ad2, ud2[-1], savepath=OUT/'python_fig_diurnal_final.png')
    wyvis_plot_phi_ct(1, td2, ud2, xi=19, savepath=OUT/'python_fig_phi_ct.png')

    def stack(ud, name):
        return np.column_stack([getattr(u, name).ravel() for u in ud])

    fields = ['S','phi_c','h','phi','u','Q','N_c','M','q','N','m','kappa']
    out = {}
    for label, td, ad_, ud_ in [('phase1', td1, ad1, ud1), ('phase2', td2, ad2, ud2)]:
        out[f'{label}_td'] = td
        out[f'{label}_x'] = ad_.x
        out[f'{label}_Z_b'] = ad_.Z_b
        out[f'{label}_Z_s'] = ad_.Z_s
        out[f'{label}_phi_b'] = ad_.phi_b
        out[f'{label}_phi_s'] = ad_.phi_s
        for f in fields:
            out[f'{label}_{f}'] = stack(ud_, f)
    np.savez(OUT/'python_example_results.npz', **out)
    print(f'saved {OUT/"python_example_results.npz"}')


if __name__ == '__main__':
    main()
