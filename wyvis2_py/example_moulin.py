"""Python port of example_moulin.m.

Flat bed and linear surface slope, conduit drainage system.
Phase 1: 1 year toward steady state with distributed input only.
Phase 2: add a moulin at x=0 with ramping flux, continue another year.
Saves PNG plots and the numerical results (.npz) next to this script.
"""
from types import SimpleNamespace
from pathlib import Path
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from wyvis_parameters import wyvis_default_parameters
from wyvis_solve import wyvis_solve
from wyvis_plot import wyvis_plot_all_c

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

    print('--- Phase 1: 1 year to steady state (no moulin) ---')
    t0 = time.time()
    td1, ad1, ud1 = wyvis_solve(
        np.arange(0, 366) * 24*60*60, ad, pd)
    print(f'phase 1 wall time: {time.time()-t0:.2f}s')
    wyvis_plot_all_c(1, ad1, ud1[-1], savepath=OUT/'python_fig_moulin_steady.png')

    print('--- Phase 2: add moulin, 1 year continuation ---')
    ad1.x_m = np.array([0.0])
    ad1.S_m = np.array([100.0])
    ty = pd.ty
    ad1.Q_m = lambda t: np.array([10.0*(1.0 - np.exp(-t/(7.0/365.0*ty)))])
    t0 = time.time()
    td2, ad2, ud2 = wyvis_solve(
        np.arange(0, 366) * 24*60*60, ad1, pd, ud=[ud1[-1]])
    print(f'phase 2 wall time: {time.time()-t0:.2f}s')
    wyvis_plot_all_c(1, ad2, ud2[-1], savepath=OUT/'python_fig_moulin_final.png')

    Q_t = np.column_stack([np.asarray(u.Q).ravel() for u in ud2])
    f3 = plt.figure(3, figsize=(20/2.54, 12/2.54)); f3.clf()
    ax = f3.add_axes([.12, .15, .78, .75])
    cs = ax.contourf(ad2.x/1e3, td2/pd.td, Q_t.T, levels=20)
    f3.colorbar(cs, ax=ax)
    ax.set_xlabel('Distance [km]')
    ax.set_ylabel('Time [d]')
    f3.savefig(OUT/'python_fig_moulin_Qt.png', dpi=150, bbox_inches='tight')

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
    out['phase2_x_m'] = np.asarray(ad2.x_m).ravel()
    out['phase2_S_m'] = np.asarray(ad2.S_m).ravel()
    out['phase2_xim'] = np.asarray(ad2.xim).ravel()
    np.savez(OUT/'python_example_moulin_results.npz', **out)
    print(f'saved {OUT/"python_example_moulin_results.npz"}')


if __name__ == '__main__':
    main()
