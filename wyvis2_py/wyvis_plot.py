"""Plot helpers mirroring the MATLAB wyvis_plot_*.m functions."""
import numpy as np
import matplotlib.pyplot as plt


def _stack(ud_list, field):
    cols = [np.asarray(getattr(u, field), dtype=float).ravel() for u in ud_list]
    return np.column_stack(cols)


def _as_list(ud):
    if isinstance(ud, (list, tuple)):
        return list(ud)
    return [ud]


# ---------- channel-only summary (example.py) ----------

def wyvis_plot_all_c(fig, ad, ud, savepath=None):
    ud_use = [_as_list(ud)[-1]]
    S = _stack(ud_use, 'S')
    phi_c = _stack(ud_use, 'phi_c')
    Q = _stack(ud_use, 'Q')
    M = _stack(ud_use, 'M')

    f = plt.figure(fig, figsize=(10/2.54, 20/2.54)); f.clf()
    fontsize = 12

    ax1 = f.add_axes([.12, .84, .82, .13])
    ax1.plot(ad.x/1e3, ad.phi_b/1e6, color=[0.8]*3, lw=1)
    ax1.plot(ad.x/1e3, ad.phi_s/1e6, color=[0.8]*3, lw=1)
    ax1.plot(ad.x/1e3, phi_c/1e6, 'b', lw=1)
    ax1.set_title('Hydraulic potential [ MPa ]', fontsize=fontsize)

    ax2 = f.add_axes([.12, .64, .82, .13])
    ax2.plot(ad.x/1e3, Q, 'b', lw=1)
    ax2.set_title(r'Discharge [ m$^3$/s ]', fontsize=fontsize)

    ax3 = f.add_axes([.12, .44, .82, .13])
    ax3.plot(ad.x/1e3, S, 'b', lw=1)
    ax3.set_title(r'Cross-section [ m$^2$ ]', fontsize=fontsize)

    ax4 = f.add_axes([.12, .24, .82, .13])
    ax4.plot(ad.x/1e3, M, 'b', lw=1)
    ax4.set_xlabel('Distance [ km ]', fontsize=fontsize)
    ax4.set_title(r'Melting rate [ m$^2$/s ]', fontsize=fontsize)

    if savepath:
        f.savefig(savepath, dpi=150, bbox_inches='tight')
    return f


# ---------- full snapshot (channel + sheet + ice vel) ----------

def wyvis_plot_all(fig, ad, ud, savepath=None):
    """Two-column plot of all variables at the last timestep."""
    last = [_as_list(ud)[-1]]
    S      = _stack(last, 'S')
    h      = _stack(last, 'h')
    phi_c  = _stack(last, 'phi_c')
    phi    = _stack(last, 'phi')
    u      = _stack(last, 'u')
    Q      = _stack(last, 'Q')
    q      = _stack(last, 'q')
    M      = _stack(last, 'M')
    m      = _stack(last, 'm')
    kappa  = _stack(last, 'kappa')

    f = plt.figure(fig, figsize=(20/2.54, 20/2.54)); f.clf()
    fontsize = 12
    x_km = ad.x/1e3
    W_col = ad.W.reshape(-1, 1)

    def add(pos, title):
        ax = f.add_axes(pos); ax.set_title(title, fontsize=fontsize); return ax

    ax = add([.06, .84, .43, .13], 'Hydraulic potential [ MPa ]')
    ax.plot(x_km, ad.phi_b/1e6, color=[0.8]*3, lw=1)
    ax.plot(x_km, ad.phi_s/1e6, color=[0.8]*3, lw=1)
    ax.plot(x_km, phi_c/1e6, 'b', lw=1)
    ax.plot(x_km, phi/1e6, 'r', lw=1)

    ax = add([.06, .64, .43, .13], r'Discharge [ m$^3$/s ]')
    ax.plot(x_km, Q, 'b', lw=1)
    ax.plot(x_km, Q + W_col*q, 'k', lw=1)

    ax = add([.06, .44, .43, .13], r'Cross-section [ m$^2$ ]')
    ax.plot(x_km, S, 'b', lw=1)

    ax = add([.06, .24, .43, .13], r'Melting rate [ m$^2$/s ]')
    ax.plot(x_km, M, 'b', lw=1)

    ax = add([.06, .04, .43, .13], 'Elevation [ m ]')
    ax.plot(x_km, ad.Z_b, 'k', lw=1)
    ax.plot(x_km, ad.Z_s, 'k', lw=1)
    ax.set_xlabel('Distance [ km ]', fontsize=fontsize)

    ax = add([.56, .84, .43, .13], r'Exchange [ m$^2$/s ]')
    ax.plot(x_km, 0*ad.x, color=[0.8]*3, lw=1)
    ax.plot(x_km, kappa, 'k', lw=1)

    ax = add([.56, .64, .43, .13], r'Discharge [ m$^2$/s ]')
    ax.plot(x_km, q, 'r', lw=1)

    ax = add([.56, .44, .43, .13], 'Sheet depth [ m ]')
    ax.plot(x_km, h, 'r', lw=1)

    ax = add([.56, .24, .43, .13], 'Melting rate [ mm/d ]')
    ax.plot(x_km, m*24*60*60*1e3, 'r', lw=1)

    ax = add([.56, .04, .43, .13], 'Ice velocity [ m/y ]')
    ax.plot(x_km, u*365*24*60*60, 'k', lw=1)
    ax.set_xlabel('Distance [ km ]', fontsize=fontsize)

    if savepath:
        f.savefig(savepath, dpi=150, bbox_inches='tight')
    return f


# ---------- time evolution of all variables ----------

def wyvis_plot_allt(fig, td, ud, xi=None, savepath=None):
    """Multi-panel time series of all variables at grid index xi (0-based)."""
    ud = _as_list(ud)
    S     = _stack(ud, 'S')
    h     = _stack(ud, 'h')
    phi_c = _stack(ud, 'phi_c')
    phi   = _stack(ud, 'phi')
    u     = _stack(ud, 'u')
    Q     = _stack(ud, 'Q')
    q     = _stack(ud, 'q')
    N_c   = _stack(ud, 'N_c')
    N     = _stack(ud, 'N')
    M     = _stack(ud, 'M')
    m     = _stack(ud, 'm')
    kappa = _stack(ud, 'kappa')

    if xi is None:
        xi = np.arange(phi_c.shape[0])
    xi = np.atleast_1d(xi)
    td_d = td/86400.0

    f = plt.figure(fig, figsize=(20/2.54, 20/2.54)); f.clf()
    fontsize = 12

    def add(pos, title):
        ax = f.add_axes(pos); ax.set_title(title, fontsize=fontsize); return ax

    ax = add([.06, .84, .43, .13], 'Effective pressure [ MPa ]')
    ax.plot(td_d, N_c[xi,:].T/1e6, 'b', lw=1)
    ax.plot(td_d, N[xi,:].T/1e6,   'r', lw=1)

    ax = add([.06, .64, .43, .13], r'Discharge [ m$^3$/s ]')
    ax.plot(td_d, Q[xi,:].T, 'b', lw=1)

    ax = add([.06, .44, .43, .13], r'Cross-section [ m$^2$ ]')
    ax.plot(td_d, S[xi,:].T, 'b', lw=1)

    ax = add([.06, .24, .43, .13], r'Melting rate [ m$^2$/s ]')
    ax.plot(td_d, M[xi,:].T, 'b', lw=1)
    ax.set_xlabel('Time [ d ]', fontsize=fontsize)

    ax = add([.56, .84, .43, .13], r'Exchange [ m$^2$/s ]')
    ax.plot(td_d, 0*td_d, color=[0.8]*3, lw=1)
    ax.plot(td_d, kappa[xi,:].T, 'k', lw=1)

    ax = add([.56, .64, .43, .13], r'Discharge [ m$^2$/s ]')
    ax.plot(td_d, q[xi,:].T, 'r', lw=1)

    ax = add([.56, .44, .43, .13], 'Sheet depth [ m ]')
    ax.plot(td_d, h[xi,:].T, 'r', lw=1)

    ax = add([.56, .24, .43, .13], 'Melting rate [ mm/d ]')
    ax.plot(td_d, m[xi,:].T*24*60*60*1e3, 'r', lw=1)

    ax = add([.56, .04, .43, .13], 'Ice velocity [ m/y ]')
    ax.plot(td_d, (365*24*60*60)*u[xi,:].T, 'k', lw=1)
    ax.set_xlabel('Time [ d ]', fontsize=fontsize)

    if savepath:
        f.savefig(savepath, dpi=150, bbox_inches='tight')
    return f


# ---------- channel-potential time series ----------

def wyvis_plot_phi_ct(fig, td, ud, xi=None, savepath=None):
    phi_c = _stack(_as_list(ud), 'phi_c')
    if xi is None:
        xi = np.arange(phi_c.shape[0])
    xi = np.atleast_1d(xi)
    f = plt.figure(fig, figsize=(20/2.54, 20/2.54)); f.clf()
    ax = f.add_axes([.12, .15, .82, .75])
    ax.plot(td/86400.0, phi_c[xi, :].T/1e6, 'b', lw=2)
    ax.set_xlabel('Time [ d ]', fontsize=12)
    ax.set_title('Hydraulic potential [ MPa ]', fontsize=12)
    if savepath:
        f.savefig(savepath, dpi=150, bbox_inches='tight')
    return f


# ---------- bed/surface geometry ----------

def wyvis_plot_geometry(fig, ad, savepath=None):
    f = plt.figure(fig, figsize=(12/2.54, 6/2.54)); f.clf()
    ax = f.add_axes([.15, .2, .8, .7])
    ax.plot(ad.x/1e3, ad.Z_b, 'k', lw=2)
    ax.plot(ad.x/1e3, ad.Z_s, 'k', lw=2)
    x_m = getattr(ad, 'x_m', None)
    if x_m is not None and np.asarray(x_m).size > 0:
        from wyvis_discretize import wyvis_nearest_gridpoint
        xi = wyvis_nearest_gridpoint(x_m, ad.x)
        ax.plot(ad.x[xi]/1e3, ad.Z_b[xi], 'bo', ms=8,
                markerfacecolor=[0.6]*3)
    ax.set_xlabel('Distance [ km ]', fontsize=12)
    ax.set_ylabel('Elevation [ m ]', fontsize=12)
    if savepath:
        f.savefig(savepath, dpi=150, bbox_inches='tight')
    return f


# ---------- annual signal (forcing helper) ----------

def wyvis_annual_signal(t, t_spr, t_aut, t_per):
    """Two-tanh annual cycle, centered on t_spr/t_aut with width t_per."""
    ty = 365*24*60*60
    return (np.tanh((np.mod(t, ty) - t_spr)/t_per)
            - np.tanh((np.mod(t, ty) - t_aut)/t_per)) / 2.0
