"""Python port of wyvis_solve.m — 1D subglacial hydrology / ice-velocity solver.

The MATLAB code expresses the system as M(Y) * Y' = F(t, Y) and integrates
with ode15s and a sparse Jacobian pattern. SciPy's solve_ivp has no mass-
matrix API, so we exploit the block structure of M (two 2x2 blocks plus
three diagonals) and solve M(Y) * Y' = F in closed form at each evaluation.

Default options match wyvis_solve.m: continuing=detected from ud input,
include_a4=1, include_a44=1, ice velocity off, etc.
"""
from types import SimpleNamespace
import numpy as np
import scipy.sparse as sp
from scipy.integrate import solve_ivp

from wyvis_discretize import wyvis_discretize, wyvis_nearest_gridpoint
from wyvis_parameters import wyvis_default_parameters


# ---------- scaling + discretization ----------

def _scale_and_discretize(ad, pd, oo):
    """Build the dimensionless state: dd, aa, pp, ps."""
    # Defaults for optional ad fields
    if not hasattr(ad, 'Z_b') and hasattr(ad, 'b'):
        ad.Z_b = ad.b
    if not hasattr(ad, 'Z_s') and hasattr(ad, 's'):
        ad.Z_s = ad.s
    x = np.asarray(ad.x, dtype=float).reshape(-1)
    if not hasattr(ad, 'W'):
        ad.W = pd.W * np.ones_like(x)
    if not hasattr(ad, 'U_b'):
        ad.U_b = pd.U_b * np.ones_like(x)
    if not hasattr(ad, 'M_in'):
        ad.M_in = lambda t: 0*x
    if not hasattr(ad, 'm_in'):
        ad.m_in = lambda t: 0*x
    if not hasattr(ad, 'x_m'):
        ad.x_m = np.array([], dtype=float)
    if not hasattr(ad, 'Q_m'):
        ad.Q_m = lambda t: 0*np.asarray(ad.x_m, dtype=float)
    if not hasattr(ad, 'S_m'):
        ad.S_m = 0*np.asarray(ad.x_m, dtype=float)

    ad.x = x
    ad.Z_b = np.asarray(ad.Z_b, dtype=float).reshape(-1)
    ad.Z_s = np.asarray(ad.Z_s, dtype=float).reshape(-1)
    ad.W   = np.asarray(ad.W, dtype=float).reshape(-1)
    ad.U_b = np.asarray(ad.U_b, dtype=float).reshape(-1)
    ad.x_m = np.asarray(ad.x_m, dtype=float).reshape(-1)
    ad.S_m = np.asarray(ad.S_m, dtype=float).reshape(-1)

    ps = SimpleNamespace()
    ps.x = 10*1e3
    ps.s = 1000
    ps.u = pd.U_b
    ps.Q = 1
    ps.W = 1e3
    ps.tau = pd.rho_i*pd.g*ps.s**2/ps.x
    ps.phi = pd.rho_i*pd.g*ps.s
    ps.S = (ps.Q/pd.K_c/(ps.phi/ps.x)**0.5)**(1.0/pd.alpha_c)
    ps.M = ps.Q*ps.phi/ps.x/pd.L/pd.rho_w
    ps.t = pd.rho_i*ps.S/ps.M/pd.rho_w
    ps.kappa = ps.Q/ps.x
    ps.q = ps.Q/ps.W
    ps.h = (ps.q/pd.K/(ps.phi/ps.x))**(1.0/pd.alpha)
    ps.N = ps.phi
    ps.m = ps.q/ps.x
    ps.M_in = ps.Q/ps.x
    ps.m_in = ps.q/ps.x
    ps.sigma = pd.sigma
    ps.sigma_c = pd.sigma_c
    ps.Sigma = ps.sigma*ps.phi/pd.rho_w/pd.g
    ps.Sigma_c = ps.sigma_c*ps.phi/pd.rho_w/pd.g
    ps.V_in = ps.S*ps.phi/pd.rho_w/pd.g

    # Dimensionless parameters
    pp = SimpleNamespace()
    pp.alpha_c = pd.alpha_c
    pp.alpha = pd.alpha
    pp.n = pd.n
    pp.p_slide = pd.p_slide
    pp.q_slide = pd.q_slide
    pp.a1  = 1
    pp.a2  = pd.rho_w*ps.M*ps.t/pd.rho_i/ps.S
    pp.a3  = pd.Atil*ps.N**pd.n*ps.t
    pp.a4  = ps.u*pd.h_rc*ps.t/ps.S
    pp.a5  = ps.S/pd.h_rc/pd.l_rc if pd.h_rc != 0 else 0.0  # h_rc=0 by default; matches MATLAB's 1/0 -> inf later guarded by max(1-a5*S,0)
    pp.a6  = ps.Q*ps.phi/ps.x/pd.L/ps.M/pd.rho_w
    pp.a7  = ps.Q*ps.phi/ps.x/pd.L/ps.M/pd.rho_w*pd.gamma*pd.rho_w*pd.c
    pp.a8  = ps.m*pd.l_rc/ps.M
    pp.a9  = pd.G/ps.m/pd.rho_w/pd.L
    pp.a10 = ps.tau*ps.u/ps.m/pd.rho_w/pd.L
    pp.a11 = ps.q*ps.phi/ps.x/ps.m/pd.rho_w/pd.L
    pp.a12 = ps.q*ps.phi/ps.x/ps.m/pd.rho_w/pd.L*pd.gamma*pd.rho_w*pd.c
    pp.a13 = 1
    pp.a14 = ps.u*pd.h_r/pd.l_r*ps.t/ps.h
    pp.a15 = ps.h/pd.h_r
    pp.a16 = pd.Ahat*ps.N**pd.n*ps.t
    pp.a17 = ps.S*ps.x/ps.t/ps.Q
    pp.a18 = 1
    pp.a19 = ps.M*ps.x/ps.Q
    pp.a20 = 1
    pp.a21 = ps.kappa*ps.x/ps.Q
    pp.a22 = ps.h*ps.x/ps.t/ps.q
    pp.a23 = ps.x/ps.t/ps.q*ps.sigma*ps.phi/pd.rho_w/pd.g
    pp.a24 = 1
    pp.a25 = ps.m*ps.x/ps.q
    pp.a26 = ps.kappa*ps.x/ps.q/ps.W
    pp.a27 = pd.lambda_*ps.N/ps.kappa
    pp.a28 = ps.Q/pd.K_c/ps.S**pd.alpha_c/(ps.phi/ps.x)**0.5
    pp.a29 = pd.K*ps.h**pd.alpha*(ps.phi/ps.x)/ps.q
    pp.a30 = pd.rho_w*ps.m*ps.t/pd.rho_i/ps.h
    pp.a31 = ps.x/ps.t/ps.Q*pd.sigma_c*ps.W*ps.phi/pd.rho_w/pd.g
    pp.a32 = pd.rho_w/pd.rho_i
    pp.a33 = ps.S/ps.t/ps.Q*ps.phi/pd.rho_w/pd.g
    pp.a34 = ps.M_in*ps.x/ps.Q
    pp.a35 = ps.m_in*ps.x/ps.q
    pp.a36 = 1
    pp.a37 = pd.p_ff/ps.phi
    pp.a38 = ps.S*ps.x/ps.t/ps.Q*pd.beta_c*ps.phi
    pp.a39 = pd.lambda2*pd.K*ps.h**pd.alpha*ps.phi
    pp.a40 = pd.p_u/ps.phi
    pp.a41 = pd.E_u/ps.h*ps.phi
    pp.a42 = ps.h*ps.x/ps.t/ps.q
    pp.a43 = ps.h*ps.x/ps.t/ps.q*pd.beta*ps.phi
    pp.a44 = pd.Abar**(-1.0/pp.n)*(ps.u/ps.x)**(1.0/pp.n)*ps.s/ps.x/ps.tau
    pp.a45 = pd.rho_i*pd.g*ps.s**2/ps.x/ps.tau
    pp.a46 = 1
    pp.a47 = pd.Re_reg*ps.u/ps.t/ps.tau
    pp.a48 = pd.C_slide*ps.u**pp.q_slide*ps.N**pp.p_slide/ps.tau
    pp.a49 = pd.eps_reg/(ps.u/ps.x)
    pp.a50 = pd.N_min/ps.phi
    pp.a51 = pd.phi_m/ps.phi
    pp.a52 = (1 - pd.rho_i/pd.rho_o)
    pp.a53 = pd.Abar**(-1.0/pp.n)*(ps.u/ps.W)**(1.0/pp.n)*ps.s/ps.W/ps.tau

    if not oo.include_a4:  pp.a4 = 0
    if not oo.include_a8:  pp.a8 = 0
    if not oo.include_a11: pp.a11 = 0
    if not oo.include_a12: pp.a12 = 0
    if not oo.include_a44: pp.a44 = 0
    if not oo.include_a53: pp.a53 = 0

    # Discretization on dimensionless grid
    dd = wyvis_discretize(x/ps.x)

    # Dimensionless geometry / inputs
    aa = SimpleNamespace()
    aa.x = x/ps.x
    aa.Z_b = ad.Z_b/ps.s
    aa.Z_s = ad.Z_s/ps.s
    aa.W = ad.W/ps.W
    aa.U_b = ad.U_b/ps.u
    aa.x_m = ad.x_m/ps.x
    Q_m_fn = ad.Q_m
    M_in_fn = ad.M_in
    m_in_fn = ad.m_in
    aa.Q_m = lambda t: np.asarray(Q_m_fn(ps.t*t), dtype=float)/ps.Q
    aa.S_m = ad.S_m/ps.S
    aa.M_in = lambda t: np.asarray(M_in_fn(ps.t*t), dtype=float).reshape(-1)/ps.M_in
    aa.m_in = lambda t: np.asarray(m_in_fn(ps.t*t), dtype=float).reshape(-1)/ps.m_in

    aa.xim = wyvis_nearest_gridpoint(aa.x_m, dd.x)
    aa.phi_s = aa.Z_s + (pp.a32 - 1)*aa.Z_b
    aa.phi_b = pp.a32*aa.Z_b
    aa.phi_ext = np.maximum(aa.phi_b, pp.a51)
    aa.phi_ext[dd.xin] = aa.phi_s[dd.xin]
    aa.tau_d = -(aa.Z_s - aa.Z_b) * (dd.ddxx @ aa.Z_s)
    aa.H = aa.Z_s - aa.Z_b
    aa.sigma = np.ones(dd.I)
    aa.sigma_c = np.ones(dd.I)
    # MATLAB sets these to NaN so sparse matvec implicitly skips the boundary
    # column. With dense matvec we need 0 here; the NaN is restored on output.
    aa.h_ext = np.zeros(dd.I)
    aa.S_ext = np.zeros(dd.I)
    aa.Q_ext = np.zeros_like(aa.x)
    aa.q_ext = np.zeros_like(aa.x)
    aa.u_ext = aa.U_b.copy()

    # Unscaled derived (mirrored back into ad for output)
    ad.xim = aa.xim
    ad.phi_s = ps.phi*aa.phi_s
    ad.phi_b = ps.phi*aa.phi_b
    ad.phi_ext = ps.phi*aa.phi_ext
    ad.tau_d = ps.tau*aa.tau_d
    ad.H = ps.s*aa.H
    ad.U_b = ps.u*aa.U_b
    ad.sigma = ps.sigma*aa.sigma
    ad.sigma_c = ps.sigma_c*aa.sigma_c

    return dd, aa, pp, ps


def _unscale_aa(aa, ps):
    """Mirror of wyvis_unscale_aa.m: produce an ad-like namespace."""
    ad = SimpleNamespace()
    ad.x = ps.x*aa.x
    ad.Z_b = ps.s*aa.Z_b
    ad.Z_s = ps.s*aa.Z_s
    ad.W = ps.W*aa.W
    ad.U_b = ps.u*aa.U_b
    ad.x_m = ps.x*aa.x_m
    Q_m_aa = aa.Q_m
    M_in_aa = aa.M_in
    m_in_aa = aa.m_in
    ad.Q_m = lambda t: ps.Q*Q_m_aa(t/ps.t)
    ad.S_m = ps.S*aa.S_m
    ad.M_in = lambda t: ps.M_in*M_in_aa(t/ps.t)
    ad.m_in = lambda t: ps.m_in*m_in_aa(t/ps.t)
    ad.xim = aa.xim
    ad.phi_s = ps.phi*aa.phi_s
    ad.phi_b = ps.phi*aa.phi_b
    ad.phi_ext = ps.phi*aa.phi_ext
    ad.tau_d = ps.tau*aa.tau_d
    ad.H = ps.s*aa.H
    ad.sigma = ps.sigma*aa.sigma
    ad.sigma_c = ps.sigma_c*aa.sigma_c
    ad.Q_ext = ps.Q*aa.Q_ext
    ad.q_ext = ps.q*aa.q_ext
    return ad


def _unscale_uu(uu_list, ps):
    out = []
    for uu in uu_list:
        u_out = SimpleNamespace()
        u_out.S = ps.S*uu.S
        u_out.phi_c = ps.phi*uu.phi_c
        u_out.h = ps.h*uu.h
        u_out.phi = ps.phi*uu.phi
        u_out.u = ps.u*uu.u
        u_out.Q = ps.Q*uu.Q
        u_out.q = ps.q*uu.q
        u_out.N_c = ps.N*uu.N_c
        u_out.N = ps.N*uu.N
        u_out.M = ps.M*uu.M
        u_out.m = ps.m*uu.m
        u_out.kappa = ps.kappa*uu.kappa
        u_out.Sigma = ps.Sigma*uu.Sigma
        u_out.Sigma_c = ps.Sigma_c*uu.Sigma   # NB: MATLAB writes ps.Sigma_c*uu.Sigma (apparent typo, preserved for parity)
        u_out.V_in = ps.V_in*uu.V_in
        u_out.Sw = ps.S*uu.Sw
        u_out.hw = ps.h*uu.hw
        u_out.hu = ps.h*uu.hu
        if hasattr(uu, 'tau_b'):
            u_out.tau_b = uu.tau_b
        out.append(u_out)
    return out


# ---------- helpers (mirror nested functions) ----------

def _fill_factor(p_w, pp):
    return 1.0 - pp.a36*(1.0 - np.tanh((p_w - pp.a37)/pp.a37))/2.0


def _uplift(p, pp):
    return pp.a41*pp.a40*np.log1p(np.exp(p/pp.a40))


def _duplift(p, pp):
    e = np.exp(p/pp.a40)
    return pp.a41*e/(1.0 + e)


def _sliding(u, N, pp):
    return pp.a48*u**pp.q_slide*np.maximum(N, pp.a50)**pp.p_slide


# ---------- unpack solution vector ----------

def _UP_full(Y, aa, dd, pp, oo):
    """Return all derived fields. Y is 1D length 5*Iin."""
    xin, xxin = dd.xin, dd.xxin
    Iin = xin.size

    S = aa.S_ext.copy()
    phi_c = aa.phi_ext.copy()
    h = aa.h_ext.copy()
    phi = aa.phi_ext.copy()
    u = aa.u_ext.copy()

    t1 = 0
    t2 = t1 + Iin;  S[xin] = Y[t1:t2]
    t1 = t2; t2 += Iin; phi_c[xin] = Y[t1:t2]
    t1 = t2; t2 += Iin; h[xin] = Y[t1:t2]
    t1 = t2; t2 += Iin; phi[xin] = Y[t1:t2]
    t1 = t2; t2 += Iin; u[xxin] = Y[t1:t2]

    if oo.include_uplift:
        hu = _uplift(phi - aa.phi_s, pp)
    else:
        hu = np.zeros_like(h)

    if oo.include_fill_factor:
        hw = _fill_factor(phi - aa.phi_b, pp) * np.abs(h + hu)
        Sw = _fill_factor(phi_c - aa.phi_b, pp) * np.abs(S)
    else:
        hw = h
        Sw = S

    hav = dd.avxxin_d @ hw
    Sav = dd.avxxin_d @ Sw
    dphi = dd.ddxx_d @ phi
    dphi_c = dd.ddxx_d @ phi_c
    dphi_b = dd.ddxx_d @ aa.phi_b

    q = aa.q_ext.copy()
    Q = aa.Q_ext.copy()
    EPS = np.finfo(float).eps
    # MATLAB returns complex for negative^non-integer and silently propagates
    # NaNs that poison the Jacobian here. Use signed abs so transient negative
    # excursions don't crash the implicit solver.
    q_tmp = -pp.a29 * np.sign(hav) * np.abs(hav)**pp.alpha * dphi
    Q_tmp = -pp.a28 * np.sign(Sav) * np.abs(Sav)**pp.alpha_c * (dphi_c**2 + EPS**2)**((1.0/2.0)*(1.0/2.0 - 1)) * dphi_c
    q[xxin] = q_tmp[xxin]
    Q[xxin] = Q_tmp[xxin]

    N = aa.phi_s - phi
    N_c = aa.phi_s - phi_c
    kappa = pp.a27*(N_c - N) + pp.a39*hav**pp.alpha*(N_c - N)

    tau_b = _sliding(u, dd.avxx_d @ N, pp)

    if oo.couple_frictional_heating:
        m = pp.a9 + pp.a10*(tau_b*u) - pp.a11*q*dphi + pp.a12*q*(dphi - dphi_b)
    else:
        m = pp.a9 + pp.a10*(aa.tau_d*aa.U_b) - pp.a11*q*dphi + pp.a12*q*(dphi - dphi_b)
    M = -pp.a6*Q*dphi_c + pp.a7*Q*(dphi_c - dphi_b) + pp.a8*m

    m = dd.avxin_d @ m
    M = dd.avxin_d @ M

    Sigma = aa.sigma*(phi - aa.phi_b)
    Sigma_c = aa.sigma_c*(phi_c - aa.phi_b)
    S_in = np.zeros_like(dd.x); S_in[aa.xim] = aa.S_m
    V_in = S_in*(phi_c - aa.phi_b)

    return (S, phi_c, h, phi, u, Q, N_c, M, q, N, m, kappa,
            Sigma, Sigma_c, V_in, Sw, hw, hu, tau_b)


def _F_full(t, Y, aa, dd, pp, oo):
    """Right-hand side, before mass matrix. Returns length-5*Iin vector."""
    (S, phi_c, h, phi, u, Q, N_c, M, q, N, m, kappa,
     _Sigma, _Sigma_c, _V_in, _Sw, _hw, _hu, tau_b) = _UP_full(Y, aa, dd, pp, oo)

    Q_in = np.zeros_like(dd.x); Q_in[aa.xim] = aa.Q_m(t)
    M_in = aa.M_in(t)
    m_in = aa.m_in(t)

    F1 = pp.a2*M - pp.a3*S*np.abs(N_c)**(pp.n - 1)*N_c + pp.a4*(dd.avx_d @ u)*np.maximum(1.0 - pp.a5*S, 0.0)
    F2 = -pp.a18*(dd.ddx_d @ Q) + pp.a19*M + pp.a20*(Q_in/dd.dx) + pp.a21*kappa + pp.a34*M_in
    F3 = pp.a14*(dd.avx_d @ u)*np.maximum(1.0 - pp.a15*h, 0.0) - pp.a16*h*np.abs(N)**(pp.n - 1)*N + pp.a30*m
    F4 = -pp.a24*(dd.ddx_d @ (aa.W*q)) + pp.a25*aa.W*m - pp.a26*kappa + pp.a35*aa.W*m_in

    # ice momentum (F5). With include_ice_velocity off, zero out.
    if oo.include_ice_velocity:
        I = dd.I
        du = dd.ddx_d[:I-1, :] @ u
        inner = pp.a44*2*aa.H[:I-1]*(du**2 + pp.a49**2)**((1.0/pp.n - 1)/2.0) * du
        last = pp.a52*pp.a45*0.5*aa.H[I-1]**2
        vec = np.concatenate([inner, [last]])
        F5 = (dd.ddxx_d @ (aa.W*vec) + pp.a45*aa.tau_d*aa.W
              - pp.a46*tau_b*aa.W)
        if pp.a53 != 0:
            # Lateral-drag term ~ u^(1/n); use signed-|u| to avoid NaN when the
            # transient drives u briefly negative (MATLAB would return complex
            # and discard, here we'd otherwise poison the Jacobian).
            u_signed_pow = np.sign(u) * np.abs(u)**(1.0/pp.n)
            F5 -= (pp.a53*((pp.n + 2)/2.0)**(1.0/pp.n)
                   * u_signed_pow / aa.W**(1.0/pp.n) * aa.H)
    else:
        F5 = np.zeros(dd.I)

    return np.concatenate([F1[dd.xin], F2[dd.xin], F3[dd.xin], F4[dd.xin], F5[dd.xxin]])


def _mass_diagonals(Y, aa, dd, pp):
    """Return (a17, diag2, a22*W, diag4, a47*W) for the block-mass solve.
    diag2 corresponds to block 2 (phi_c eq) diagonal, diag4 to block 4 (phi eq).
    a1 and a13 are constant scalars and applied in _apply_mass_inv."""
    xin = dd.xin
    Iin = xin.size
    S_in = np.zeros_like(dd.x); S_in[aa.xim] = aa.S_m
    S = aa.S_ext.copy(); S[xin] = Y[0:Iin]
    h = aa.h_ext.copy(); h[xin] = Y[2*Iin:3*Iin]
    phi = aa.phi_ext.copy(); phi[xin] = Y[3*Iin:4*Iin]
    dhu = _duplift(phi - aa.phi_s, pp)

    diag2 = pp.a31*aa.sigma_c[xin]*aa.W[xin] + pp.a38*S[xin] + pp.a33*S_in[xin]/dd.dx[xin]
    diag4 = pp.a23*aa.sigma[xin]*aa.W[xin] + pp.a42*dhu[xin]*aa.W[xin] + pp.a43*h[xin]*aa.W[xin]
    return diag2, diag4


def _rhs(t, Y, aa, dd, pp, oo):
    """Return Y' = M(Y)^-1 F(t, Y) using the block-triangular structure."""
    F = _F_full(t, Y, aa, dd, pp, oo)
    Iin = dd.xin.size
    diag2, diag4 = _mass_diagonals(Y, aa, dd, pp)

    dY = np.empty_like(Y)
    # Block 1: a1 * dS/dt = F1
    dY[0:Iin] = F[0:Iin] / pp.a1
    # Block 2: a17 * dS/dt + diag2 * dphi_c/dt = F2
    dY[Iin:2*Iin] = (F[Iin:2*Iin] - pp.a17*dY[0:Iin]) / diag2
    # Block 3: a13 * dh/dt = F3
    dY[2*Iin:3*Iin] = F[2*Iin:3*Iin] / pp.a13
    # Block 4: a22*W * dh/dt + diag4 * dphi/dt = F4
    dY[3*Iin:4*Iin] = (F[3*Iin:4*Iin] - pp.a22*aa.W[dd.xin]*dY[2*Iin:3*Iin]) / diag4
    # Block 5: a47*W * du/dt = F5  (W indexed on edges = xxin)
    dY[4*Iin:5*Iin] = F[4*Iin:5*Iin] / (pp.a47*aa.W[dd.xxin])
    return dY


# ---------- initial condition ----------

def _initial(t, aa, dd, pp):
    Smin = 1e-3
    Q_in = np.zeros_like(dd.x); Q_in[aa.xim] = aa.Q_m(t)
    M_in = aa.M_in(t)
    m_in = aa.m_in(t)
    m = pp.a9 + pp.a10*aa.tau_d*aa.U_b
    M = np.zeros_like(Q_in)
    m = dd.avxin_d @ m
    M = dd.avxin_d @ M
    Q = np.cumsum(pp.a19*M + pp.a20*Q_in + pp.a34*M_in*dd.dx)
    phi_c = aa.phi_ext.copy()
    dphi_c = dd.ddxx_d @ phi_c
    EPS = np.finfo(float).eps
    S = np.abs(-Q/pp.a28 / ((dphi_c**2 + EPS**2)**(0.5*(0.5 - 1)) * dphi_c))**(1.0/pp.alpha_c)
    S = np.maximum(S, Smin)
    q = np.cumsum(pp.a25*m*dd.dx + pp.a35*m_in*dd.dx)
    phi = aa.phi_ext.copy()
    dphi = dd.ddxx_d @ phi
    h = np.abs(-q/pp.a29/dphi)**(1.0/pp.alpha)
    h = np.minimum(h, 1.0/pp.a15)
    u = aa.u_ext.copy()
    return S, phi_c, h, phi, u


def _pack(S, phi_c, h, phi, u, dd):
    return np.concatenate([S[dd.xin], phi_c[dd.xin], h[dd.xin], phi[dd.xin], u[dd.xxin]])


def _jac_pattern(dd):
    """5x5 block sparsity pattern. Each block is tri-diagonal (matching MATLAB)."""
    Iin = dd.xin.size
    rows = np.concatenate([np.arange(Iin), np.arange(Iin-1), np.arange(1, Iin)])
    cols = np.concatenate([np.arange(Iin), np.arange(1, Iin), np.arange(Iin-1)])
    data = np.ones(rows.size)
    tri = sp.csr_matrix((data, (rows, cols)), shape=(Iin, Iin))
    return sp.bmat([[tri]*5]*5, format='csr')


# ---------- main entry ----------

def wyvis_solve(td, ad, pd=None, ud=None, oo=None, verbose=True):
    """Mirror of wyvis_solve.m.

    td : 1D array of (dimensional) times to record output at.
    ad : SimpleNamespace with .x, .Z_b, .Z_s, .M_in (callable t -> array), etc.
    pd : parameters namespace (defaults from wyvis_default_parameters).
    ud : optional initial-condition list/namespace; if provided uses ud[-1].
    oo : options namespace.
    Returns (td, ad, ud_list).
    """
    if pd is None:
        pd = wyvis_default_parameters()
    if ud is None:
        ud = []
    if oo is None:
        oo = SimpleNamespace()

    ud_init = None
    if isinstance(ud, list) and len(ud) > 0:
        ud_init = ud[-1]
    elif hasattr(ud, 'phi'):
        ud_init = ud
    continuing = ud_init is not None

    defaults = dict(
        display_parameters=0,
        separate_timesteps=10,
        small_output=0,
        include_fill_factor=0,
        include_uplift=0,
        include_a4=1,
        include_a8=0,
        include_a11=0,
        include_a12=0,
        include_a44=1,
        include_a53=0,
        include_ice_velocity=0,
        couple_frictional_heating=0,
        # Integrator choice: 'BDF' (default, fast), 'Radau' (more robust on stiff
        # ice-velocity coupling), or 'LSODA' (auto stiff/non-stiff switch).
        method='BDF',
        rtol=1e-4,
        atol=1e-6,
    )
    for k, v in defaults.items():
        if not hasattr(oo, k):
            setattr(oo, k, v)

    if verbose:
        print('wyvis_solve: Discretizing ...')
    dd, aa, pp, ps = _scale_and_discretize(ad, pd, oo)

    td = np.asarray(td, dtype=float).reshape(-1)
    tt = td/ps.t

    if continuing:
        S0 = ud_init.S/ps.S
        phi_c0 = ud_init.phi_c/ps.phi
        h0 = ud_init.h/ps.h
        phi0 = ud_init.phi/ps.phi
        u0 = ud_init.u/ps.u
    else:
        S0, phi_c0, h0, phi0, u0 = _initial(tt[0], aa, dd, pp)

    Y = _pack(S0, phi_c0, h0, phi0, u0, dd)

    if verbose:
        print('wyvis_solve: Solving ...')

    jac_sparsity = _jac_pattern(dd)
    rtol, atol = oo.rtol, oo.atol
    # LSODA does not accept jac_sparsity; the other implicit methods do.
    ivp_kwargs = {} if oo.method == 'LSODA' else {'jac_sparsity': jac_sparsity}
    if getattr(oo, 'first_step', None) is not None:
        ivp_kwargs['first_step'] = oo.first_step
    if getattr(oo, 'max_step', None) is not None:
        ivp_kwargs['max_step'] = oo.max_step

    def rhs_call(t, y):
        return _rhs(t, y, aa, dd, pp, oo)

    nis = oo.separate_timesteps if oo.separate_timesteps > 0 else len(tt) - 1
    if nis == 0:
        # all at once
        sol = solve_ivp(rhs_call, (tt[0], tt[-1]), Y, t_eval=tt, method=oo.method,
                        rtol=rtol, atol=atol, **ivp_kwargs)
        if not sol.success:
            raise RuntimeError(f'solve_ivp failed: {sol.message}')
        tt = sol.t
        YY = sol.y
    else:
        import time as _time
        tt_out = [tt[0]]
        YY_cols = [Y.copy()]
        Y_curr = Y
        ti = 0
        while ti < len(tt) - 1:
            t0 = _time.time()
            if ti + 2*nis <= len(tt):
                seg = tt[ti:ti + nis + 1]
            else:
                seg = tt[ti:]
            sol = solve_ivp(rhs_call, (seg[0], seg[-1]), Y_curr, t_eval=seg,
                            method=oo.method, rtol=rtol, atol=atol, **ivp_kwargs)
            if not sol.success:
                raise RuntimeError(f'solve_ivp failed at t={seg[0]}: {sol.message}')
            tt_out.extend(sol.t[1:].tolist())
            YY_cols.extend([sol.y[:, k] for k in range(1, sol.y.shape[1])])
            Y_curr = sol.y[:, -1]
            ti += len(seg) - 1
            if verbose:
                dt = _time.time() - t0
                print(f'wyvis_solve:  t = {ps.t*seg[-1]/pd.td:g} / '
                      f'{ps.t*tt[-1]/pd.td:g} [{dt:.3f}s]')
        tt = np.array(tt_out)
        YY = np.array(YY_cols).T

    if verbose:
        print('wyvis_solve: Finalizing ...')

    uu = []
    for i in range(tt.size):
        y_i = YY[:, i]
        ui = SimpleNamespace()
        if oo.small_output:
            xin, xxin = dd.xin, dd.xxin
            Iin = xin.size
            S = aa.S_ext.copy(); S[xin] = y_i[0:Iin]
            phi_c = aa.phi_ext.copy(); phi_c[xin] = y_i[Iin:2*Iin]
            h = aa.h_ext.copy(); h[xin] = y_i[2*Iin:3*Iin]
            phi = aa.phi_ext.copy(); phi[xin] = y_i[3*Iin:4*Iin]
            u = aa.u_ext.copy(); u[xxin] = y_i[4*Iin:5*Iin]
            ui.S, ui.phi_c, ui.h, ui.phi, ui.u = S, phi_c, h, phi, u
        else:
            (ui.S, ui.phi_c, ui.h, ui.phi, ui.u, ui.Q, ui.N_c, ui.M, ui.q,
             ui.N, ui.m, ui.kappa, ui.Sigma, ui.Sigma_c, ui.V_in,
             ui.Sw, ui.hw, ui.hu, ui.tau_b) = _UP_full(y_i, aa, dd, pp, oo)
        # restore MATLAB convention: NaN at Dirichlet/Neumann boundaries
        ui.S[dd.xext] = np.nan
        ui.h[dd.xext] = np.nan
        if hasattr(ui, 'Sw'): ui.Sw[dd.xext] = np.nan
        if hasattr(ui, 'hw'): ui.hw[dd.xext] = np.nan
        uu.append(ui)

    td_out = ps.t*tt
    ad_out = _unscale_aa(aa, ps)
    ud_out = _unscale_uu(uu, ps)
    if verbose:
        print('wyvis_solve: Done')
    return td_out, ad_out, ud_out
