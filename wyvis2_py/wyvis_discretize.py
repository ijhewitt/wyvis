"""1D finite-volume discretization operators on a non-uniform node grid.

Port of wyvis_discretize from wyvis_solve.m (nested subfunction).
All indices are 0-based here; MATLAB's 1-based xin/xxin/xext/xxext are
shifted by -1.

Nodes:   x[0..I-1]     (cell centers)
Edges:   xx[0..I-1]    (cell faces; xx[0] = x[0], xx[i] = (x[i-1]+x[i])/2)
Cell widths dx[i] :    distance between edge i and edge i+1 (last cell width handled specially)
Node spacings dxx[i] : distance from previous node to node i (dxx[0] = x[0])
"""
from types import SimpleNamespace
import numpy as np
import scipy.sparse as sp


def wyvis_discretize(x):
    x = np.asarray(x, dtype=float).reshape(-1)
    I = x.size
    dd = SimpleNamespace()
    dd.I = I
    dd.x = x

    # xx[0] = x[0]; xx[i] = (x[i-1] + x[i])/2 for i=1..I-1
    xx = np.empty(I)
    xx[0] = x[0]
    xx[1:] = 0.5*(x[:-1] + x[1:])
    dd.xx = xx

    # dx[i] = xx[i+1]-xx[i] for i=0..I-2 ; dx[I-1] = x[-1]-xx[I-1]
    dx = np.empty(I)
    dx[:-1] = xx[1:] - xx[:-1]
    dx[-1] = x[-1] - xx[-1]
    dd.dx = dx

    # dxx[i] = x[i] - x[i-1] for i=1..I-1; dxx[0] = x[0]
    dxx = np.empty(I)
    dxx[0] = x[0]
    dxx[1:] = x[1:] - x[:-1]
    dd.dxx = dxx

    # ddx: divergence-like operator. ddx[i,i] = -1/dx[i] for i=0..I-1,
    # ddx[i,i+1] = +1/dx[i] for i=0..I-2.
    rows = np.concatenate([np.arange(I), np.arange(I-1)])
    cols = np.concatenate([np.arange(I), np.arange(1, I)])
    vals = np.concatenate([-1.0/dx, 1.0/dx[:-1]])
    dd.ddx = sp.csr_matrix((vals, (rows, cols)), shape=(I, I))

    # ddxx: edge-gradient operator with extrapolation on left boundary.
    # ddxx[0,0] = -1/dxx[1], ddxx[0,1] = +1/dxx[1]
    # ddxx[i,i-1] = -1/dxx[i], ddxx[i,i] = +1/dxx[i] for i=1..I-1
    dxx_used = dxx.copy()
    dxx_used[0] = dxx[1]  # extrapolated using dxx[1] for the first edge
    rows = np.concatenate([np.arange(I), np.arange(I)])
    cols = np.concatenate(
        [np.array([0]), np.arange(0, I-1), np.array([1]), np.arange(1, I)]
    )
    vals = np.concatenate([-1.0/dxx_used, 1.0/dxx_used])
    dd.ddxx = sp.csr_matrix((vals, (rows, cols)), shape=(I, I))

    # avx: average from interior edges to nodes
    # avx[i,i] = 0.5 for i=0..I-2; avx[I-1,I-1] = 1; avx[i,i+1] = 0.5 for i=0..I-2
    rows = np.concatenate([np.arange(I), np.arange(I-1)])
    cols = np.concatenate([np.arange(I), np.arange(1, I)])
    diag_vals = np.concatenate([0.5*np.ones(I-1), [1.0]])
    off_vals = 0.5*np.ones(I-1)
    vals = np.concatenate([diag_vals, off_vals])
    dd.avx = sp.csr_matrix((vals, (rows, cols)), shape=(I, I))

    # avxx: average from interior nodes to edges
    # avxx[i,i-1] = 0.5 for i=1..I-1; avxx[0,0] = 1; avxx[i,i] = 0.5 for i=1..I-1
    rows = np.concatenate([np.arange(1, I), np.arange(I)])
    cols = np.concatenate([np.arange(I-1), np.arange(I)])
    sub_vals = 0.5*np.ones(I-1)
    diag_vals = np.concatenate([[1.0], 0.5*np.ones(I-1)])
    vals = np.concatenate([sub_vals, diag_vals])
    dd.avxx = sp.csr_matrix((vals, (rows, cols)), shape=(I, I))

    # Boundary index sets (0-based)
    dd.xext = np.array([I-1])           # Dirichlet nodes
    dd.xxext = np.array([0])            # Neumann edges
    dd.xin = np.setdiff1d(np.arange(I), dd.xext)
    dd.xxin = np.setdiff1d(np.arange(I), dd.xxext)

    # avxin: normalize avx so each row sums to 1 using interior edges only
    tmp = np.asarray(dd.avx[:, dd.xxin].sum(axis=1)).ravel()
    inv_tmp = np.where(tmp != 0, 1.0/np.where(tmp == 0, 1, tmp), 0.0)
    avxin = sp.diags(inv_tmp).dot(dd.avx).tolil()
    avxin[:, dd.xxext] = 0
    dd.avxin = avxin.tocsr()

    tmp = np.asarray(dd.avxx[:, dd.xin].sum(axis=1)).ravel()
    inv_tmp = np.where(tmp != 0, 1.0/np.where(tmp == 0, 1, tmp), 0.0)
    avxxin = sp.diags(inv_tmp).dot(dd.avxx).tolil()
    avxxin[:, dd.xext] = 0
    dd.avxxin = avxxin.tocsr()

    # Dense copies for the hot RHS path — for the small grids the model
    # targets (~40 nodes) dense matvec beats sparse by ~2.5x.
    dd.ddx_d    = dd.ddx.toarray()
    dd.ddxx_d   = dd.ddxx.toarray()
    dd.avx_d    = dd.avx.toarray()
    dd.avxx_d   = dd.avxx.toarray()
    dd.avxin_d  = dd.avxin.toarray()
    dd.avxxin_d = dd.avxxin.toarray()

    return dd


def wyvis_nearest_gridpoint(x, xd):
    """Indices (0-based) of nearest grid points in xd to each entry of x."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    xd = np.asarray(xd, dtype=float)
    if x.size == 0:
        return np.array([], dtype=int)
    return np.array([int(np.argmin((xi - xd)**2)) for xi in x])
