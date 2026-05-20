"""Default dimensional parameters for the wyvis subglacial hydrology model.

Direct port of wyvis_default_parameters.m. Returns a plain dataclass-like
namespace so attribute access matches MATLAB struct usage (pd.rho_w etc).
"""
from types import SimpleNamespace


def wyvis_default_parameters(pd=None):
    if pd is None:
        pd = SimpleNamespace()

    def setdef(name, value):
        if not hasattr(pd, name):
            setattr(pd, name, value)

    setdef('rho_w', 1e3)        # water density [kg/m^3]
    setdef('rho_i', 0.9e3)      # ice density [kg/m^3]
    setdef('g', 9.8)            # gravity [m/s^2]  (matlab sets pd.g under "d" branch)
    setdef('L', 3.35e5)         # latent heat [J/kg]
    setdef('c', 4.2e3)          # specific heat capacity
    setdef('gamma', 7.5e-8)     # melting point pressure gradient [K/Pa]
    setdef('G', 6e-2)           # geothermal heat [W/m^2]
    setdef('td', 24*60*60)
    setdef('ty', 365*24*60*60)
    setdef('alpha_c', 4/3)
    setdef('K_c', 0.05)
    setdef('alpha', 3)
    setdef('K', 1e-6)
    setdef('n', 3)
    setdef('A', 6.8e-24)
    setdef('Atil', 10*2*pd.A/pd.n**pd.n)
    setdef('Ahat', 10*2*pd.A/pd.n**pd.n)
    setdef('h_r', 1.0)
    setdef('l_r', 10.0)
    setdef('h_rc', 0.0)
    setdef('l_rc', 10.0)
    setdef('sigma', 0.0)
    setdef('sigma_c', 0.0)
    setdef('lambda_', 0.0)      # 'lambda' is reserved in python
    setdef('U_b', 1e2/pd.ty)
    setdef('W', 1e3)
    setdef('beta_c', 5e-10)
    setdef('beta', 5e-10)
    setdef('p_ff', 10e4)
    setdef('lambda2', 0.0)
    setdef('p_u', 1e4)
    setdef('E_u', 0.0)
    setdef('p_slide', 1)
    setdef('q_slide', 1)
    setdef('C_slide', 3.2e4)
    setdef('Abar', pd.A/10)
    setdef('eps_reg', 1e-13)
    setdef('Re_reg', 1e12)
    setdef('N_min', 0*10**4)
    setdef('phi_m', 0.0)
    setdef('rho_o', pd.rho_w)
    return pd
