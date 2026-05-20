function pd = wyvis_default_parameters()
% pd = wyvis_default_parameters()
% default parameter values

if nargin<1, pd = struct; end

% DIMENSIONAL PARAMETERS
if ~isfield(pd,'rho_w'), pd.rho_w = 1e3; end                % water density [kg/m^3]
if ~isfield(pd,'rho_i'), pd.rho_i = .9e3; end               % ice density [kg/m^3]
if ~isfield(pd,'d'), pd.g = 9.8; end                        % gravity [m/s^2]
if ~isfield(pd,'L'), pd.L = 3.35e5; end                     % latent heat [J/kg]
if ~isfield(pd,'c'), pd.c = 4.2e3; end                      % specific heat capacity [J/kg/K]
if ~isfield(pd,'gamma'), pd.gamma = 7.5e-8; end             % melting point pressure gradient [K/Pa]
if ~isfield(pd,'G'), pd.G = 6e-2; end                       % geothermal heat [W/m^2]
if ~isfield(pd,'td'), pd.td = 24*60*60; end                 % time in a day [s]
if ~isfield(pd,'ty'), pd.ty = 365*24*60*60; end             % time in a year [s]
if ~isfield(pd,'alpha_c'), pd.alpha_c = 4/3; end            % flow exponent in channel 
if ~isfield(pd,'K_c'), pd.K_c = 0.05; end                   % flow coefficient in channel [kg^(-1/2)m^(4-2*alpha_c) ]
if ~isfield(pd,'alpha'), pd.alpha = 3; end                  % flow exponent in sheet
if ~isfield(pd,'K'), pd.K = 1e-6; end                       % flow coefficient in sheet [kg^(-1)m^(4-alpha)s]
if ~isfield(pd,'n'), pd.n = 3; end                          % ice viscosity exponent
if ~isfield(pd,'A'), pd.A = 6.8e-24; end                    % ice flow law coefficent [Pa^(-n)/s]
if ~isfield(pd,'Atil'), pd.Atil = 10*2*pd.A/pd.n^pd.n; end  % creep coefficient in channel [Pa^(-n)/s]
if ~isfield(pd,'Ahat'), pd.Ahat = 10*2*pd.A/pd.n^pd.n; end  % creep coefficient in sheet [Pa^(-n)/s]
if ~isfield(pd,'h_r'), pd.h_r = 1; end                      % bed roughness height [m]
if ~isfield(pd,'l_r'), pd.l_r = 10; end                     % bed roughness wavelength [m]
if ~isfield(pd,'h_rc'), pd.h_rc = 0; end                    % channel roughness height [m]
if ~isfield(pd,'l_rc'), pd.l_rc = 10; end                    % channel roughness length [m]
if ~isfield(pd,'sigma'), pd.sigma = 0; end                  % storage coefficient connected to sheet
if ~isfield(pd,'sigma_c'), pd.sigma_c = 0; end              % storage coefficient connected to channel
if ~isfield(pd,'lambda'), pd.lambda = 0; end                % exchange coefficient [m^2/s/Pa]
if ~isfield(pd,'U_b'), pd.U_b = 1e2/pd.ty; end              % sliding speed [m/s]
if ~isfield(pd,'W'), pd.W = 1e3; end                        % flowline width [m]
if ~isfield(pd,'beta_c'), pd.beta_c = 5e-10; end            % channel compressibility [Pa^(-1)] (numerical artifice)
if ~isfield(pd,'beta'), pd.beta = 5e-10; end                % sheet compressibility [Pa^(-1)] (numerical artifice)
if ~isfield(pd,'p_ff'), pd.p_ff = 10e4; end                 % pressure scale for fill factor regularization [Pa] (numerical artifice)
if ~isfield(pd,'lambda2'), pd.lambda2 = 0; end              % exchange coefficient [m^(-1)]
if ~isfield(pd,'p_u'), pd.p_u = 1e4; end                    % regularization pressure range for uplift [Pa]
if ~isfield(pd,'E_u'), pd.E_u = 0; end                      % displacement of uplift [m/Pa]
if ~isfield(pd,'p_slide'), pd.p_slide = 1; end              % exponent of N in sliding law
if ~isfield(pd,'q_slide'), pd.q_slide = 1; end              % exponent of u in sliding law
if ~isfield(pd,'C_slide'), pd.C_slide = 3.2e4; end          % constant in sliding law [s/m]
if ~isfield(pd,'Abar'), pd.Abar = pd.A/10; end                 % average Glen's law coefficient [Pa^(-n)/s]
if ~isfield(pd,'eps_reg'), pd.eps_reg = 1e-13; end          % regularizing strain rate in viscosity [1/s]
if ~isfield(pd,'Re_reg'), pd.Re_reg = 1e12; end             % density in ice momentum equation (artificial inertial term) [kg/m^3]
if ~isfield(pd,'N_min'), pd.N_min = 0*10^4; end               % minimum effective pressure in sliding law [Pa]
if ~isfield(pd,'phi_m'), pd.phi_m = 0; end                  % minimum hydraulic potential (corresponding to sea level) [Pa]
if ~isfield(pd,'rho_o'), pd.rho_o = pd.rho_w; end           % ocean density [kg/m^3]

end