% topography taken from flow line in greenland,
% coupled drainage system, fed by constant input to channel

% CLEAR WORKSPACE
clear; clc;

% DIMENTIONAL PARAMETERS
pd = wyvis_default_parameters;     % load default parameters
pd.h_rc = 0.01;                    % roughness height for channels [m]
pd.sigma = 1e-3;                   % connected englacial void fraction []
pd.lambda = 1e-9;                  % conduit-sheet exchange coefficient [m^2/s/Pa]
pd.K = 1e-4;                       % flow coefficient in sheet [kg^(-1)m^(4-alpha)s]
pd.C_slide = 3e4;                  % constant in sliding law [s/m]
pd.N_min = 1e5;                    % minimum effective pressure in sliding law [Pa]
disp(pd);                          % display parameters

% DIMENSIONAL GEOMETRY AND INPUTS
load('greenland_slice2');
fi = 2; ts = fl(fi).ts; bs = fl(fi).bs; ss = fl(fi).ss;
x = linspace(ts(1),ts(end),100)'; b = interp1(ts,bs,x); s = interp1(ts,ss,x); 
N = 3; smooth = spdiags(ones(length(x),1+2*N),-N:N,length(x),length(x)); smooth = spdiags(sum(smooth,2).^(-1),0,length(x),length(x))*smooth;
b = smooth*b; s = smooth*s; % smooth
ad.x = x(end)-flip(x);     % x grid points [m]
ad.Z_b = b;                % basal elevation [m]
ad.Z_s = s;                % surface elevation [m]  
ad.M_in = @(t) ad.x.^0*(1e-6*t.^0);       % distributed input [m^2/s] 

wyvis_plot_geometry( 3 , ad ); shg;
% return;

% SOLVE TOWARDS STEADY STATE
oo.include_ice_velocity = 1;                            % include ice velocity in calculation
oo.include_a53 = 0;                                     % include lateral drag
[td,ad,ud] = wyvis_solve( [0:1:2*365]*24*60*60 , ad , pd , [] , oo); % solve model
wyvis_plot_all( 1 , ad , ud(end) );                     % plot all variables at final timestep
wyvis_plot_allt( 2 , td , ud );                       % plot all variables through time
ud0 = ud(end);                                          % remember solution for later use

return;

% RUN FOR ANNUAL CYCLE
ad.M_in = @(t) [ad.x.^0]*( 1e-6+1e-3*wyvis_annual_signal(t,135*pd.td,244*pd.td,21*pd.td) );  % update distributed conduit input 
[td,ad,ud] = wyvis_solve( pd.td*[0:7:1*365] , ad , pd , ud0  ,oo );   % solve model

% PLOT RESULTS
wyvis_plot_allt( 2 , td , ud );       % plot all variables through time

% contour plot of ice velocity in m/y
figure(4); 
    contourf( [ad.x]/1e3 , td/pd.td , [ud.u]'*pd.ty ); colorbar;           
    xlabel('Distance [km]');
    ylabel('Time [d]');
