% flat bed and linear surface slope, driven by diurnal water input
% coupled drainage system, includes calculation of sliding velociy

% CLEAR WORKSPACE
clear; clc;

% DIMENTIONAL PARAMETERS
pd = wyvis_default_parameters;     % load default parameters
pd.lambda = 1e-9;                  % conduit-sheet exchange coefficient [m^2/s/Pa]
pd.h_rc = 0*0.01;                    % roughness height for channels [m]
disp(pd);                          % display parameters

% DIMENSIONAL GEOMETRY AND INPUTS
ad.x = linspace(0,10e3,100)';                            % x grid points [m]
b0 = 5*100;
ad.Z_b = 0*(ad.x-ad.x(end))-b0*ad.x.^0;                              % basal elevation [m]
ad.Z_s = 0.01*(ad.x(end)-ad.x)+(pd.rho_w/pd.rho_i-1)*b0;    % surface elevation [m]
ad.M_in = @(t) [ones(size(ad.x))]*(1e3);            % distributed conduit input [m^2/s]

ad.M_in = @(t) 1*[ones(size(ad.x))]*(1e-4);            % distributed conduit input [m^2/s]
ad.x_m = 5e3;                            % moulin poision [m]
ad.Q_m = @(t) 10; %*(1-exp(-t/(10*pd.td)));   % moulin source [m^3/s]

% SOLVE FOR 1 YEAR TOWARDS STEADY STATE
oo.include_ice_velocity = 1;                            % include ice velocity in calculation
[td,ad,ud] = wyvis_solve( [0:1:365]*24*60*60 , ad , pd , [] , oo); % solve model
wyvis_plot_all( 1 , ad , ud(end) );                     % plot all variables at final timestep

return;

% SOLVE FOR 10 DAYS WITH SMALL DIURNAL OSCILLATION
ad.M_in = @(t) [ones(size(ad.x))]*1e-6*max(0, 1 + 1*sin(2*pi*t/24/60/60) );    % update distributed conduit input [m^2/s]
[td,ad,ud] = wyvis_solve( [0:1/24:10]*24*60*60 , ad , pd , ud(end) , oo );  % solve model
wyvis_plot_all( 1 , ad , ud(end) );                     % plot all variables at final timestep

% PLOT RESULTS
wyvis_plot_allt( 2 , td , ud );       % plot all variables through time
% wyvis_movie_u( 3 , ad , ud , td );  % movie of ice velocity

% contour plot of ice velocity in m/y
figure(4); 
    contourf( [ad.x]/1e3 , td/pd.td , [ud.u]'*pd.ty ); colorbar;           
    xlabel('Distance [km]');
    ylabel('Time [d]');
