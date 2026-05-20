% flat bed and linear surface slope, driven by diurnal water input, 
% conduit drainage system

% CLEAR WORKSPACE
clear; clc;

% DIMENSIONAL PARAMETERS
pd = wyvis_default_parameters;                  % load default parameters
disp(pd);                                       % display parameters

% DIMENSIONAL GEOMETRY AND INPUTS
ad.x = linspace(0,10000,40)';                   % x grid points [m]
ad.Z_b = 0*ad.x;                                % bed elevation [m]
ad.Z_s = 0.1*(ad.x(end)-ad.x);                  % surface elevation [m]
ad.M_in = @(t) [ones(size(ad.x))]*1e-3*t.^0;    % distributed conduit input [m^2/s]

% SOLVE FOR 1 YEAR TOWARDS STEADY STATE
[td,ad,ud] = wyvis_solve( [0:1:365]*24*60*60 , ad , pd );   % solve model
wyvis_plot_all_c( 1 , ad , ud(end) );           % plot channel variables at final timestep
% figure(2); plot( td , [ud.phi_c] );            % check time evolution

% return;

% SOLVE FOR 10 DAYS WITH SMALL DIURNAL OSCILLATION
ad.M_in = @(t) [ones(size(ad.x))]*1e-3*( 1 + 0.5*sin(2*pi*t/24/60/60) );    % update distributed conduit input [m^2/s]
[td,ad,ud] = wyvis_solve( [0:1/24:10]*24*60*60 , ad , pd , ud(end) );       % solve model
wyvis_plot_all_c( 1 , ad , ud(end) );        % plot channel variables at final timestep

% PLOT RESULTS
wyvis_plot_phi_ct( 1 , td , ud , 20 );      % plot channel potential through time at particular grid point
% wyvis_movie_phi_c( 2 , ad , ud , td );      % movie of channel potential through time
