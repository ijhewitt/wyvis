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

% ADD MOULIN
ad.x_m = 0;          % moulin location
ad.S_m = 100;         % moulin area
ad.Q_m = @(t) 10*(1-exp(-t/(7/365*pd.ty)));           % moulin flux (ramping from 0 to steady value over time)

% SOLVE FOR 1 YEAR TOWARDS STEADY STATE
[td,ad,ud] = wyvis_solve( [0:1:365]*24*60*60 , ad , pd , ud(end) );   % solve model
wyvis_plot_all_c( 1 , ad , ud(end) );           % plot channel variables at final timestep
% figure(2); plot( td , [ud.phi_c] );            % check time evolution

% PLOT CHANNEL DISCHARGE OVER TIME
figure(3); clf; 
    tmp = contourf( [ad.x]/1e3 , td/pd.td , [ud.Q]' , 'linestyle','none' ); colorbar;       
    xlabel('Distance [km]');
    ylabel('Time [d]');