% glacier scale flat bed and sqrt surface elevation, 
% seasonal water input, coupled drainage system

% CLEAR WORKSPACE
clear; clc;

% DIMENSIONAL PARAMETERS
pd = wyvis_default_parameters;      % default parameters
pd.h_rc = 0.01;                     % roughness height for channels [m]
pd.lambda = 1e-9;                   % conduit-sheet exchange coefficient [m^2/s/Pa]
pd.sigma = 0*1e-3;                    % connected englacial void fraction []
pd.K = 1e-5;                        % sheet permeability [kg^(-1)m^(4-alpha)s]
pd.W = 10e3;                        % width
disp(pd);                           % display parameters

% DIMENSIONAL GEOMETRY AND INPUTS
ad.x = linspace(0,40*1e3,80)';            % x grid points [m]            
ad.Z_b = -400*ad.x.^0;                       % bed elevation [m]  
ad.Z_s = -(pd.rho_w/pd.rho_i-1)*ad.Z_b(end)+1000*(1-ad.x/ad.x(end)).^(1/2);  % surface elevation [m]
ad.M_in = @(t) ad.x.^0*(1e-8*t.^0);       % distributed conduit input [m^2/s]

% RUN FOR 1 YEAR TOWARDS STEADY STATE
oo.include_ice_velocity = 1;     % include calculation of ice velocty
oo.couple_frictional_heating = 1;   % include frictional heat from velocity
oo.include_a53 = 0;              % include lateral drag
[td,ad,ud] = wyvis_solve( [0:10:1*365]*pd.td , ad , pd , [] , oo);  % solve model
ud0 = ud(end);                   % remember solution for later use

% ad.x_m = [70*1e3]; 
% ad.S_m = 0;
% ad.Q_m = @(t) [ad.x_m.^0]*( 1*wyvis_annual_signal(t,135*pd.td,244*pd.td,21*pd.td) );  % update distributed conduit input 

% RUN FOR SEVERAL YEARS
ad.M_in = @(t) [ad.x.^0]*( 1e-8+1e-3*wyvis_annual_signal(t,135*pd.td,244*pd.td,21*pd.td) );  % update distributed conduit input 
[td,ad,ud] = wyvis_solve( pd.td*[0:7:2*365] , ad , pd , ud0  ,oo );   % solve model
td1 = td; ud1 = ud;             % remember solution for later use

% PLOT RESULTS
wyvis_plot_allt(1,td,ud);       % plot of all variables through time
% wyvis_movie_phi(2,ad,ud,td);  % movie of potential 

% plot of total input and output
Q_in = wyvis_total_input(td,ad);
Q_out = [ud.Q]; Q_out = Q_out(end,:);
q_out = [ud.q]; q_out = ad.W(end)*q_out(end,:);
figure(2); clf; plot(td/pd.td,Q_in,'k',td/pd.td,Q_out,'ro-',td/pd.td,q_out,'bo-');
    xlabel('Time [ d ]','interpreter','latex'); ylabel('Discharge [ m${}^3$/s ]','interpreter','latex');

% contour plot of effective pressure
figure(3); clf; imagesc(ad.x/1e3,td/pd.td,[ud.N]'/1e6); colorbar; set(gca,'YDir','normal');
  xlabel('Distance [ km ]','interpreter','latex'); ylabel('Time [ d ]','interpreter','latex'); title('Effective pressure [ MPa ]','interpreter','latex');
  
% contour plot of channel area
figure(4); clf; imagesc(ad.x/1e3,td/pd.td,[ud.S]'); colorbar; set(gca,'YDir','normal');
  xlabel('Distance [ km ]','interpreter','latex'); ylabel('Time [ d ]','interpreter','latex'); title('Cross-sectional area [ m${}^2$ ]','interpreter','latex');

% contour plot of ice velocity
figure(5); clf; imagesc(ad.x/1e3,td/pd.td,([ud.u]'*(365*24*60*60))); colorbar; set(gca,'YDir','normal');
  xlabel('Distance [ km ]','interpreter','latex'); ylabel('Time [ d ]','interpreter','latex'); title('Ice velocity [ m/y ]','interpreter','latex');

