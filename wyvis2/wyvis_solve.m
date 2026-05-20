function [td,ad,ud] = wyvis_solve(td,ad,pd,ud,oo)
% [td,ad,ud] = wyvis_solve(td,ad,pd,ud,oo)
% main component of wyvis model 
% INPUTS
%   td      time range and timesteps for output [1-by-L]
%   ad      struct containing geometry and inputs 
%   pd      struct containing parameters
%   ud      struct containing initial condition (optional)
%   oo      struct containing run options (optional)
% OUTPUTS
%   td      time range and timesteps for output [1-by-L]
%   ad      struct containing geometry and inputs including derived fields
%   ud      struct containing output [1-by-L]
%
% IJH 21 August 2013 
%     9 Sept 2016 updated to include ice velocity
%     11 Sept 2019 updated to include lateral drag
%     26 May 2023 updated to allow coupling of frictional heat

    % OPTIONS    
    if nargin<3, pd = wyvis_default_parameters(); end
    if nargin<4, ud = struct; end
    if nargin<5, oo = struct; end
    if isfield(ud,'phi'), oo.continuing = 1; else oo.continuing = 0; end    % initial condition provided
    if ~isfield(oo,'display_parmaters'), oo.display_parameters = 0; end     % display dimensionless parameters and scaling
    if ~isfield(oo,'separate_timesteps'), oo.separate_timesteps = 10; end   % solve this number of timesteps at each call to ode15s (0 for all timesteps at once, which should be quickest)
    if ~isfield(oo,'small_output'), oo.small_output = 0; end                % output only primary variables S,phi_c,h,phi
    if ~isfield(oo,'include_fill_factor'), oo.include_fill_factor = 0; end  % fill factor regularization for atmospheric pressure
    if ~isfield(oo,'include_uplift'), oo.include_uplift = 0; end            % uplift above overburden
    if ~isfield(oo,'include_a4'), oo.include_a4 = 1; end                    % sliding opening of channel
    if ~isfield(oo,'include_a8'), oo.include_a8 = 0; end                    % sheet contribution to channel melting
    if ~isfield(oo,'include_a11'), oo.include_a11 = 0; end                  % dissipation in sheet 
    if ~isfield(oo,'include_a12'), oo.include_a12 = 0; end                  % pressure melting term in sheet
    if ~isfield(oo,'include_a44'), oo.include_a44 = 1; end                  % longitudinal stress
    if ~isfield(oo,'include_a53'), oo.include_a53 = 0; end                  % lateral drag
    if ~isfield(oo,'include_ice_velocity'), oo.include_ice_velocity = 0; end% solve for changes in ice velocity
    if ~isfield(oo,'couple_frictional_heating'), oo.couple_frictional_heating = 0; end% calculate frictional heat using solved-for ice velocity and stress (rather than driving stress and aa.U_b)

    % SCALE AND DISCRETIZE
    disp('wyvis_solve: Discretizing ...');
    [dd,aa,pp,ps] = wyvis_scale_and_discretize(ad,pd,oo);
    tt = td/ps.t; 
    if oo.display_parameters, disp('wyvis_solve: Scales:'); disp(ps); disp('wyvis_solve: Parameters:'); disp(pp); end
    
    % INITIAL CONDITION
    if oo.continuing
        Y = wyvis_P(ud(end).S/ps.S,ud(end).phi_c/ps.phi,ud(end).h/ps.h,ud(end).phi/ps.phi,ud(end).u/ps.u);  
    else 
        [S,phi_c,h,phi,u] = wyvis_initial(tt(1),aa,dd);
        Y = wyvis_P(S,phi_c,h,phi,u); 
    end
    
    % SOLVE
    disp('wyvis_solve: Solving ...');
    mass = @(t,Y) wyvis_mass(Y,aa,dd,pp);     
    opts = odeset('Mass',mass,'MvPattern',wyvis_mvpattern(dd),'JPattern',wyvis_jpattern(dd),'Vectorized','on','Reltol',1e-4,'abstol',1e-6); % ode15s options  
    if oo.separate_timesteps==1,
    t2 = tt(1); Y2 = Y; Y1 = Y;
    for ti = 1:length(tt)-1,
        tic
        [t1,Y1] = ode15s(@wyvis_F,linspace(tt(ti),tt(ti+1),3),Y1,opts);
        t1 = t1(end); Y1 = Y1(end,:)'; 
        t2 = [t2 t1]; Y2 = [Y2 Y1];
        comp_time = toc;
        disp(['wyvis_solve:  t = ',num2str(ps.t*t1/pd.td),' / ',num2str(ps.t*tt(end)/pd.td),' [',num2str(comp_time),'s]']);
    end
    tt = t2; YY = Y2;
    elseif oo.separate_timesteps>1
    t2 = tt(1); Y2 = Y; Y1 = Y; ti = 1;
    nis = oo.separate_timesteps;
    while ti<length(tt)
        tic;
        if ti+2*nis<=length(tt), t1 = tt(ti:(ti+nis));
        else nis = length(tt)-ti; t1 = tt(ti:(ti+nis)); 
        end
        [t1,Y1] = ode15s(@wyvis_F,t1,Y1,opts);
        t2 = [t2 t1(2:end)']; Y2 = [Y2 Y1(2:end,:)'];
        ti = ti+length(t1)-1; t1 = t1(end); Y1 = Y1(end,:)';  
        comp_time = toc;
        disp(['wyvis_solve:  t = ',num2str(ps.t*t1/pd.td),' / ',num2str(ps.t*tt(end)/pd.td),' [',num2str(comp_time),'s]']);
    end
    tt = t2; YY = Y2;
    else
    [tt,YY] = ode15s(@wyvis_F,tt,Y,opts);
    tt = tt'; YY = YY';
    end

    % EXTRACT VARIABLES
    disp('wyvis_solve: Finalizing ...');
    for ti = 1:length(tt)
        uui = struct;
        if oo.small_output, [uui.S,uui.phi_c,uui.h,uui.phi,uui.u] = wyvis_UP(YY(:,ti),aa,dd);
        else [uui.S,uui.phi_c,uui.h,uui.phi,uui.u,uui.Q,uui.N_c,uui.M,uui.q,uui.N,uui.m,uui.kappa,uui.Sigma,uui.Sigma_c,uui.V_in,uui.Sw,uui.hw,uui.hu,uui.tau_b] = wyvis_UP_full(YY(:,ti),aa,dd,pp);
        end
        if ti==1, uu = uui;
        else uu(ti) = uui;
        end
    end
    
    % UNSCALE
    td = ps.t*tt;
    ad = wyvis_unscale_aa(aa,ps);
    ud = wyvis_unscale_uu(uu,ps);
    
    disp('wyvis_solve: Done');
    
    
%% SUBFUNCTIONS
function [dd,aa,pp,ps] = wyvis_scale_and_discretize(ad,pd,oo)
    
% FILL IN OPTIONAL FIELDS WITH DEFAULTS   
if ~isfield(ad,'Z_b') && isfield(ad,'b'), ad.Z_b = ad.b; end  % old bed elevation
if ~isfield(ad,'Z_s') && isfield(ad,'s'), ad.Z_s = ad.s; end  % old surface elevation
if ~isfield(ad,'W'), ad.W = pd.W*ad.x.^0; end               % flowline width [m]
if ~isfield(ad,'U_b'), ad.U_b = pd.U_b*ad.x.^0; end         % sliding speed [m/s]
if ~isfield(ad,'M_in'), ad.M_in = @(t) 0*ad.x.^0*t; end     % channel source [m^2/s]
if ~isfield(ad,'m_in'), ad.m_in = @(t) 0*ad.x.^0*t; end     % sheet source [m/s]
if ~isfield(ad,'x_m'), ad.x_m = []; end                     % moulin positions [m]
if ~isfield(ad,'Q_m'), ad.Q_m = @(t) 0*ad.x_m.^0*t; end     % moulin inputs [m^3/s]
if ~isfield(ad,'S_m'), ad.S_m = 0*ad.x_m.^0; end            % moulin areas [m^2]

% SCALING PARAMETERS
ps.x = 10*1e3;
ps.s = 1000;
ps.u = pd.U_b;
ps.Q = 1;
ps.W = 1e3;
ps.tau = pd.rho_i*pd.g*ps.s^2/ps.x;
ps.phi = pd.rho_i*pd.g*ps.s;
ps.S = (ps.Q/pd.K_c/(ps.phi/ps.x)^(1/2))^(1/pd.alpha_c);
ps.M = ps.Q*ps.phi/ps.x/pd.L/pd.rho_w;
ps.t = pd.rho_i*ps.S/ps.M/pd.rho_w;
ps.kappa = ps.Q/ps.x;
ps.q = ps.Q/ps.W;
ps.h = (ps.q/pd.K/(ps.phi/ps.x))^(1/pd.alpha);
ps.N = ps.phi;
ps.m = ps.q/ps.x;
ps.M_in = ps.Q/ps.x;
ps.m_in = ps.q/ps.x;
ps.sigma = pd.sigma;
ps.sigma_c = pd.sigma_c;
ps.Sigma = ps.sigma*ps.phi/pd.rho_w/pd.g;
ps.Sigma_c = ps.sigma_c*ps.phi/pd.rho_w/pd.g;
ps.V_in = ps.S*ps.phi/pd.rho_w/pd.g;

% DIMENSIONLESS PARAMETERS
pp.alpha_c = pd.alpha_c;
pp.alpha = pd.alpha;
pp.n = pd.n;
pp.p_slide = pd.p_slide;
pp.q_slide = pd.q_slide;
pp.a1 = 1;
pp.a2 = pd.rho_w*ps.M*ps.t/pd.rho_i/ps.S;
pp.a3 = pd.Atil*ps.N^pd.n*ps.t;
pp.a4 = ps.u*pd.h_rc*ps.t/ps.S;
pp.a5 = ps.S/pd.h_rc/pd.l_rc;
pp.a6 = ps.Q*ps.phi/ps.x/pd.L/ps.M/pd.rho_w;
pp.a7 = ps.Q*ps.phi/ps.x/pd.L/ps.M/pd.rho_w*pd.gamma*pd.rho_w*pd.c;
pp.a8 = ps.m*pd.l_rc/ps.M;
pp.a9 = pd.G/ps.m/pd.rho_w/pd.L;
pp.a10 = ps.tau*ps.u/ps.m/pd.rho_w/pd.L;
pp.a11 = ps.q*ps.phi/ps.x/ps.m/pd.rho_w/pd.L;
pp.a12 = ps.q*ps.phi/ps.x/ps.m/pd.rho_w/pd.L*pd.gamma*pd.rho_w*pd.c;
pp.a13 = 1;
pp.a14 = ps.u*pd.h_r/pd.l_r*ps.t/ps.h;
pp.a15 = ps.h/pd.h_r;
pp.a16 = pd.Ahat*ps.N^pd.n*ps.t;
pp.a17 = ps.S*ps.x/ps.t/ps.Q;
pp.a18 = 1;
pp.a19 = ps.M*ps.x/ps.Q;
pp.a20 = 1;
pp.a21 = ps.kappa*ps.x/ps.Q;
pp.a22 = ps.h*ps.x/ps.t/ps.q;
pp.a23 = ps.x/ps.t/ps.q*ps.sigma*ps.phi/pd.rho_w/pd.g;
pp.a24 = 1;
pp.a25 = ps.m*ps.x/ps.q;
pp.a26 = ps.kappa*ps.x/ps.q/ps.W;
pp.a27 = pd.lambda*ps.N/ps.kappa;
pp.a28 = ps.Q/pd.K_c/ps.S^pd.alpha_c/(ps.phi/ps.x)^(1/2);
pp.a29 = pd.K*ps.h^pd.alpha*(ps.phi/ps.x)/ps.q;
pp.a30 = pd.rho_w*ps.m*ps.t/pd.rho_i/ps.h;
pp.a31 = ps.x/ps.t/ps.Q*pd.sigma_c*ps.W*ps.phi/pd.rho_w/pd.g;
pp.a32 = pd.rho_w/pd.rho_i;
pp.a33 = ps.S/ps.t/ps.Q*ps.phi/pd.rho_w/pd.g;
pp.a34 = ps.M_in*ps.x/ps.Q;
pp.a35 = ps.m_in*ps.x/ps.q;
pp.a36 = 1;
pp.a37 = pd.p_ff/ps.phi;
pp.a38 = ps.S*ps.x/ps.t/ps.Q*pd.beta_c*ps.phi;
pp.a39 = pd.lambda2*pd.K*ps.h^pd.alpha*ps.phi;
pp.a40 = pd.p_u/ps.phi;
pp.a41 = pd.E_u/ps.h*ps.phi;
pp.a42 = ps.h*ps.x/ps.t/ps.q;
pp.a43 = ps.h*ps.x/ps.t/ps.q*pd.beta*ps.phi;
pp.a44 = pd.Abar^(-1/pp.n)*(ps.u/ps.x).^(1/pp.n)*ps.s/ps.x/ps.tau;  
pp.a45 = pd.rho_i*pd.g*ps.s^2/ps.x/ps.tau;
pp.a46 = 1;
pp.a47 = pd.Re_reg*ps.u/ps.t/ps.tau;  % artificially enhanced Reynold's number [ numerical convenience ]
pp.a48 = pd.C_slide*ps.u^(pp.q_slide).*ps.N^(pp.p_slide)/ps.tau;
pp.a49 = pd.eps_reg/(ps.u/ps.x);
pp.a50 = pd.N_min/ps.phi; 
pp.a51 = pd.phi_m/ps.phi;
pp.a52 = (1-pd.rho_i/pd.rho_o);
pp.a53 = pd.Abar^(-1/pp.n)*(ps.u/ps.W).^(1/pp.n)*ps.s/ps.W/ps.tau;

% TURN OFF CERTAIN TERMS IN EQUATIONS
if ~oo.include_a4, pp.a4 = 0; end 
if ~oo.include_a8, pp.a8 = 0; end 
if ~oo.include_a11, pp.a11 = 0; end     
if ~oo.include_a12, pp.a12 = 0; end  
if ~oo.include_a44, pp.a44 = 0; end  
if ~oo.include_a53, pp.a53 = 0; end  

% DISCRETIZATION
dd = wyvis_discretize(ad.x/ps.x);

% DIMENSIONLESS GEOMETRY AND INPUTS
aa.x = ad.x/ps.x;
aa.Z_b = ad.Z_b/ps.s;
aa.Z_s = ad.Z_s/ps.s;
aa.W = ad.W/ps.W;
aa.U_b = ad.U_b/ps.u;
aa.x_m = ad.x_m/ps.x;
aa.Q_m = @(t) ad.Q_m(ps.t*t)/ps.Q;
aa.S_m = ad.S_m/ps.S;
aa.M_in = @(t) ad.M_in(ps.t*t)/ps.M_in;
aa.m_in = @(t) ad.m_in(ps.t*t)/ps.m_in;

% DIMENSIONLESS DERIVED OR PRESCRIBED FIELDS
aa.xim = wyvis_nearest_gridpoint(aa.x_m,dd.x);
aa.phi_s = aa.Z_s + (pp.a32-1)*aa.Z_b;
aa.phi_b = pp.a32*aa.Z_b;
aa.phi_ext = max(aa.phi_b,pp.a51);
aa.phi_ext(dd.xin) = aa.phi_s(dd.xin);
aa.tau_d = -(aa.Z_s-aa.Z_b).*(dd.ddxx*aa.Z_s);
aa.H = aa.Z_s-aa.Z_b;
aa.sigma = ones(dd.I,1);
aa.sigma_c = ones(dd.I,1);
aa.h_ext = NaN*ones(dd.I,1);
aa.S_ext = NaN*ones(dd.I,1);
aa.Q_ext = 0*aa.x;
aa.q_ext = 0*aa.x;
aa.u_ext = aa.U_b;

% UNSCALED DERIVED OR PRESCRIBED FIELDS
ad.xim = aa.xim;
ad.phi_s = ps.phi*aa.phi_s;
ad.phi_b = ps.phi*aa.phi_b;
ad.phi_ext = ps.phi*aa.phi_ext;
ad.tau_d = ps.tau*aa.tau_d;
ad.H = ps.s*aa.H;
ad.U_b = ps.u*aa.U_b;
ad.sigma = ps.sigma*aa.sigma;
ad.sigma_c = ps.sigma_c*aa.sigma_c;

end

function ad = wyvis_unscale_aa(aa,ps)
% unscale prescribed field structure aa
    ad.x = ps.x*aa.x;
    ad.Z_b = ps.s*aa.Z_b;
    ad.Z_s = ps.s*aa.Z_s;
    ad.W = ps.W*aa.W;
    ad.U_b = ps.u*aa.U_b;
    ad.x_m = ps.x*aa.x_m;
    ad.Q_m = @(t) ps.Q*aa.Q_m(t/ps.t);
    ad.S_m = ps.S*aa.S_m;
    ad.M_in = @(t) ps.M_in*aa.M_in(t/ps.t);
    ad.m_in = @(t) ps.m_in*aa.m_in(t/ps.t);
    
    ad.xim = aa.xim;
    ad.phi_s = ps.phi*aa.phi_s;
    ad.phi_b = ps.phi*aa.phi_b;
    ad.phi_ext = ps.phi*aa.phi_ext;
    ad.tau_d = ps.tau*aa.tau_d;
    ad.sigma = ps.sigma*aa.sigma;
    ad.Q_ext = ps.Q*aa.Q_ext;
    ad.q_ext = ps.q*aa.q_ext;
end

function ud = wyvis_unscale_uu(uu,ps)
% unscale solution structure uu
    ud = uu;
    for i = 1:length(uu)
        ud(i).S = ps.S*uu(i).S;
        ud(i).phi_c = ps.phi*uu(i).phi_c;
        ud(i).h = ps.h*uu(i).h;
        ud(i).phi = ps.phi*uu(i).phi;
        ud(i).u = ps.u*uu(i).u;
        ud(i).Q = ps.Q*uu(i).Q;
        ud(i).q = ps.q*uu(i).q;
        ud(i).N_c = ps.N*uu(i).N_c;
        ud(i).N = ps.N*uu(i).N;
        ud(i).M = ps.M*uu(i).M;
        ud(i).m = ps.m*uu(i).m;
        ud(i).kappa = ps.kappa*uu(i).kappa;
        ud(i).Sigma = ps.Sigma*uu(i).Sigma;
        ud(i).Sigma_c = ps.Sigma_c*uu(i).Sigma;
        ud(i).V_in = ps.V_in*uu(i).V_in;
        ud(i).Sw = ps.S*uu(i).Sw;
        ud(i).hw = ps.h*uu(i).hw;
        ud(i).hu = ps.h*uu(i).hu;
    end
end

function dd = wyvis_discretize(x)
% produce struct of discretization details
    I = length(x);
    dd.I = I;
    dd.x = x;
    dd.xx = [dd.x(1); (dd.x(1:I-1)+dd.x(2:I))/2];
    dd.dx = [dd.xx(2:I); dd.x(end)]-dd.xx;
    dd.dxx = dd.x(1:I)-[0; dd.x(1:I-1)];
    
    dd.ddx = sparse([1:I 1:I-1],[1:I 1+(1:I-1)],[-dd.dx.^(-1); dd.dx(1:I-1).^(-1)],I,I); % divergence operator
    dd.ddxx = sparse([2:I 1:I],[(2:I)-1 (1:I)],[-dd.dxx(2:I).^(-1); dd.dxx.^(-1)],I,I); dd.ddxx(isinf(dd.ddxx)) = 0; % gradient operator
    dd.ddxx = sparse([1:I 1:I],[1 (2:I)-1 2 (2:I)],[-dd.dxx([2 2:I]).^(-1); dd.dxx([2 2:I]).^(-1)],I,I); % gradient operator [ with extrapolated value on left edge ]
    dd.avx = sparse([1:I 1:I-1],[1:I 1+(1:I-1)],[0.5*ones(I-1,1); 1; 0.5*ones(I-1,1)],I,I); % averaging operator on nodes
    dd.avxx = sparse([2:I 1:I],[(2:I)-1 (1:I)],[0.5*ones(I-1,1); 1; 0.5*ones(I-1,1)],I,I); % averaging operator on edges
    
    dd.xext = I;    % Dirichlet nodes
    dd.xxext = 1;   % Neumann nodes
    dd.xin = setdiff(1:I,dd.xext); % intetior nodes
    dd.xxin = setdiff(1:I,dd.xxext); % interior edges
    
    tmp = dd.avx(:,dd.xxin)*ones(length(dd.xxin),1); 
    dd.avxin = spdiags(tmp.^(-1),0,I,I)*dd.avx; dd.avxin(:,dd.xxext) = 0; % averaging operator on nodes using interior edges
    tmp = dd.avxx(:,dd.xin)*ones(length(dd.xin),1); 
    dd.avxxin = spdiags(tmp.^(-1),0,I,I)*dd.avxx; dd.avxxin(:,dd.xext) = 0; % averaging operator on edges using interior nodes
end

function xi = wyvis_nearest_gridpoint(x,xd)
% indices of nearest grid points to x
    xi = NaN*x;
    for i = 1:length(x)
        [~,tmp] = min((x(i)-xd).^2);
        xi(i) = tmp;
    end
end

function mass = wyvis_mass(Y,aa,dd,pp)
% mass matrix for equations
    xin = dd.xin;
    Iin = length(xin);
    S_in = 0*dd.x; S_in(aa.xim) = aa.S_m;  % assign moulins to S_in vector
    S = aa.S_ext; S(xin,:) = Y(1:Iin,:);
    h = aa.h_ext; h(xin,:) = Y(2*Iin+(1:Iin),:);
    phi = aa.phi_ext; phi(xin,:) = Y(3*Iin+(1:Iin)); 
    dhu = wyvis_duplift(phi-aa.phi_s,pp); 
    
    mass = sparse([1:Iin                Iin+(1:Iin)          Iin+(1:Iin)                                                                                2*Iin+(1:Iin)           3*Iin+(1:Iin)        3*Iin+(1:Iin)      4*Iin+(1:Iin)], ...
                   [1:Iin                1:Iin                Iin+(1:Iin)                                                                                2*Iin+(1:Iin)           2*Iin+(1:Iin)        3*Iin+(1:Iin)      4*Iin+(1:Iin)], ...
                    [pp.a1*ones(1,Iin)    pp.a17*ones(1,Iin)   pp.a31*aa.sigma_c(xin)'.*aa.W(xin)'+pp.a38*S(xin)'+pp.a33*S_in(xin)'.*dd.dx(xin)'.^(-1)    pp.a13*ones(1,Iin)      pp.a22*aa.W(xin)'    pp.a23*aa.sigma(xin)'.*aa.W(xin)'+pp.a42*dhu(xin)'.*aa.W(xin)'+pp.a43*h(xin)'.*aa.W(xin)'        pp.a47*aa.W(xin)'], ...
                     5*Iin , 5*Iin );
%     if pp.a31>0, 
%          pause;
%     end

    % enlarge
end

function mvpattern = wyvis_mvpattern(dd)
% sparsity pattern for mass matrix
    xin = dd.xin;
    Iin = length(xin);
    tmp = speye(Iin);
    mvpattern = [tmp 0*tmp 0*tmp 0*tmp 0*tmp; tmp tmp 0*tmp 0*tmp 0*tmp; 0*tmp 0*tmp tmp 0*tmp 0*tmp; 0*tmp 0*tmp tmp tmp 0*tmp; 0*tmp 0*tmp 0*tmp 0*tmp tmp];
end

function jpattern = wyvis_jpattern(dd)
% sparsity pattern for jacobian of equations (conservative pattern -
% jacobian is likely more sparse than this)
    xin = dd.xin;
    Iin = length(xin);
    tmp = sparse([1:Iin 1:Iin-1 2:Iin],[1:Iin 2:Iin 1:Iin-1],ones(1,3*Iin-2),Iin,Iin);
    jpattern = [tmp tmp tmp tmp tmp; tmp tmp tmp tmp tmp; tmp tmp tmp tmp tmp; tmp tmp tmp tmp tmp; tmp tmp tmp tmp tmp];
end

function F = wyvis_F(t,Y)
% right hand side of equations 
    [S,~,h,~,u,Q,N_c,M,q,N,m,kappa,~,~,~,~,~,~,tau_b] = wyvis_UP_full(Y,aa,dd,pp);
    
    ext = ones(1,size(Y,2));
    Q_in = 0*dd.x; Q_in(aa.xim) = aa.Q_m(t);   % assign moulins to Q_in vector
    M_in = aa.M_in(t);
    m_in = aa.m_in(t);
    
    F1 = pp.a2*M - pp.a3*S.*abs(N_c).^(pp.n-1).*N_c + pp.a4*(dd.avx*u).*max(1-pp.a5*S,0);
    F2 = -pp.a18*(dd.ddx*Q) + pp.a19*M + pp.a20*(Q_in.*dd.dx.^(-1))*ext + pp.a21*kappa + pp.a34*(M_in*ext);
    F3 = pp.a14*(dd.avx*u).*max(1-pp.a15*h,0) - pp.a16*h.*abs(N).^(pp.n-1).*N + pp.a30*m;
    F4 = -pp.a24*(dd.ddx*((aa.W*ext).*q)) + pp.a25*(aa.W*ext).*m - pp.a26*kappa + pp.a35*(aa.W*ext).*(m_in*ext);
    F5 = dd.ddxx*( (aa.W*ext).*[ pp.a44*2*(aa.H(1:dd.I-1)*ext).*((dd.ddx(1:dd.I-1,:)*u).^2+pp.a49^2).^((1/pp.n-1)/2).*(dd.ddx(1:dd.I-1,:)*u); pp.a52*pp.a45*1/2*aa.H(dd.I)^2*ext] ) + pp.a45*(aa.tau_d*ext).*(aa.W*ext) - pp.a46*tau_b.*(aa.W*ext) - pp.a53*((pp.n+2)/2)^(1/pp.n).*u.^(1/pp.n)./aa.W.^(1/pp.n).*aa.H;
    if ~oo.include_ice_velocity, F5 = 0*F5; end % ice velocity derivatives zero if no ice velocity
    F = [F1(dd.xin,:); F2(dd.xin,:); F3(dd.xin,:); F4(dd.xin,:); F5(dd.xxin,:)];
end

function [S,phi_c,h,phi,u] = wyvis_initial(t,aa,dd)
% estimate for initial condition using potential phi_ext and prescribed inputs
    Smin = 1e-3;
    Q_in = 0*dd.x; Q_in(aa.xim) = aa.Q_m(t);   % assign moulins to Q_in vector
    M_in = aa.M_in(t);
    m_in = aa.m_in(t);        
    m = pp.a9 + pp.a10*aa.tau_d.*aa.U_b + 0;    % ignore q dependent terms in m
    M = 0*Q_in;     % ignore M
    m = dd.avxin*m;
    M = dd.avxin*M;
    Q = cumsum(pp.a19*M+pp.a20*Q_in+pp.a34*(M_in.*dd.dx)); % integrate inputs
    phi_c = aa.phi_ext;
    dphi_c = dd.ddxx*phi_c;
    S = abs(-Q/pp.a28./((dphi_c.^2+eps^2).^(1/2*(1/2-1)).*dphi_c)).^(1/pp.alpha_c);
    S = max(S,Smin);
    q = cumsum(pp.a25*(m.*dd.dx)+pp.a35*(m_in.*dd.dx)); % integrate inputs
    phi = aa.phi_ext;
    dphi = dd.ddxx*phi;
    h = abs(-q/pp.a29./(dphi)).^(1/pp.alpha);
    h = min(h,1/pp.a15);
    u = aa.u_ext;
end

function Y = wyvis_P(S,phi_c,h,phi,u)
% pack together variables into solution vector Y
    Y = [S(dd.xin); phi_c(dd.xin); h(dd.xin); phi(dd.xin); u(dd.xxin)];
end

function [S,phi_c,h,phi,u] = wyvis_UP(Y,aa,dd)
% unpack variables from solution vector Y
    xin = dd.xin;
    Iin = length(xin);
    ext = ones(1,size(Y,2));
    
    % inialize primary variables with boundary conditions
    S = aa.S_ext*ext;
    phi_c = aa.phi_ext*ext;
    h = aa.h_ext*ext;
    phi = aa.phi_ext*ext;
    u = aa.u_ext*ext;
    
    % assign interior values from solution vector
    tmp1 = 1; tmp2 = tmp1 + length(xin)-1;
    S(xin,:) = Y(tmp1:tmp2,:); tmp1 = tmp2+1; tmp2 = tmp1+Iin-1;
    phi_c(xin,:) = Y(tmp1:tmp2,:); tmp1 = tmp2+1; tmp2 = tmp1+Iin-1;
    h(xin,:) = Y(tmp1:tmp2,:); tmp1 = tmp2+1; tmp2 = tmp1+Iin-1;
    phi(xin,:) = Y(tmp1:tmp2,:); tmp1 = tmp2+1; tmp2 = tmp1+Iin-1;
    u(xxin,:) = Y(tmp1:tmp2,:);

end

function [S,phi_c,h,phi,u,Q,N_c,M,q,N,m,kappa,Sigma,Sigma_c,V_in,Sw,hw,hu,tau_b] = wyvis_UP_full(Y,aa,dd,pp)
% unpack variables and derived functions from solution vector Y
    xin = dd.xin;
    xxin = dd.xxin;
    Iin = length(xin);
    ext = ones(1,size(Y,2));
    
    % inialize primary variables with boundary conditions
    S = aa.S_ext*ext;
    phi_c = aa.phi_ext*ext;
    h = aa.h_ext*ext;
    phi = aa.phi_ext*ext;
    u = aa.u_ext*ext;
    
    % assign interior values from solution vector
    tmp1 = 1; tmp2 = tmp1 + length(xin)-1;
    S(xin,:) = Y(tmp1:tmp2,:); tmp1 = tmp2+1; tmp2 = tmp1+Iin-1;
    phi_c(xin,:) = Y(tmp1:tmp2,:); tmp1 = tmp2+1; tmp2 = tmp1+Iin-1;
    h(xin,:) = Y(tmp1:tmp2,:); tmp1 = tmp2+1; tmp2 = tmp1+Iin-1;
    phi(xin,:) = Y(tmp1:tmp2,:); tmp1 = tmp2+1; tmp2 = tmp1+Iin-1;
    u(xxin,:) = Y(tmp1:tmp2,:);
    
    % additional h due to uplift
    if oo.include_uplift
        hu = wyvis_uplift(phi-(aa.phi_s*ext),pp);
    else
        hu = 0*h;
    end

    % water filled fraction [ abs since may otherwise cause an issue ]
    if oo.include_fill_factor
        hw = wyvis_fill_factor(phi-(aa.phi_b*ext),pp).*abs(h+hu); 
        Sw = wyvis_fill_factor(phi_c-(aa.phi_b*ext),pp).*abs(S); 
    else 
        hw = h;
        Sw = S;
    end
    
    % fluxes
    hav = dd.avxxin*hw;
    Sav = dd.avxxin*Sw;
    dphi = dd.ddxx*phi;
    dphi_c = dd.ddxx*phi_c;
    dphi_b = dd.ddxx*(aa.phi_b*ext);
    
    q = aa.q_ext*ext;
    Q = aa.Q_ext*ext;
    tmp = -pp.a29*hav.^pp.alpha.*dphi;
    q(xxin,:) = tmp(xxin,:);
    tmp = -pp.a28*Sav.^pp.alpha_c.*(dphi_c.^2+eps^2).^(1/2*(1/2-1)).*dphi_c;
    Q(xxin,:) = tmp(xxin,:);
    
    % effective pressure and exchange term
    N = aa.phi_s*ext - phi;
    N_c = aa.phi_s*ext - phi_c;
    kappa = pp.a27*(N_c-N) + pp.a39*hav.^(pp.alpha).*(N_c-N);
    
    % basal shear stress
    tau_b = wyvis_sliding(u,dd.avxx*(N),pp);
    
    % melting rates
    if oo.couple_frictional_heating
    m = pp.a9 + pp.a10*(tau_b.*u) ...
        - pp.a11*q.*dphi ...
        + pp.a12*q.*(dphi-dphi_b);
    else
    m = pp.a9 + pp.a10*(aa.tau_d.*aa.U_b)*ext ...
        - pp.a11*q.*dphi ...
        + pp.a12*q.*(dphi-dphi_b);
%     m = pp.a9 + pp.a10*dd.avx*(tau_b.*u) ...
%         - pp.a11*q.*dphi ...
%         + pp.a12*q.*(dphi-dphi_b);
    end
    M = - pp.a6*Q.*dphi_c ...
        + pp.a7*Q.*(dphi_c-dphi_b) ...
        + pp.a8*m ;
    
    m = dd.avxin*m;
    M = dd.avxin*M;
    
    % storage terms
    Sigma = (aa.sigma*ext).*(phi-aa.phi_b*ext);
    Sigma_c = (aa.sigma_c*ext).*(phi_c-aa.phi_b*ext);
    S_in = 0*dd.x; S_in(aa.xim) = aa.S_m;  % assign moulins to S_in vector
    V_in = (S_in*ext).*(phi_c-aa.phi_b*ext);
end

function f = wyvis_fill_factor(p_w,pp)
% fill factor as function of water pressure
    f = 1-pp.a36*(1-tanh((p_w-pp.a37)/pp.a37))/2;
end

function h = wyvis_uplift(p,pp)
% uplift as function of pressure excess
    h = pp.a41*pp.a40*log(1+exp(p/pp.a40));
end

function h = wyvis_duplift(p,pp)
% derivative of uplift with respect to pressure
    h = pp.a41*(1+exp(p/pp.a40)).^(-1).*exp(p/pp.a40);
end

function tau = wyvis_sliding(u,N,pp)
% sliding law
    tau = pp.a48*u.^pp.q_slide.*max(N,pp.a50).^pp.p_slide;
end

end