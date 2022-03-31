function [T,M] = iceHydrates_flowModel(hMH,totMH,tLength,d,T_m,T_t,eta_0,...
    rhoFunc,viscFunc,kappaFunc,JExp,R_m,grav,res)
% Authors: Evan Carnahan (evan.carnahan@utexas.edu), Steve Vance, Marc Hesse, 
%          Baptiste Jorneuax, Christophe Sotin
% Date: 3/30/2022

% calculate parameters values from inputs
fMH = totMH-hMH; %km or km^3/km^2 as area flux, over life time of convection
DT = T_m-T_t; % K
c_p = 2.108e3; %J kg^-1 K^-1

%% Characteristic scale and dimensionless property relationships
% nondimensional thermal conductivity
k_c = kappaFunc(1,0);
kappaPFunc = @(nonT,MCon) kappaFunc(nonT,MCon)./k_c;

% nondimensional density
rho_c = rhoFunc(1,0);
rho_t = rhoFunc(0,0);
DRho = rho_c - rho_t;
rhoPFunc = @(nonT,MCon) (rhoFunc(nonT,MCon) - rho_c)./DRho; 

% nondimensional viscosity
viscPFunc = @(nonT,MCon) viscFunc(nonT,MCon,JExp)./eta_0;
eta_max = 6.5e14; %Pa s viscosity of maximum tidal dissipation

D_T = k_c/rho_c/c_p; % thermal diffusion, m^2/s
t_c = d^2/D_T; % time scale, s

%% Dimensionless numbers
pi_tide = R_m*t_c/(rho_c*c_p*DT); % dimensionless tidal heating
Ra = DRho*grav*d^3/(eta_0*D_T); % Rayliegh number calc

%% Grid construction
anis = 2;
lam = 2;
Gridp.xmin = 0; Gridp.xmax = 1*lam; Gridp.Nx = res*lam;
Gridp.ymin = 0; Gridp.ymax = 1; Gridp.Ny = res*anis;

Grid  = build_stokes_grid(Gridp);
Zp = zeros(Grid.p.N);
Ip = speye(Grid.p.N);

tempPert = cos(2*pi/lam*Grid.p.xc);
fullPert = repmat(tempPert',Grid.p.Ny,1);
[X,Y] = meshgrid(Grid.p.xc,Grid.p.yc);

[D,Edot,Dp,Gp,I,Gyy]=build_stokes_ops(Grid); 

linInds = find(Gyy > 0);
[row,~] = ind2sub(size(Gyy),linInds);
nyyVecT = zeros(Grid.y.Nfy,1);
nyyVecM = zeros(Grid.y.Nfy,1);

%% Simple intital condition - linear temperature profile
Tplot = 1- Y + 0.001*fullPert;
T = Tplot(:);

%% Build temperature boundary condition
T1 = 1;
T0 = 0;
Param.T.dof_dir = [Grid.p.dof_ymin;Grid.p.dof_ymax];
Param.T.dof_f_dir = [Grid.p.dof_f_ymin;Grid.p.dof_f_ymax];
Param.T.g = [T1*ones(length(Grid.p.dof_ymin),1); T0*ones(size(Grid.p.dof_ymax))];
Param.T.dof_neu = [Grid.p.dof_xmin;Grid.p.dof_xmax];
Param.T.dof_f_neu = [Grid.p.dof_f_xmin;Grid.p.dof_f_xmax];
Param.T.qb = 0*[Grid.p.dof_f_xmin;Grid.p.dof_f_xmax];
[BT,NT,~] = build_bnd(Param.T,Grid.p,Ip);
Param.T.dof_out = [Grid.p.dof_ymin];

%% Methane clathare initial condition
Mplot = zeros(Grid.p.Ny,Grid.p.Nx);
Mplot(flip(Grid.p.yc < hMH/d),:) = 1;
M = Mplot(:);

%% Build methane flux boundary condition 
M0 = M(end); % top
volMHFlux = fMH/tLength/3.154e7; %m/s/yr or m km^2/km^2/s as area flux, over life time of convection
volMHFlux_non = volMHFlux*t_c/d; %non-dimensional volume flux

Param.M.dof_dir = [Grid.p.dof_ymax];
Param.M.dof_f_dir = [Grid.p.dof_f_ymax];
Param.M.g = [M0*ones(size(Grid.p.dof_ymax))];

Param.M.dof_neu = [Grid.p.dof_xmin;Grid.p.dof_xmax; Grid.p.dof_ymin];
Param.M.dof_f_neu = [Grid.p.dof_f_xmin;Grid.p.dof_f_xmax;Grid.p.dof_f_ymin];
Param.M.qb = [0*Grid.p.dof_f_xmin; 0*Grid.p.dof_f_xmax;volMHFlux_non*ones(size(Grid.p.dof_f_ymin))];
[BM,NM,fn_M] = build_bnd(Param.M,Grid.p,Ip);
Param.M.dof_out = [Grid.p.dof_ymin];

%% Build boundary conditions
Param.dof_dir =  [...
                  Grid.x.dof_xmax;...           %set x_max x-vel
                  Grid.x.dof_xmin;...           %set x_min x-vel
                  Grid.x.N+Grid.y.dof_ymin;...  %set y_min y-vel
                  Grid.x.N+Grid.y.dof_ymax;...  %set y_max y-vel
                  Grid.p.Nf+1];                 %set pressure
                         
              
Param.g = 0*Param.dof_dir;
Param.g(end) = 0;
B = I([Param.dof_dir],:);
N = I;
N(:,[Param.dof_dir]) = [];
theta = 0.5; % set Crank-Nichelson theta value

%% Dissociation plot values
p0 = 1.5e5; % Pa from Sohl, 1995
rho_c = 917; % kg/m^3
grav = 1.352; % m/s^2

p = p0 + rho_c*grav*(Grid.p.yc)*d; %kg s^-2 m^-1 or Pa
TDis_func = @(p) 264.395+21.105*log10(p/1e6)-0.0424805 * (log10(p/1e6)).^2; %calculated in M Pa
TDis = TDis_func(p);
% reduce dissociation temp if impurities present
TDis_m20 = TDis - 20; %shift based off of Prasad, 2019

%% Time evolution of system
dtVec = [];
i = 1;
warning off % matrix is illconditioned

while ~(sum(dtVec) > tLength || any(M > 1.01))
    Tplot= reshape(T,Grid.p.Ny,Grid.p.Nx);
    Mplot= reshape(M,Grid.p.Ny,Grid.p.Nx);

    % Bouyancy forcing
    rhoPrime = rhoPFunc(T,M);
    rhoPlot = reshape(rhoPrime,Grid.p.Ny,Grid.p.Nx);
    rhoFace = diag(comp_mean(rhoPlot,1,1,Grid.p));
    rhoY = rhoFace(Grid.p.Nfx+1:Grid.p.Nf);
    fsVec = Ra*rhoY;
    fs = [zeros(Grid.p.Nfx,1); fsVec; zeros(Grid.p.N,1)];
    
    % Temperature-dependent viscosity
    nxxVec = zeros(Grid.x.Nfx,1);
    nxxVec(Grid.x.Ny+Grid.p.dof) = Tplot;
    nyyVecT(row) = Tplot;
    ncVec = comp_mean_corners(Tplot,-1,Grid.p);
    tempVec = [nxxVec;nyyVecT;ncVec];
    
    nxxVecM = zeros(Grid.x.Nfx,1);
    nxxVecM(Grid.x.Ny+Grid.p.dof) = Mplot;
    nyyVecM(row) = Mplot;
    ncVecM = comp_mean_corners(Mplot,-1,Grid.p);
    MVec = [nxxVecM;nyyVecM;ncVecM];
    viscVec = viscPFunc(tempVec,MVec);
    viscMat = spdiags(viscVec,0,length(viscVec),length(viscVec));

    % Form strain matrix
    tau = D*2*viscMat*Edot;

    % Form linear operator
    L = [tau,Gp;
        Dp,Zp];

    % Solve linear operator
    u = solve_lbvp(L,fs,B,Param.g,N);
    
    % Extract velocities
    vx = u(1:Grid.p.Nfx);
    vy = u(Grid.p.Nfx+1:(Grid.p.Nfx+Grid.p.Nfy));
    vm = [vx;vy];
    vmax= max((abs(vm)));
    dt = min([0.5*Grid.p.dx^2, Grid.p.dx/vmax,0.5*Grid.p.dy^2, Grid.p.dy/vmax])*0.8;
    
    % Tidal forcing
    tideVisc = viscPFunc(T,M);
    R_prime = 2*((eta_0/eta_max)*tideVisc)./(1+((eta_0/eta_max)*tideVisc).^2);
    R_tide = R_prime*pi_tide;
    
    % Thermal conductivity matrix
    thermCond = kappaPFunc(T,M);
    thermCondPlot = reshape(thermCond,Grid.p.Ny,Grid.p.Nx);
    thermCondFace = comp_mean(thermCondPlot,1,1,Grid.p);
    
    % Temperature advection-diffusion equation
    TOld = T; % create for heat flux calculation
    AT = build_adv_op(vm,T,dt,Gp,Grid.p,Param.T,'mc');
    L_T_E = Ip - (1-theta)*dt*Dp*(AT-thermCondFace*Gp);
    L_T_I = Ip + (theta)*dt*Dp*(AT-thermCondFace*Gp);
    RHS_T = L_T_E*T + dt*R_tide;
    T = solve_lbvp(L_T_I,RHS_T,BT,Param.T.g,NT);
    
    % Methane clathrate advection equation
    AM = build_adv_op(vm,M,dt,Gp,Grid.p,Param.M,'gd');
    L_M_E = Ip - (1-theta)*dt*Dp*(AM);
    L_M_I = Ip + (theta)*dt*Dp*(AM);
    RHS_M = L_M_E*M + fn_M*dt;
    M = solve_lbvp(L_M_I,RHS_M,BM,Param.M.g,NM);
        
    %% Heat flux calculation
    res_trans = @(T,TOld) (L_T_I*T - L_T_E*TOld)/dt - R_tide;
    flux_trans = @(T,TOld) (AT - Gp)*(theta*TOld+(1-theta)*T);
    q_non = comp_flux_gen(flux_trans,res_trans,T,Grid.p,Param.T,TOld);
    
    % redimensionalizing heat fluxes, Watts
    surfHeatF_dim = q_non(Grid.p.dof_f_ymax) * (k_c * DT / d);%W/m^2
    botHeatF_dim = q_non(Grid.p.dof_f_ymin)  * (k_c * DT / d);%W/m^2
    tidePiDim = R_prime * R_m; % W/m^3
    strDim = (T - TOld)/dt * (rho_c * c_p * DT / t_c); %W/m^3

    % total energy fluxes over time
    totBotHeatFlux_dim(i) = sum(botHeatF_dim .* Grid.p.A(Grid.p.dof_f_ymin)) * d * 1; % W
    totSurfHeatFlux_dim(i) = sum(surfHeatF_dim .* Grid.p.A(Grid.p.dof_f_ymax)) * d * 1; % W
    totTide_dim(i) = sum(tidePiDim.* Grid.p.V) * d * d * 1; % W
    totStr_dim(i) = sum(strDim.*Grid.p.V) * d * d * 1; % W
    totFlux_dim(i) = totTide_dim(i) + totBotHeatFlux_dim(i) - totSurfHeatFlux_dim(i)- totStr_dim(i);
    
    % keeping time
    dtVec(i) = dt*t_c/3.154e7; % yr

    %% plotting
    i
    if mod(i,100) == 0 || i == 1 || (sum(dtVec) > tLength || any(M > 1.01))
        figure(5);
        set(gcf, 'Position', [50 50 1500 900])
        [PSI,~,~] = comp_streamfun(vm,Grid.p);
        TMean = mean(Tplot,2);
        [TMax,~] = max(Tplot,[],2);
        
        subplot(3,3,1)
        cla;
        hold on
        axis equal
        contourf(X*d/1e3,Y*d/1e3,Tplot*DT+T_t,40,'linestyle','none'),view(2),hold on
        contour(X*d/1e3,Y*d/1e3,Tplot*DT+T_t,'r','LevelList',T_m)
        c = colorbar('NorthOutside');
        cStrVal = linspace(min(PSI(:)),max(PSI(:)),10);
        contour(Grid.p.xf*d/1e3,Grid.p.yf*d/1e3,PSI,'k','LevelList',cStrVal);
        caxis([100 273]);
        c.Label.String = 'Temperature, K';
        xlabel('x-dir, km')
        ylabel('z-dir, km')
        
        subplot(3,3,2)
        cla;
        axis equal
        cenVisc = reshape(eta_0*viscVec(Grid.p.Ny+1:Grid.p.Ny+Grid.p.N),Grid.p.Ny,Grid.p.Nx);
        contourf(X*d/1e3,Y*d/1e3,log10(cenVisc),40,'linestyle','none'),view(2),hold on
        c = colorbar('NorthOutside');
        c.Label.String = 'Log. Viscosity, Pa s';
        xlabel('x-dir, km')
        ylabel('z-dir, km')
        
        subplot(3,3,3)
        cla;
        hold on
        axis equal 
        contourf(X*d/1e3,Y*d/1e3,reshape(R_tide*rho_c*c_p*DT/t_c,Grid.p.Ny,Grid.p.Nx),40,'linestyle','none'),view(2),hold on
        c = colorbar('NorthOutside');
        xlabel('x-dir, km')
        ylabel('z-dir, km')
        c.Label.String = 'Tidal energy, W/m^3';
       
        subplot(3,3,4)
        cla;
        hold on
        axis equal
        contourf(X*d/1e3,Y*d/1e3,reshape(thermCondPlot*k_c,Grid.p.Ny,Grid.p.Nx),40,'linestyle','none'),view(2),hold on
        c = colorbar('NorthOutside');
        xlabel('x-dir, km')
        ylabel('z-dir, km')
        c.Label.String = 'Thermal conductivity, W m^{-1} K^{-1}';
        
        subplot(3,3,5)
        cla;
        hold on
        axis equal
        contourf(X*d/1e3,Y*d/1e3,Mplot,40,'linestyle','none'),view(2),hold on
        c = colorbar('NorthOutside');
        c.Label.String = 'Methane clathrate hydrate volume fraction, 1';
        xlabel('x-dir, km')
        ylabel('z-dir, km')
        
        subplot(3,3,6)
        cla;
        hold on
        axis equal
        contourf(X*d/1e3,Y*d/1e3,reshape(rhoPlot*DRho+rho_c,Grid.p.Ny,Grid.p.Nx),40,'linestyle','none'),view(2),hold on
        c = colorbar('NorthOutside');
        xlabel('x-dir, km')
        ylabel('z-dir, km')
        c.Label.String = 'Density, kg/m^3';
        
        subplot(3,3,7)
        cla;
        hold on
        plot(TDis,flip(Grid.p.yc*d/1e3),'-r');
        plot(TDis_m20,flip(Grid.p.yc*d/1e3),'-k');
        plot(TMean*DT+T_t,Grid.p.yc*d/1e3,'-b');
        plot(TMax*DT+T_t,Grid.p.yc*d/1e3,'-g');
        legend('Disassociation temp','Dis temp w/ impurities (Prasad, 2019)','Mean temperature','Max temperature','Location','SouthWest')
        xlabel('Temperature, K');
        ylabel('z-dir, km');
        
        subplot(3,3,8)
        cla;
        hold on
        plot(cumsum(dtVec)/1e6,totBotHeatFlux_dim/(2*d)*1e3)
        plot(cumsum(dtVec)/1e6,-totSurfHeatFlux_dim/(2*d)*1e3)
        plot(cumsum(dtVec)/1e6,-totStr_dim/(2*d)*1e3)
        plot(cumsum(dtVec)/1e6,totTide_dim/(2*d)*1e3)
        plot(cumsum(dtVec)/1e6,totFlux_dim/(2*d)*1e3)
        ylabel('Heat flux, mW/m$^2$')
        legend('Bottom flux','Surface flux','Stored','Tidal','Total','Location','NorthWest');
        xlabel('MYrs');
        ylim([-30,30]);
    end
    i = i + 1;
end
return



