% Wrapper to run simple mixed clathrate-ice shell convection demo
% Author: Evan Carnahan
% Date: 3/30/2022
clear all; close all

%% model parameters
% Titan-ish scenario
T_t = 95; % Surface temperature, K
grav = 1.354; % Gravity, m/s^2
R_m = 5e-8; % Maximum tidal heating, W/m^3
d = 100e3; % Ice shell thickness, km
hMH = 10e3; % Surface clathrate layer, m

% variable parameters
eta_0 = 1e15; % Basal viscosity, Pa s
J = 0.5; % Viscosity mixing exponenet

grRes = 35; % grid resolution in x-direction for simulations

%% Short-ish demo 
tLen = 4e7; % Simulation length, years
totMH = hMH + 15e3; % Total clathrate released, m 

%% Long demo - full planetary life simulation
% tLen = 4e9; % Simulation length, years
% totMH = 40e3; % Total clathrate released, m 

%% Pressure dependent melting temperature
par_Ih =[0.119539337e7,0.808183159e5,0.333826860e4,3,0.257500e2,0.103750e3,611.657e-9,273.16];
res=.01;
MS=0;
TmIh = 251.165-MS:res:273.16;
PmIh=par_Ih(7)*(1+ par_Ih(1).*(1-(TmIh./par_Ih(8)).^(par_Ih(4))) + par_Ih(2).*(1-(TmIh./par_Ih(8)).^(par_Ih(5))) + par_Ih(3).*(1-(TmIh./par_Ih(8)).^(par_Ih(6)))); % G Pa
p0 = 1.5e5; % Pa from Sohl, 1995
rho = 917; % kg/m^3
g = 1.352; % m/s^2
p = p0 + rho*g*d; %kg s^-2 m^-1 or Pa
pGPa = p/1e9;

T_m = interp1(PmIh,TmIh,pGPa);
DT = T_m - T_t;
%% Density
% ice
alpha = 1.6e-4; %1/K
rho_b = ice_density(T_m,5e6);
rhoIce_func = @(nonT) rho_b*(1-alpha*DT*(nonT-1));

% methane clathrate
rhoMH_func = @(nonT) Helgerud_sI(5,nonT*DT+T_t) + (Helgerud_sI(5,T_m) - Helgerud_sI(5,T_m));
rhoMix_func = @(nonT,MCon) MCon .* rhoMH_func(nonT) + (1-MCon) .* rhoIce_func(nonT);
rho_func = rhoMix_func;

%% Viscosity 
R = 8.314; % J K^-1 mol^-1
% ice
E_a = 50e3;% J mol^-1
Apar = E_a/R/T_m;
viscIce_func = @(nonT) eta_0*exp(Apar*(T_m./(DT.*nonT+T_t)-1));
viscPrimeIce_func = @(nonT) exp(Apar*(T_m./(DT.*nonT+T_t)-1));

% methane clathrate
E_aMH = 90e3;
T_mMH = T_m; 
AParMH = E_aMH/R/T_mMH;
viscMH_func = @(nonT) (20*eta_0)*exp(AParMH*(T_mMH./(nonT*DT+T_t)-1));

% mixture
viscMix_func = @(nonT,MCon,JExp) (MCon .* viscMH_func(nonT).^JExp + (1-MCon) .* viscIce_func(nonT).^JExp).^(1/JExp);
visc_func = viscMix_func;

%% Thermal conductivity
% ice
kappaIce_func = @(nonT) 612./(nonT*DT+T_t);

% methane clathrate 
kappaMH_func = @(nonT) 0.5*ones(size(nonT));

% mixture
kappaMix_func = @(nonT,MCon) MCon .* kappaMH_func(nonT) + (1-MCon) .* kappaIce_func(nonT);
kappa_func = kappaMix_func;

%% Model simulation
[T,M] = iceHydrates_flowModel(hMH,totMH,tLen,d,T_m,T_t,eta_0,...
    rho_func,visc_func,kappa_func,J,R_m,grav,grRes);
