function cp = ice_heat_capacity(T,p)
% Specific isobaric heat capacity of ice as a function of temperature and
% pressure. Equations from the Gibbs energy equation of state.
%
% Inputs:
% Temperature   (K)
% Pressure      (Pa)
%
% Outputs:
% Heat Capacity (J/kg/K)
%
% Source: 
% IAPWS R10-06(2009)
% http://www.iapws.org/relguide/Ice-Rev2009.pdf
%
% Range of Validity:
% 0 K <= T <= 273.16 K
% 0 < p <= 210 MPa

%% Column Vectors
flip = false;
if ~iscolumn(T)
    T = T';
    flip = true;
end

if ~iscolumn(p)
    p = p';
    flip = true;
end

%% Convert to Kelvin
if min(T)<0
    T = T+273.15;
end

%% Coefficients of the equation of state (Gibbs potential function)
t1 = 0.368017112855051e-1+1i*0.510878114959572e-1;
t2 = 0.337315741065416+1i*0.335449415919309;
tk = [t1 t2];

r1 = 0.447050716285388e2+1i*0.656876847463481e2; % J/kg/K
r20 = -0.725974574329220e2+1i*-0.781008427112870e2; % J/kg/K
r21 = -0.557107698030123e-4+1i*0.464578634580806e-4; % J/kg/K
r22 = 0.234801409215913e-10+1i*-0.285651142904972e-10; % J/kg/K
r2k = [r20 r21 r22];

%% Constants
Tt = 273.16; % K
pt = 611.657; % Pa
pii0 = 101325/pt;

%% Dimensionless variables
tau = T/Tt;
pii = p/pt;

%% r2
r2 = 0;
for k = 0:2
    r2 = r2 + r2k(k+1)*(pii-pii0).^k;
end

%% Second derivative of Gibbs energy, gtt(T,p)
rk = [r1*ones(size(r2)) r2];
Re = 0;
for n = 1:2
    Re = Re + real(rk(:,n).*((1./(tk(n)-tau))+(1./(tk(n)+tau))-(2/tk(n))));
end
gtt = Re./Tt;

%% cp
cp = -T.*gtt; 

if flip
    cp = cp';
end
end