function k = ice_thermal_conductivity(T)
% Best fit to the measurents of Jakob & Erk (1929) and Dean & Timmerhaus
% (1963).
%
% Inputs:
% Temperature   (K)
%
% Outputs:
% Thermal Conductivity    (W/m/K)
%
% Source:
% Rabin (2000)
%
% Range of Validity:
% T > 80 K

%% Column Vectors
flip = false;
if ~iscolumn(T)
    T = T';
    flip = true;
end

%% Convert to Kelvin
if min(T)<0
    T = T+273.15;
end

%% Thermal Conductivity
a = 2135; % W/m
b = 1.235;
k = a./T.^b; % W/m/K

if flip
    k = k';
end
end