function [fHighCut,rho_air,TI,TL] = data_Wind_State(V_10,z_hub)
% Wind State
fHighCut = 0.5;             % cut-off frequency for Kaimal spectrum
rho_air = 1.22;             % air density [kg/m^3]
%% method for calculating longitudinal wind standard deviation, then converting to TI
% from page 35, 36 of IEC 61400-3-1_2019
Ac = 0.011;                  % Charnock's Constant
g = 9.81;                   % constant acceleration due to gravity [m/s^2]
k = 0.4;                    % Van karmen's constant
z0_guess = 0;               % roughness length scale initial guess [m]
z0 = 0.002;                 % initial set roughness length scale [m]
% find z0
while abs((z0_guess - z0)/z0)>10^-4
    z0_guess = z0;
    z0 = (Ac/g)*((k*V_10)/log(z_hub/z0_guess))^2;
end
I15 = 0.16;                 % Class Ib turbulence intensity
d = 4;                      % [m/s]
% Wind Standard Deviation
Sig1 = V_10/log(z_hub/z0) + 1.28*d*I15;
% Turbulence Intensity
TI = Sig1/V_10;
TL = 340.2;                 %*** Kaimal wind spectra turbulence length scale [m]