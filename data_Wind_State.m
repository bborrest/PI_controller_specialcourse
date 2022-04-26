function [fHighCut,rho_air,TI,TL] = data_Wind_State
% Wind State
fHighCut = 0.5;             % cut-off frequency for Kaimal spectrum
rho_air = 1.22;             % air density [kg/m^3]
TI = 0.14;                  %*** Kaimal wind spectra turbulence intensity
TL = 340.2;                 %*** Kaimal wind spectra turbulenec length scale [m]