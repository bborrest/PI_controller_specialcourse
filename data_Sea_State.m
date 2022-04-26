function [Hs,Tp,gammaJS,df,rho_H20] = data_Sea_State(Tdur)
% Sea State
Hs = 2;                     % linear wave amplitude and significant wave height [m]
Tp = 6;                     % linear wave period and significant wave period [s]
gammaJS = 3.3;              %*** JONSWAP peak enhancement factor
df = 1/Tdur;                % JONSWAP frequency spectra time step
rho_H20 = 1025;             % salt water density [kg/m^3]