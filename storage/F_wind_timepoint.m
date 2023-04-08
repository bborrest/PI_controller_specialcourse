function [F_Tauwind,CT_hub] = F_wind_timepoint(dxdt_hub,tindex)
%% Description
% This function calculates the forcing from wind, and can compute for
% steady wind or for unsteady time series.
% This function takes the V_10 time average, V_hub, a and b parameters, and
% the hub motion to calculate the relative velocity, corrected CT, and
% determine thrust.
% This function returns the forcing and torque for a given point in time.
%% Implementation
global rho_air A_r V_rated CT_0 V_10 V_hub aCT bCT z_hub
% Find V_rel
V_rel = V_hub(tindex)-dxdt_hub;
% Set values for CT_hub
if V_rel <= V_rated
   CT_hub = CT_0;
else
   CT_hub = CT_0 *exp(-aCT*(V_rel-V_rated)^bCT);
end
% Set values for CT_10
if V_10 <= V_rated
   CT_10 = CT_0;
else
   CT_10 = CT_0 *exp(-aCT*(V_10-V_rated)^bCT);
end
% Calculate mean aero force
Fwind_m = 0.5*rho_air*A_r*CT_10*V_10^2;
% Calculate time-varying aero force
Fwind_t = 0.5*rho_air*A_r*CT_hub*V_rel^2;
% Calculate reduction factor
if V_10<V_rated
   f_red = 0.54;
elseif V_10
   f_red = 0.54 + 0.027*(V_10-V_rated);
end
Fwind = Fwind_m + f_red*(Fwind_t-Fwind_m);
Tauwind = Fwind*z_hub;
F_Tauwind=[Fwind;Tauwind];