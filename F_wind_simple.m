function [F_Tauwind] = F_wind_simple(dxdt_hub,tindex)
%% Description
% This function calculates the forcing from wind, and can compute for
% steady wind or for unsteady time series.
% This function takes the V_10 time average, V_hub, a and b parameters, and
% the hub motion to calculate the relative velocity, corrected CT, and
% determine thrust.
% This function returns the forcing and torque for a given point in time
% for Region 1 operation.
%% Implementation
global rho_air A_r V_10 V_hub z_hub CT_0
% Find V_rel
V_rel = V_hub(tindex)-dxdt_hub;
% Calculate mean aero force using first order expansion of thrust around
% V_10
Fwind = 0.5*rho_air*A_r*CT_0*V_10^2 + (V_rel - V_10)*rho_air*A_r*V_rel*CT_0;
Tauwind = Fwind*z_hub;
% 7 DOF, assuming no yaw
F_Tauwind=[Fwind;0;0;0;Tauwind;0;0];    % 7th term for 0 controller input
% This needs to be updated to take torque into account for unsteady wind