function [F_Tauwind] = F_wind(q8,q12,tindex)
%% Description
% This function calculates the forcing from wind, and can compute for
% steady wind or for unsteady time series.
% This function takes the V_10 time average, V_hub, a and b parameters, and
% the hub motion to calculate the relative velocity, corrected CT, and
% determine thrust.
% This function returns the forcing and torque for a given point in time.
%% Implementation
global rho_air A_r V_10 CT_10 V_hub z_hub pTpV pQpU0 dtheta ki3 kp3
dtheta(tindex+1) = kp3*q12 + ki3*q8;
% if V_rel >= V_rated
    % the blade pitch and previous rotor speed are being stored in the state vector system
    % 1st order expansion of thrust force (F_op + CF_v + CF_om + CF_th) from Jiatian masters thesis
    %Fwind = 0.5*rho_air*A_r*V_10^2*CT_10 + pTpOm0*q14 - pTpV*(dxdt_hub) + pTpTh0*dtheta(tindex+1) + pTpV*(V_hub(tindex) - V_10);
    Fwind = 0.5*rho_air*A_r*V_10^2*CT_10 + pTpV*(V_hub(tindex) - V_10);
    Tauwind = Fwind*z_hub;
    % no yaw yet!!
    %F_Tauwind=[Fwind;0;0;0;Tauwind;0;0];    % 7th term for 0 controller input
    F_Tauwind=[Fwind;0;0;0;Tauwind;0;pQpU0*(V_hub(tindex) - V_10)];
% elseif V_rel < V_rated
%     % Calculate mean aero force using first order expansion of thrust around
%     % V_10
%     Fwind = 0.5*rho_air*A_r*CT_0*V_10^2 + (V_rel - V_10)*rho_air*A_r*V_rel*CT_0;
%     Tauwind = Fwind*z_hub;
%     % 7 DOF, assuming no yaw
%     F_Tauwind=[Fwind;0;0;0;Tauwind;0;0];    % 7th term for 0 controller input
%     % This needs to be updated to take torque into account for unsteady wind
% end