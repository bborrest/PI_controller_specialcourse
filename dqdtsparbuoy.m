function dq = dqdtsparbuoy(tode4,q,M,B,C,forcing)
%% Description
% This function takes a mass matrix, M (add added mass before), damping
% matrix, B, restoring matrix, C, and initial
% conditions, q
%
% The function then calculates the relevant forcing vector, depending on
% the forcing input case, forcing
%
% The function then returns a generalized motion, dq
%
% This function is developed for solving the surge and pitch motions of a
% floating, moored spar-buoy wind turbine system
%
% Future development will include making this function robust for various
% degree-of-freedom considerations
%
% As it is currently implemented, M,B,C are all 6x6 matrices; F is a 6x1
% forcing vector with Force in 1-3 and torque in 4-6; q is a 4x1 vector that is
% split into two 2x1 vectors for solving the system. The return, dq, is 4x1
%% Implementation
global t z_hub
% Find the proper time for forcing
[~ ,index]= min(abs(tode4-t));
% Calculate forcing
if forcing == 1
    % no external forcing
    % updated for 6 DOF + controller
    F = [0;0;0;0;0;0;0];
elseif forcing == 2
    % hydro forcing only
    % updated for 6 DOF, assuming 0 forcing applied to heave and yaw DOF
    F = hydroforcing(index,q(8),q(9),q(11),q(12));
elseif forcing == 3
    % wind only forcing (with PI controller)
    qhub = q(8)+z_hub*q(12);                            % dxdt of the hub
    F = F_wind_Region3(qhub,q(7),q(14),index);
elseif forcing == 4
    % hydro plus wind forcing with a PI controller implemented
    qhub = q(8)+z_hub*q(12);                            % dxdt of the hub
    [Windforcing] = F_wind_Region3(qhub,q(7),q(14),index);    % time-delay controller
    F = hydroforcing(index,q(7),q(8),q(10),q(11)) + Windforcing;                        % hydro plus steady wind forcing
end
dq8_14 = (M)\(F - B*q(8:14) - C*q(1:7));        % calculate qdot 8 to 14
dq=[q(8:14);dq8_14];             % return qdot vector 14x1