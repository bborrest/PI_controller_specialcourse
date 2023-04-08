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
% As it is currently implemented, M,B,C are all 2x2 matrices; F is a 2x1
% forcing vector with Force in 1 and torque in 2; q is a 4x1 vector that is
% split into two 2x1 vectors for solving the system. The return, dq, is 4x1
%% Implementation
% Find the proper time for forcing
global t z_hub gammaCT
[~ ,index]= min(abs(tode4-t));
% Dummy CT_vrel if not using controller
CT_vrel = 0;
% Calculate forcing
if forcing == 1
    F = [0;0];                                                              % no external forcing
elseif forcing == 2
    F = hydroforcing(index,q(3),q(4));                                      % hydro forcing only
elseif forcing == 3
    q3hub = q(3)+z_hub*q(4);                            % dxdt of the hub
    [Windforcing,~] = F_wind_timepoint(q3hub,index);
    F = hydroforcing(index,q(3),q(4)) + Windforcing;                        % hydro plus steady wind forcing
elseif forcing == 4
    q3hub = q(3)+z_hub*q(4);                            % dxdt of the hub
    [Windforcing,CT_vrel] = F_wind_controller(q3hub,q(5),index);    % time-delay controller
    F = hydroforcing(index,q(3),q(4)) + Windforcing;                        % hydro plus steady wind forcing
end
dq34 = (M)^-1*(F-B*q(3:4)-C*q(1:2));        % calculate qdot 3 and qdot 4
dq5 = -gammaCT*(q(5)-CT_vrel);
dq=[q(3);q(4);dq34(1);dq34(2);dq5];             % return qdot vector 4x1