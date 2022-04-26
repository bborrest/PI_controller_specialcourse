function [K,KPg,KIg] = controller_initialization()
%% Description
% All of the initial values for the turbine's control scheme are made here
%% Region 1: Optimal Torque Control
global rho_air A_r
% rho_air           % air density
% A_r               % Rotor swept area
eta = 0.965;        % generator efficiency
R = A_r^0.5/pi;     % Rotor radius
CPopt =0.489;       % optimal CP value
% Ng        % gearbox ratio
9.0;                % optimal tip speed ratio
K = (eta*rho_air*A_r*R^3*CPopt)/(2*Ng*lambdaOPT^3);
%% Region 2: Torque Control/Speed Regulation
w_omegag = [0.05,0.03]; % natural frequency ~0.1 Hz (0.05 for torque, 0.03 for
% pitch)
zeta_omegag = 0.7;  % damping ratio [0.6,0.7] (0.7 for torque and pitch)
Ir = 1.66639*10^7;  % rotor rotational moment of inertia [kgm^2]*** needs to be updated/found with mass dist.
% Ig            % generator rotational moment of inerta
KPg = 2*eta*(zeta_omegag)*(w_omegag)*(Ir + Ng^2*Ig);
KIg = eta*(Ir + Ng^2*Ig)*(w_omegag)^2;
%% Region 3: Blade Pitch Control
% The region 3 controller values are based on KPg and KIg but are updated
% for each time step. This process is handled in dqdtsparbuoy.m