function [R_r,A_r,V_rated,CT_0,CP_opt,TSR_opt,M_nacellerotor,M_nacelle,M_rotor,z_hub,xCG_nacelle,zCG_nacelle,xCG_rotor] = data_IEA_Turbine()
% IEA Wind Turbine data
R_r = 120;                              % 15 MW rotor radius [m]
A_r = pi*(240/2)^2;                     % 15 MW rotor area [m^2]
V_rated = 10.59;                        % 15 MW rated wind speed [m/s]
CT_0 = 0.799;                           % 15 MW nominal thrust coefficient
CP_opt = 0.489;                         % 15 MW design "optimal" CP
TSR_opt = 9.0;                          % 15 MW design "optimal" tip speed ratio
M_nacellerotor = 1.017*10^6;            % nacelle-rotor assembly mass [kg]**** discrepancy between IEA and COREWIND
M_nacelle = 820888;                     % nacelle mass [kg]
M_rotor = M_nacellerotor - M_nacelle;   % rotor mass [kg]
z_hub = 135;                            % hub height [m]
xCG_nacelle = -5.486;                    % nacelle x center of gravity [m] (adjusted from 5.486 to match CG)
zCG_nacelle = (3.978-5.462)+135;        % nacelle z center of gravity [m]
xCG_rotor = -11.42;                     % hub xCG 10.604