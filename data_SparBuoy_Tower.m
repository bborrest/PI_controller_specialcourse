function [z_buoy,CD,Cm_cyl,Cm_sph,IxxCM_spartower,IyyCM_spartower,IzzCM_spartower,M_spartower,z_bot,zCM_spartower,D_spar,h_moor] = data_SparBuoy_Tower
% Windcrete sparbuoy-tower data
z_buoy = -77.29;                % center of buoyancy force [m] (-113.08)
CD = 0.6;                       %*** spar-buoy coefficient of drag
Cm_cyl = 1.0;                   % added mass coefficient for morrison equation for cylinder
Cm_sph = 0.7357;                % added mass coefficient for heave
IxxCM_spartower = 1.5536*10^11; % spar-tower moment of inertia around mass center in the x-axis [kg*m^2]
IyyCM_spartower = 1.5536*10^11; % spar-tower moment of inertia around mass center in the y-axis [kg*m^2]
IzzCM_spartower = 1.9025*10^9;  % spar-tower moment of inertia around mass center in the z-axis [kg*m^2]
%M_spar = 3.655*10^7;           % total mass of the spar [kg]
M_spartower = 3.9805*10^7;      % total mass of the spar-tower [kg]
z_bot = -155;                   % bottom depth of the spar [m]
%zCM_spar = -98.41;             % spar center of mass location [m]
zCM_spartower = -98.41;         % spar-tower center of mass location [m]
D_spar = 13.60;                 % spar diamter at mean sea level [m]
h_moor = 110;                   % mooring fairlead height about seabed for 200m depth [m]