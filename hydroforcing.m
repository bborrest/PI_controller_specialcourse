function Fvec = hydroforcing(t_index,q8,q9,q11,q12)
%% Description
% This function can evaluate the forces and moments from hydrodynamic
% forcing for still water or for a given sea state. This function is
% designed for use with dqdtsparbuoy and ode4.
% NOTE: The diameter used for the spar is unique to the WINDCRETE design
%% Inputs:
% t is the current time state and will be used to pull the pre-determined
% sea state
% q8 is d_x1/dt (surge)
% q9 is d_x2/dt (sway)
% q11 is d_x4/dt (roll)
% q12 is d_x5/dt (pitch)
%% Outputs:
% Hydrodynamic force and moment for a point in time, as a 6x1 vector
%% Implementation
global rho_H2O Cm_cyl CD D_spar u udot v vdot z_spar
% intialize F and Tau (assume no z forcing, yaw forcing)
F_surge = 0;
F_sway = 0;
Tau_roll = 0;
Tau_pitch = 0;
for i=1:length(z_spar)
    % define appropriate diameter
    if z_spar(i)<(-145.7)
        D = 2*sqrt(9.3^2 - (145.7 + z_spar(i))^2);   % spherical cap
    elseif z_spar(i)>=(-145.7) && z_spar(i)<=(-10)
        D = 18.60;  % cylindrical section
    elseif z_spar(i)>(-10)
        D = D_spar - 0.54*z_spar(i); % conical section; D_spar at mean sea level
    end
    % calculate area
    A = D^2*pi/4;
    % forces and stuff
    dfsurge = rho_H2O*((Cm_cyl+1)*A*udot(i,t_index) + 0.5*CD*D_spar*(u(i,t_index)-q8-z_spar(i)*q12)*abs(u(i,t_index)-q8-z_spar(i)*q12));
    dfsway = rho_H2O*((Cm_cyl+1)*A*vdot(i,t_index) + 0.5*CD*D_spar*(v(i,t_index)-q9-z_spar(i)*q11)*abs(v(i,t_index)-q9+z_spar(i)*q11));
    dtau_pitch = dfsurge*z_spar(i);
    dtau_roll = -dfsway*z_spar(i);
    F_surge = F_surge + dfsurge;
    F_sway = F_sway + dfsway;
    Tau_pitch = Tau_pitch + dtau_pitch;
    Tau_roll = Tau_roll + dtau_roll;
end
Fvec=[F_surge; F_sway; 0; Tau_roll; Tau_pitch; 0; 0];