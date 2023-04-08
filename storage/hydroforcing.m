function Fvec = hydroforcing(t_index,q3,q4)
%% Description
% This function can evaluate the forces and moments from hydrodynamic
% forcing for still water or for a given sea state. This function is
% designed for use with dqdtsparbuoy and ode4.
%% Inputs:
% t is the current time state and will be used to pull the pre-determined
% sea state
% q3 is dx0/dt
% q4 is d_theta/dt
%% Outputs:
% Hydrodynamic force and moment for a point in time
%% Implementation
global rhow Cm CD D_spar u udot z
A = D_spar^2*pi/4;
F = 0;
Tau = 0;
for i=1:length(z)
    df = rhow*((Cm+1)*A*udot(i,t_index) + 0.5*CD*D_spar*(u(i,t_index)-q3-z(i)*q4)*abs(u(i,t_index)-q3-z(i)*q4));
    dtau = df*z(i);
    F = F + df;
    Tau = Tau + dtau;
end
Fvec=[F;Tau];