function [u,a,eta]=LinearWaveKinematics() 
% This function calculates the kinematics of regular waves
% Inputs
% f= wave frequency 
% h= depth of water or spar
% g= gravity
% rho = density of water
% U= Horizontal Velocity
global g z_bot Hs Tp z_spar t
% Given constants
    h = -z_bot;
    f = 1/Tp;
    H = Hs;
% Calculated Constants
    w = 2*pi*f;
    k = wave_number(f,g,h);
    % Pre calculations
    x3 = 0;
    u = zeros(length(z_spar),length(t));
    a = zeros(length(z_spar),length(t));
    for j = 1:length(z_spar)
        for i = 1:length(t)          
            u(j,i) = w*H/2 * cosh(k*(z_spar(j)+h)) / sinh(k*h) * cos(w*t(i)-k*x3);
            a(j,i) = -w^2*H/2 * cosh(k*(z_spar(j)+h)) / sinh(k*h) * sin(w*t(i)-k*x3);    
         end
    end
    eta = H/2 * cos(w*t-k*x3);