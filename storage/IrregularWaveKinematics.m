function [u,a]=IrregularWaveKinematics(fvec,amp) 
% This function calculates the velocity and acceleration at various heights
% for a floating spar for irregular waves
% Inputs
% fvec = wave frequency 
% h = depth of water or spar
% g = gravity
% amp = wave amplitude
% rho = density of water
% U= Horizontal Velocity
global g z_bot z t
% Given constants
for gc = 1:1
    h = -z_bot;
end
% Calculated Constants and Initializations
for cc = 1:1
    k = zeros(size(fvec));
    x = 0;
    u = zeros(length(z),length(t));
    a = zeros(length(z),length(t));
end
% Frequency domain inputs for wave number
for ifreq = 1:length(fvec)
    k(ifreq) = wave_number(fvec(ifreq),g,h);
end
% random error
random = 2*pi*rand(1,length(fvec));
% calculating acceleration and velocity
for iz=1:length(z)
    for it=1:length(t)
        uj = 0;
        aj = 0;
        for ifreq = 2:length(fvec)
            omega = 2*pi*fvec(ifreq);
            uj = uj + amp(ifreq) *omega* cosh(k(ifreq)*(z(iz)+h)) / sinh(k(ifreq)*h)*cos(omega*t(it)-(k(ifreq)*x)+(random(ifreq)));
            aj = aj - omega^2*amp(ifreq) * cosh(k(ifreq)*(z(iz)+h)) / sinh(k(ifreq)*h) * sin(omega*t(it)-(k(ifreq)*x)+(random(ifreq))); 
        end
        u(iz,it) = uj;
        a(iz,it) = aj;
    end
end