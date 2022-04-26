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
global g z_bot z_spar t
% Given constants
for gc = 1:1
    h = -z_bot;
end
% Calculated Constants and Initializations
for cc = 1:1
    k = zeros(size(fvec));
    x = 0;
    u = zeros(length(z_spar),length(t));
    a = zeros(length(z_spar),length(t));
end
% Frequency domain inputs for wave number
for ifreq = 1:length(fvec)
    k(ifreq) = wave_number(fvec(ifreq),g,h);
end
% random error
random = 2*pi*rand(1,length(fvec));
% some frequency parameter
omega2 = 2*pi*fvec(2:end);
% calculating acceleration and velocity
for iz=1:length(z_spar)
    for it=1:length(t)
%         uj = 0;
%         aj = 0;
%         for ifreq = 2:length(fvec)
%             omega = 2*pi*fvec(ifreq);
%             uj = uj + amp(ifreq) *omega* cosh(k(ifreq)*(z_spar(iz)+h)) / sinh(k(ifreq)*h)*cos(omega*t(it)-(k(ifreq)*x)+(random(ifreq)));
%             aj = aj - omega^2*amp(ifreq) * cosh(k(ifreq)*(z_spar(iz)+h)) / sinh(k(ifreq)*h) * sin(omega*t(it)-(k(ifreq)*x)+(random(ifreq))); 
%         end
%         omega2 = 2*pi*fvec(2:end);
%         u(iz,it) = uj;
%         a(iz,it) = aj;
        u(iz,it) = sum(amp(2:end).*omega2.*(cosh(k(2:end)*(z_spar(iz) + h))./sinh(k(2:end)*h)).*cos(omega2*t(it) - k(2:end)*x + random(2:end)));
        a(iz,it) = sum(-omega2.^2.*amp(2:end).*(cosh(k(2:end)*(z_spar(iz) + h))./sinh(k(2:end)*h)).*sin(omega2*t(it) - k(2:end)*x + random(2:end)));

    end
end