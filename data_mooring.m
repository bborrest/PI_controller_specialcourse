function [Xmoor,hmoor,amoor,wmoor,Zm,phimoor,moorxy,TH,TV] = data_mooring(rho_H20,h)
%% Mooring data
ld = 43.9;                                      % delta lines projected length (43.9) (47.1)
hd = 0;                                        % delta lines projected height of connection adjusted (2.48)
hmoor = 110 - hd;                               % 200m - 90m + dh
Xmoor = [600 - 9.3;600 - 9.3;600 - 9.3];        % [m]
Zm = (hmoor + hd) - h;                                 % moment arm for the mooring [m]
l = 565 + ld;                                   % line length
wmoor = 561.25*(1 - (rho_H20*pi*0.08^2)/561.25)*9.81;   % N/m
% For now, Xmoor(1) bc we are assuming 0-point equilibrium
[TH,TV,~] = inelastic_catenary_line(hmoor,Xmoor(1),l,wmoor);
TH = [TH;TH;TH];
TV = [TV;TV;TV];

amoor = TH/wmoor;

phimoor = [0; 2*pi/3; 4*pi/3];  % according to x facing downwind
moorxy = [600*cos(phimoor),600*sin(phimoor)];