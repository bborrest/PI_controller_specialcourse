function [S_W,V,f] = Kaimal_Timeseries(I,L,fHigh)
%% Description
% This function calculates wind power density using the Kaimal spectrum. It
% takes wind parameter inputs:
% I: turbulence intensity
% V_10: 10-minute average wind speed [m/s]
% L: turbulence length scale [m]
% fHigh: cut-off frequency
% t: time space [s]
% and returns:
% the spectral density function, S_W
% the velocity time series, V, and 
% the freqeuency vector, f.
%% Implementation
% initializations
global t V_10
df = 1/t(end);
f = [df:df:fHigh];
rng(1)
ep = 2*pi*rand(1,length(f));
bp = zeros(1,length(f));
wp = 2*pi*f;
S_W = zeros(size(t));
V = zeros(size(t));
% calculate spectrum
for p = 1: length(f)
    S_W(p) = 4*I^2*V_10*L/(1+6*(f(p)*L/V_10))^(5/3); 
    bp(p) = sqrt(2*S_W(p)*df); 
end
% velociy time series
for i = 1:length(t)
    V(i) = V_10 + sum(bp.*cos(wp*t(i)+ep));
end
