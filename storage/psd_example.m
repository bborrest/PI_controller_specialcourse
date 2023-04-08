% PSD function class example
clear all
clc
%% Run the PSD
t = [0:0.1:2000];
Hs = 8.238; %[m]
Tp = 12.8; % [s]
fHighCut = 0.2;
df = 0.1*fHighCut/(t(end)-t(1)); 
gammaJS = 3.3;
[fvec,a,S_JS] = jonswap(Hs,Tp,df,fHighCut,gammaJS);

% S_JS shouldn't be here but eta should
PSD(t,S_JS,fvec)