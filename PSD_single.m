function [] = PSD_single(t,signal,fHighCut,wref,analysisType,ylabelstr)
%% Description
% This function takes a timeseries (t) and the signal response for that
% time series as an input, and return a plot of the time  series and
% frequency domain. It also takes the cutoff frequency for plotting.
% This function has been adapted to plot one signal for each call.
% The user also designates the y-axis labels for the signal being plotted
% with a matrix of strings, size 2 x 1.
%% Important PSD information for plotting only steady-state response
% when plotting time decay, set timestartpos to 1. when plotting
% forced response, set timestartpos to 10001
timestartpos = 1;       % steady state time position 1
%% Implementation
global V_10
% Plot timeseries
subplot(2,1,1), plot(t,signal,'LineWidth',1.25), grid on
hold on
xlabel('Time [s]')
ylabel(ylabelstr(1))
title(analysisType,[num2str(V_10),'m/s wind'])
% Generate frequency domain
df = 1/(t(end)-t(timestartpos));       % Frequency resolution
fpsd = df*(0:length(t)-timestartpos);  % Frequency vector starts from 0 for length t

signalhat = fft(signal(timestartpos:end))/length(t(timestartpos:end));          % Fourier amplitudes
signalhat(1) = 0;                           % Discard first value (mean)
signalhat(round(length(fpsd)/2):end) = 0;   % Discard all above Nyquist fr.
signalhat = 2*signalhat;                    % Make amplitude one-sided
psd = abs(signalhat).^2/2/df;               % Calculate spectrum

% Determine maximum x-axis for plotting
xMaxLim = 0.4*fHighCut;
% if 1.5*max(wref) < fHighCut
%     xMaxLim = 1.5*max(wref);
% else
%     xMaxLim = fHighCut;
% end

% Plot frequency domain
subplot(2,1,2), plot(fpsd,psd,'LineWidth',1.25),
xline(wref,'--r')
grid on
hold on, xlim([0 xMaxLim])
legend('',['f_{s}=',num2str(wref)],'Location','best')
xlabel('Frequency [Hz]')
ylabel(ylabelstr(2))
