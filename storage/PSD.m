function [] = PSD(t,signal,fHighCut,ylabelstr)
%% Description
% This function takes a timeseries (t) and the signal response for that
% time series as an input, and return a plot of the time  series and
% frequency domain. It also takes the cutoff frequency for plotting.
% The dimension of the signal tells the function how many subplots to make.
% The user also designates the y-axis labels for the signal being plotted
% with a matrix of strings, size 2 x (number of signals).
    %% Important PSD information for plotting only steady-state response
    % when plotting time decay, set timestartpos to 1. when plotting
    % forced response, set timestartpos to 10001
    timestartpos = 1;       % steady state time position 1
%% Implementation
% get the number of signals to be plotted, and the relevant direction of
% the signal input matrix
global w1 w5
[numbersignals,mindim] = min(size(signal));
% if needed, transpose the signal matrix
if mindim==1
    signal=signal';
end
% create subplots
for numbersubplots=1:numbersignals
    % Plot timeseries
    subplot(numbersignals,2,2*numbersubplots-1), plot(t,signal(:,numbersubplots),'LineWidth',1.25), grid on
    hold on
    % xlabel only if the last plot
    if numbersubplots==numbersignals
        xlabel('Time [s]')
    end
    ylabel(ylabelstr(1,numbersubplots))

    df = 1/(t(end)-t(timestartpos));       % Frequency resolution
    fpsd = df*(0:length(t)-timestartpos);  % Frequency vector starts from 0 for length t

    signalhat = fft(signal(timestartpos:end,numbersubplots))/length(t(timestartpos:end));          % Fourier amplitudes
    signalhat(1) = 0;                           % Discard first value (mean)
    signalhat(round(length(fpsd)/2):end) = 0;   % Discard all above Nyquist fr.
    signalhat = 2*signalhat;                    % Make amplitude one-sided
    psd = abs(signalhat).^2/2/df;               % Calculate spectrum

    % Plot frequency domain
    subplot(numbersignals,2,2*numbersubplots), plot(fpsd,psd,'LineWidth',1.25),
    xline(w1,'--r')
    xline(w5,'--m')
    grid on
    hold on, xlim([0 fHighCut])
    legend('','f_{s}=0.0083 Hz','f_{p}=0.0326 Hz')
    if numbersubplots==numbersignals
        xlabel('Frequency [Hz]')
    end
    ylabel(ylabelstr(2,numbersubplots))
end
