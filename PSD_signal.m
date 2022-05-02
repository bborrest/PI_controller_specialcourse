function [] = PSD_signal(t,waveFSE,windVel,signals,fHighCut,wref,analysisType,analysisCondtions,ylabelstr)
%% Description
% This function takes a timeseries (t) and the signal response for that
% time series as an input, and return a plot of the time  series and
% frequency domain. It also takes the cutoff frequency for plotting. The
% function has been adapted specifically for use with outputs from ode4-
% meaning that the inputs (wind, waves) are double the size of the outputs
% to accomodate for the half-step in the runge-kutta method in ode4. For
% this reason, the wind and waves have to be handled separately.
% The user also designates the y-axis labels for the signal being plotted
% with a matrix of strings, size N x 1.
%% Important PSD information for plotting only steady-state response
% when plotting time decay, set timestartpos to 1. when plotting
% forced response, set timestartpos to 10001
timestartpos = 1;       % steady state time position 1
%% Inputs
% t:            timeseries from ode4 (normal length)
% waveFSE:      OPTIONAL input (put a scalar if not using), if the waves 
%               are included. The first row/column is the timeseries, the 
%               second is the free surface elevation.
% windVel:      OPTIONAL input (put a scalar if not using), if the wind is
%               included. The first row/column is the timeseries, the
%               second is the wind velocity.
% signals:      From output from ode4/dqdtsparbuoy, size (number of
%               signals) by (timeseries, normal length)
% fHighCut:     High cut frequency; x limit for PSD plots
% wrefall:      The reference natural frequency velocities (set to 1 if
%               ignoring)
% analysisType:         string input, describing the type of analysis
% analysisConditions:   string input, describing the analysis conditions
% ylabelstr:            string input, list of y-axis plot labels this needs
%                       to be the same length as the number of signals plus
%                       the wind and waves frequency input if included
%% Implementation
% get the number of signals to be plotted, and the relevant direction of
% the signal input matrix
[numbersignals,mindim] = min(size(signals));
% if needed, transpose the signal matrix
if mindim==1
    signals=signals';
end
% check if waves and or wind time signals are included; size the subplots
start = 1;
if length(waveFSE) > 1
    % if needed, transpose the signal matrix
    [~,mindim] = min(size(waveFSE));
    if mindim==2
        waveFSE=waveFSE';
    end
    numbersignals = numbersignals + 1;
    start = start + 1;
end
if length(windVel) > 1
    % if needed, transpose the signal matrix
    [~,mindim] = min(size(windVel));
    if mindim==2
        windVel=windVel';
    end
    numbersignals = numbersignals + 1;
    start = start + 1;
end
%% plotting
% plotting the wind and wave conditions
if (start - 1) > 0
    for numbersubplots = 1:(start-1)
        % if the waves free surface elevation is included, make first subplot
        % (needs to be done separately because it has a different time series)
        if length(waveFSE) > 1
            subplot(numbersignals,2,2*numbersubplots-1), plot(waveFSE(1,:),waveFSE(2,:),'LineWidth',1.25), grid on
            yMax = max(waveFSE(2,:));
            yMin = min(waveFSE(2,:));
            % catch for still water
            if yMax == yMin
               yMax = yMax + 1;
               yMin = yMin - 1;
            end
            ylim([yMin - 0.2*abs(yMax - yMin),yMax + 0.2*abs(yMax - yMin)]);
            hold on
            ylabel(ylabelstr(numbersubplots))
            % Generate frequency domain
            df = 1/(waveFSE(1,end)-waveFSE(1,timestartpos));       % Frequency resolution
            fpsd = df*(0:length(waveFSE)-timestartpos);  % Frequency vector starts from 0 for length t
            % +1 in the index here to get the free surface elevation signal
            % instead of time
            signalhat = fft(waveFSE(2,timestartpos:end))/length(waveFSE(1,timestartpos:end));          % Fourier amplitudes
            signalhat(1) = 0;                           % Discard first value (mean)
            signalhat(round(length(fpsd)/2):end) = 0;   % Discard all above Nyquist fr.
            signalhat = 2*signalhat;                    % Make amplitude one-sided
            psd = abs(signalhat).^2/2/df;               % Calculate spectrum
            % Determine maximum x-axis for plotting
            xMaxLim = fHighCut;
            % Plot frequency domain
            subplot(numbersignals,2,2*numbersubplots), plot(fpsd,psd,'LineWidth',1.25),
            grid on
            hold on, xlim([0 xMaxLim])
            %ylabel(ylabelstr(numbersubplots))
            ylabel('PSD')
        end
        % if the wind velocity is included, make the first/second subplot
        % create subplots (why numbersubplots == start is included)
        if length(windVel) > 1 && numbersubplots == (start-1)
            subplot(numbersignals,2,2*numbersubplots-1), plot(windVel(1,:),windVel(2,:),'LineWidth',1.25), grid on
            yMax = max(windVel(2,:));
            yMin = min(windVel(2,:));
            % catch for constant wind
            if yMax == yMin
               yMax = yMax*1.2;
               yMin = yMin*0.8;
            end
            ylim([yMin - 0.2*abs(yMax - yMin),yMax + 0.2*abs(yMax - yMin)]);
            hold on
            ylabel(ylabelstr(numbersubplots))
            % Generate frequency domain
            df = 1/(windVel(1,end)-windVel(1,timestartpos));       % Frequency resolution
            fpsd = df*(0:length(windVel)-timestartpos);  % Frequency vector starts from 0 for length t
            signalhat = fft(windVel(2,timestartpos:end))/length(windVel(timestartpos:end,1));          % Fourier amplitudes
            signalhat(1) = 0;                           % Discard first value (mean)
            signalhat(round(length(fpsd)/2):end) = 0;   % Discard all above Nyquist fr.
            signalhat = 2*signalhat;                    % Make amplitude one-sided
            psd = abs(signalhat).^2/2/df;               % Calculate spectrum
            % Determine maximum x-axis for plotting
            xMaxLim = fHighCut;
            % Plot frequency domain
            subplot(numbersignals,2,2*numbersubplots), plot(fpsd,psd,'LineWidth',1.25),
            grid on
            hold on, xlim([0 xMaxLim])
            %ylabel(ylabelstr(numbersubplots))
            ylabel('PSD')
        end
    end
end
% plotting the sparbuoy response
for numbersubplots = start:numbersignals
    signalsIndex = size(signals,2) + numbersubplots - numbersignals;
    subplot(numbersignals,2,2*numbersubplots-1), plot(t,signals(:,signalsIndex),'LineWidth',1.25), grid on
    yMax = max(signals(:,signalsIndex));
    yMin = min(signals(:,signalsIndex));
    ylim([yMin - 0.2*abs(yMax - yMin),yMax + 0.2*abs(yMax - yMin)]);
    hold on
    % xlabel only if the last plot
    if numbersubplots==numbersignals
        xlabel('Time [s]')
    end
    ylabel(ylabelstr(numbersubplots))
    
    % Generate frequency domain
    df = 1/(t(end)-t(timestartpos));       % Frequency resolution
    fpsd = df*(0:length(t)-timestartpos);  % Frequency vector starts from 0 for length t

    signalhat = fft(signals(timestartpos:end,signalsIndex))/length(t(timestartpos:end));          % Fourier amplitudes
    signalhat(1) = 0;                           % Discard first value (mean)
    signalhat(round(length(fpsd)/2):end) = 0;   % Discard all above Nyquist fr.
    signalhat = 2*signalhat;                    % Make amplitude one-sided
    psd = abs(signalhat).^2/2/df;               % Calculate spectrum

    % Determine maximum x-axis for plotting
    xMaxLim = fHighCut;

    % Plot frequency domain
    wrefIndex = numbersubplots + length(wref) - numbersignals;
    subplot(numbersignals,2,2*numbersubplots), plot(fpsd,psd,'LineWidth',1.25),
    xline(wref(wrefIndex),'--r')
    grid on
    hold on, xlim([0 xMaxLim])
    legend('',['f_{s}=',num2str(wref(wrefIndex))],'Location','best')
    if numbersubplots==numbersignals
        xlabel('Frequency [Hz]')
    end
    %ylabel(ylabelstr(numbersubplots))
    ylabel('PSD')
end
plotTitle = strcat(analysisType,", ",analysisCondtions);
sgtitle(plotTitle)