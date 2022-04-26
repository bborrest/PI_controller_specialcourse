function [fvec,a,S_JS] = jonswap(Hs,Tp,df,fHighCut,gammaJS)
    % This function calculates the JONSWAP distribution for waves,
    % frequency and amplitude
    % Inputs: Hs, Period, frequency step, max frequency considered, gamma
    % Outputs: time-varying frequency, time-varying wave ampltiude, Jonswap
    % frequency spectra
    fvec = [0 : df : fHighCut];
    fp= 1/Tp;
    for i =1: length(fvec)
        if fvec(i) <= fp
            sigma = 0.07;
        else
            sigma = 0.09;
        end
        gammaexp = exp(-0.5*(((fvec(i)/fp)-1)/sigma)^2);
        S_JS(i) = 0.3125* Hs^2 *Tp * (fvec(i)/fp)^(-5)* exp(-1.25*(fvec(i)/fp)^(-4))*(1-0.287*log(gammaJS))*gammaJS^gammaexp;
        a(i) = sqrt(2*S_JS(i)*df);
    end
return