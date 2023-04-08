% Assignment 5Josh & Varun
clc; clear all
% Loading the given constants for the floater
global rhow Cm CD D_spar z_bot g z_hub
load('model_constants.mat');
%% General Constants
global t Hs Tp rho_air A_r V_rated CT_0 aCT bCT gammaCT
for ig=1:1
    fHighCut = 0.5;         % cut-off frequency
    Tdur = 1600;            % total duration: 1000s transient + 600s response
    dtode = 0.1;            % time-step for ode4 solver
    dt = 0.05;              % general time step
    tode4 = [0:dtode:Tdur-dtode];
    t = [0:dt:Tdur-dt];
    z_buoy = z_bot/2;       % center of buoyancy force
    Hs = 6;                 % linear wave amplitude and significant wave height [m]
    Tp = 10;                % linear wave period and significant wave period [s]
    rho_air = 1.22;         % air density [kg/m^3]
    A_r = 24885;            % 10 MW rotor area [m^2]
    V_rated = 11.4;         % 10 MW rated wind speed [m/s]
    CT_0 = 0.81;            % 10 MW nominal thrust coefficient
    aCT = 0.5;              % 10 MW "a" thrust parameter
    bCT = 0.65;             % 10 MW "b" thrust paremeter
    gammaCT = 2;            % initial controller value
    gammaJS = 3.3;          % JONSWAP peak enhancement factor
    df = 1/Tdur;            % JONSWAP frequency spectra time step
    TI = 0.14;              % Kaimal wind spectra turbulence intensity
    TL = 340.2;             % Kaimal wind spectra turbulenec length scale [m]
end
%% Part 1, Model Formulation
% Question 1 : Mtot, zCMtot and IOtot calculation
for i1=1:1
    M_tot  = M_floater+ M_tower+M_nacelle+M_rotor;
    zCM_tot=(M_floater*zCM_floater + M_tower*zCM_tower + (M_nacelle+M_rotor)*z_hub)/ M_tot;
    IO_tot = (ICM_floater+ M_floater*zCM_floater^2)+(ICM_tower+M_tower*zCM_tower^2)+(M_nacelle+M_rotor)*z_hub^2; 
end
% Question 2,3,4,5 : written in the report
% Question 6 : Building ODE system
for i6=1:1
    I11A = pi/64*D_spar^4;                                                  % second moment of area at the waterplane
    M = [M_tot, M_tot*zCM_tot; 
        M_tot*zCM_tot, IO_tot];                      % mass matrix
    A = [-pi/4*rhow*D_spar^2*Cm*z_bot, -pi/8*rhow*D_spar^2*Cm*z_bot^2; 
        -pi/8*rhow*D_spar^2*Cm*z_bot^2, -pi/12*rhow*D_spar^2*Cm*z_bot^3];   % added mass matrix
    C = [K_moor, K_moor*z_moor;
        K_moor*z_moor, rhow*g*I11A+M_tot*g*(z_buoy-zCM_tot)+K_moor*z_moor^2];
end
% Question 7: Estimating natural frequencies of surge and pitch
global w1 w5
for i7=1:1
    [~,w15]=eig((M+A)^-1*C);
    w1=w15(1,1)^0.5/(2*pi);             % Hz
    w5=w15(2,2)^0.5/(2*pi);             % Hz
end
%% Part 2, Dynamic Analysis
% Question 8: written in the report
% Question 9: preparing for ode4 solver
for i9=1:1
    B = [2*10^5, 0; 0, 0];                              % Damping matrix
end
% Question 10: implementing ode4/dqdt algorithm
for i10=1:1
    q01 = [1;0;0;0;0];                                  % first initial condition
    q02 = [0;0.1;0;0;0];                                % second initial condition
    Y10_1 = ode4(@dqdtsparbuoy,tode4,q01,M+A,B,C,1);    % unforced response 1
    Y10_2 = ode4(@dqdtsparbuoy,tode4,q02,M+A,B,C,1);    % unforced response 2
    Y10_1(:,2)=Y10_1(:,2)*180/pi;                       % convert pitch to degrees
    Y10_2(:,2)=Y10_2(:,2)*180/pi;
    % plotting PSD
    labels10=["Surge, x_0 [m]","Pitch, \theta [deg]";"PSD [m^2/Hz]","PSD [deg^2/Hz]"];
    figure
    PSD(tode4,Y10_1(:,1:2),fHighCut,labels10)
    figure
    PSD(tode4,Y10_2(:,1:2),fHighCut,labels10)
end
% Question 11: with hyrdodynamic forcing, assuming no current u = 0
global u udot z
for i11=1:1
    % part a
    z = [0:-1:z_bot];
    u=zeros(length(z),length(t));
    udot=u;
    Y11_1 = ode4(@dqdtsparbuoy,tode4,q01,M+A,B,C,2);%hydro forcing response 1
    Y11_2 = ode4(@dqdtsparbuoy,tode4,q02,M+A,B,C,2);%hydro forcing response 2
    Y11_1(:,2)=Y11_1(:,2)*180/pi;                   % convert pitch to degrees
    Y11_2(:,2)=Y11_2(:,2)*180/pi;
    % plotting PSD
    figure
    PSD(tode4,Y11_1(:,1:2),fHighCut,labels10)
    figure
    PSD(tode4,Y11_2(:,1:2),fHighCut,labels10)
    % part b
    q03 = [0;1;0;0;0];
    Y11_3 = ode4(@dqdtsparbuoy,tode4,q03,M+A,B,C,1);%no forcing response
    Y11_4 = ode4(@dqdtsparbuoy,tode4,q03,M+A,B,C,2);%hydro forcing response
    Y11_3(:,2)=Y11_3(:,2)*180/pi;                   % convert pitch to degrees
    Y11_4(:,2)=Y11_4(:,2)*180/pi;                   % convert pitch to degrees
    % plotting PSD
    figure
    PSD(tode4,Y11_3(:,1:2),fHighCut,labels10)
    PSD(tode4,Y11_4(:,1:2),fHighCut,labels10)
end
% Question 12: Linear Wave Sea-State Hydrodynamic forcing
for i12=1:1
    [u,udot]=LinearWaveKinematics();
    q00=[0;0;0;0;0];                                  % zero initial condition
    Y12 = ode4(@dqdtsparbuoy,tode4,q00,M+A,B,C,2);  % hydro forcing response 0
    Y12(:,2)=Y12(:,2)*180/pi;                       % convert pitch to degrees
    % plotting PSD
    figure
    PSD(tode4(10001:16000),Y12(10001:16000,1:2),fHighCut,labels10)       % make sure to eliminate transient
%     ylim([-1,1])
end
% Question 13: V10 Wind and Linear Wave Forcing
global V_10 V_hub
for i13=1:1
    V_10 = 8;
    V_hub = 8*ones(size(t));
    Y13 = ode4(@dqdtsparbuoy,tode4,q00,M+A,B,C,3);  % hydro forcing response 0
    Y13(:,2)=Y13(:,2)*180/pi;                       % convert pitch to degrees
    % plotting PSD
    figure
    PSD(tode4(10001:16000),Y13(10001:16000,1:2),fHighCut,labels10)       % make sure to eliminate transient
end
% Question 14: V10 Wind and Jonswap Wave Forcing
for i14=1:1
    [fvec,amp,S_JS] = jonswap(Hs,Tp,df,fHighCut,gammaJS);
    [u,udot]=IrregularWaveKinematics(fvec,amp) ;
    Y14 = ode4(@dqdtsparbuoy,tode4,q00,M+A,B,C,3);  % hydro forcing response 0
    Y14(:,2) = Y14(:,2)*180/pi;                     % convert pitch to degrees
    % plotting PSD
    figure
    PSD(tode4(10001:16000),Y14(10001:16000,1:2),fHighCut,labels10)       % make sure to eliminate transient
end
% Question 15: Kaimal Wind and Jonswap Wave Forcing
for i15=1:1
    [~,V_hub,~] = Kaimal_Timeseries(TI,TL,fHighCut);
    Y15 = ode4(@dqdtsparbuoy,tode4,q00,M+A,B,C,3);  % hydro forcing response 0
    Y15(:,2) = Y15(:,2)*180/pi;                     % convert pitch to degrees
    % plotting PSD
    figure
    PSD(tode4(10001:16000),Y15(10001:16000,1:2),fHighCut,labels10)       % make sure to eliminate transient
end
%% Part 3: Adaptation of pitch control for dynamic stability
% Adapting the CT pitch controller for a floating configuration
% Question 16: Steady wind and no waves
for i16=1:1
    % 10 m/s case no waves
    V_10 = 10;
    V_hub = 10*ones(size(t));
    u=zeros(length(z),length(t));
    udot=u;
    Y16_1 = ode4(@dqdtsparbuoy,tode4,q00,M+A,B,C,3);  % hydro forcing response 0
    Y16_1(:,2)=Y16_1(:,2)*180/pi;                       % convert pitch to degrees
    % 16 m/s case no waves
    V_10 = 16;
    V_hub = 16*ones(size(t));
    Y16_2 = ode4(@dqdtsparbuoy,tode4,q00,M+A,B,C,3);  % hydro forcing response 0
    Y16_2(:,2)=Y16_2(:,2)*180/pi;                       % convert pitch to degrees
    % plotting PSD
    figure
    PSD(tode4(10001:16000),Y16_1(10001:16000,1:2),fHighCut,labels10)       % make sure to eliminate transient
    figure
    PSD(tode4(10001:16000),Y16_2(10001:16000,1:2),fHighCut,labels10)       % make sure to eliminate transient
end
% Question 17: CT adapted controller model
for i17=1:1
    % 16 m/s case no waves
    gammaCT = 0.1;
    q00contr = [0;0;0;0;0];
    Y17 = ode4(@dqdtsparbuoy,tode4,q00contr,M+A,B,C,4);  % hydro forcing response 0
    Y17(:,2) = Y17(:,2)*180/pi;                       % convert pitch to degrees
    % plotting PSD
    figure
    PSD(tode4(10001:16000),Y17(10001:16000,1:2),fHighCut,labels10)       % make sure to eliminate transient
end
% Question 18: Proper gamma for CT controller model
for i18=1:1
    gammaCTrange = [0.1,0.5,1.5];
    Y18 = zeros(length(Y17),length(gammaCTrange));
    figure
    for igam = 1:length(gammaCTrange)
        gammaCT = gammaCTrange(igam);
        Y18i = ode4(@dqdtsparbuoy,tode4,q00contr,M+A,B,C,4);    % hydro forcing response 0
%         Y18(:,igam) = Y18i(:,2);                                % looking only at surge
        Y18i(:,2) = Y18i(:,2)*180/pi;                           % convert pitch to degrees
        PSD(tode4,Y18i(:,1:2),fHighCut,labels10)
        hold on
    end
end