%% PI Controller for 6-DOF floating wind turbine model
% Author: Joshua Forrest, s212447@student.dtu.dk
% Date: 28/03/2022
% Description: This code is for tuning a PI controller to be used with the
% IEA 15 MW turbine situated on the windcrete spar-buoy floater. The model
% has 6 DOF, uses a single steady rotor thrust resolved at the hub, linear
% unsteady wave forcing, linearized mooring forces, PI controller,... 
% [update with progress]
clear all
clc
%% TO - DO
% finalize controller model
% controller tuning
% incorporate wind yaw
% clean up code/inputs
%% Control Conditions
global V_10
V_10 = 16;      % V10 min value
Kaimal = 0; 	% Kaimal Wind Spectrum (1) or Constant Wind (0)
Jonswap = 0;    % Jonswap Wave Distribution (1) or Still Water (0) or Regular Waves (2)
YawAngle = 0;   % Wind yaw angle
WaveAngle = 0;  % Wave angle relative to the rotor [degrees]
% IC for Decay [m,m,m,rad,rad,rad,controller,m/s,m/s,m/s,rad/s,rad/s,rad/s,dcontroller]
q0 = [0; 0; 0; 0.0; 0.0; 0.0; 0.0; 0; 0; 0; 0; 0; 0; 0];       % 14x1 (7 x, 7 dx, 6 rigid-body-DOF + azimuth)
% Controller Tuning
ctrl_damping = 0.7;     % 70% critical is the report value
ctrl_omega = 0.1885;    % 0.03 Hz is the report value [rad/s]
% Simulations to run
freeResponse = 0;   % 1 = on, 0 = off
Waves = 0;          % Just waves, 1 = on, 0 = off
Wind = 0;           % Just wind, 1 = on, 0 = off
WindWaves = 1;      % Wind + Waves, 1 = on, 0 = off
%% Constants/Inputs
global t dtode g
global rho_H2O Cm_cyl CD D_spar z_spar z_bot Tp Hs
global rho_air A_r V_rated V_hub CT_0 z_hub
for ig=1:1
% Time
    Tdur = 1600;                % total duration: 1000s transient + 600s response
    dtode = 0.1;                % time-step for ode4 solver
    dt = 0.05;                  % general time step
    tode4 = 0:dtode:Tdur-dtode;
    t = 0:dt:Tdur-dt;
% IEA Wind Turbine
    [R_r,A_r,V_rated,CT_0,CP_opt,TSR_opt,M_nacellerotor,M_nacelle,M_rotor,z_hub,xCG_nacelle,zCG_nacelle,xCG_rotor] = data_IEA_Turbine();
% Wind State
    [fHighCut,rho_air,TI,TL] = data_Wind_State();
% Windcrete sparbuoy-tower
    [z_buoy,CD,Cm_cyl,Cm_sph,IxxCM_spartower,IyyCM_spartower,IzzCM_spartower,M_spartower,z_bot,zCM_spartower,D_spar,hmoor] = data_SparBuoy_Tower();
    % Sea State
    [Hs,Tp,gammaJS,df,rho_H2O] = data_Sea_State(Tdur);
% Physical
    g = 9.81;                   % G [m/s]
    h = 200;                    % water depth [m]
    z_spar = 0:-1:z_bot;        % spar buoy integration over z for forcing
end
%% Controller Implementation
global kp3 ki3 CT_10 pTpTh0 pTpOm0 pTpV pQpU0
for ic = 1:1
    % Load Rotor Torque Derivatives and Data
    load dCPdTh
    load dCPdTSR
    load dCTdTh
    load dCTdTSR
    load IEA_15MW_HWIND_Land_Based_hs2.mat
    % Load and calculate operating points
    op0 = readmatrix('operation_equinor.txt');
    % create tsr-cp table for pQpU0
    op0_tsrCP = [op0(:,3)*R_r*2*pi/60./op0(:,1), op0(:,4)*10^3./(0.5*rho_air*A_r*op0(:,1).^3)];
    % determine operation values of tsr and pitch for given conditions
    tsrOP = R_r*(interp1(op0(:,1),op0(:,3),V_10)*2*pi/60)/V_10;
    pitchOP = interp1(op0(:,1),op0(:,2),V_10);
    % check for min operating pitch
    if pitchOP < 0.000534949676951124
        pitchOP = 0.000534949676951124;
    end
    % calculate CT_10, or operating CT
    CT_10 = interp1(turbine.pitchList,interp1(turbine.tsrList,turbine.Ct,tsrOP),pitchOP);
    % calculate operating point rates of change for CP, CT
    dCPdTh0 = interp1(turbine.tsrList,interp1(turbine.pitchList(1:24),dCPdTh',pitchOP),tsrOP);
    dCPdTSR0 = interp1(turbine.tsrList(1:55),interp1(turbine.pitchList,dCPdTSR',pitchOP),tsrOP);
    dCTdTh0 = interp1(turbine.tsrList,interp1(turbine.pitchList(1:24),dCTdTh',pitchOP),tsrOP);
    dCTdTSR0 = interp1(turbine.tsrList(1:55),interp1(turbine.pitchList,dCTdTSR',pitchOP),tsrOP);
    % Calculate Operating Partial Derivatives
    pQpTh0 = 0.5*rho_air*A_r*V_10^3*dCPdTh0/(tsrOP*V_10/R_r);       % variation of aerodynamic torque with pitch angle, linearized around 0 degree pitch [N*m/deg]
    pTpTh0 = 0.5*rho_air*A_r*V_10^2*dCTdTh0;                        % variation of the thrust with pitch, linearized around 0 degree pitch [N/deg]
    pQpOm0 = 0.5*rho_air*A_r*V_10^3*(dCPdTSR0*R_r/V_10)/(tsrOP*V_10/R_r);      % variation of aerodynamic torque with rotor speed, linearized around 0 degree pitch [N*m*s/rad]
    pTpOm0 = 0.5*rho_air*A_r*V_10^2*(dCTdTSR0*R_r/V_10);                       % variation of the thrust with rotor speed [N/rad]
    % Refer to notes for equation, but the parentheses include (CP/tsr^3 - CP2/tsr2^3)
    % but first find CP1 and CP2
    CP1 = interp1(op0(:,1),op0_tsrCP(:,2),V_10);
    CP2 = CP1 + dCPdTh0*((R_r*(interp1(op0(:,1),op0(:,3),V_10+0.1)*2*pi/60)/(V_10+0.1)) - tsrOP);  % TSR change assuming +0.1 m/s in wind change (CP1 + pTSRpCP * deltaTSR)
    pQpU0 = 0.5*rho_air*A_r*R_r^3*(tsrOP*V_10/R_r)^2*(CP1/(tsrOP^3) - CP2/(tsrOP*V_10/(V_10 + 0.1))^3);        % variation of generator torque with rotor speed , linearized around 0 degree pitch [N*m*s/rad]
    % Calculating change in T for change in wind speed
    pTpV = 0.5*rho_air*A_r*V_10^2*(-dCTdTSR0*tsrOP/V_10);
    % Constants
    eta_gen = 0.965;    % generator efficiency
    Ng = 1;             % generator gearbox ratio (1 for direct drive)
    Irg = 3.14655*10^8; % Drivetrain (Rotor + Generator) Inertia [kgm^2]
    % Region 1, below rated wind speed (Optimal Torque Control)
    K1 = (eta_gen*rho_air*A_r*R_r^3*CP_opt)/(2*Ng*TSR_opt);
    % Region 3, above rated wind speed (Blade Pitch for Pitch Stability Control)
    kp3 = ctrl_damping*ctrl_omega*(2*Irg)/(-pQpTh0);   % possibly some small error
    ki3 = ctrl_omega^2*(Irg)/(-pQpTh0);
end
%% 6DOF Model Formation
for im = 1:1
    % Mass Matrix
    [m_tot,CG,Ibxx,Ibyy,Ibzz,Ib13] = Mtot_CG_Ib_calculations(M_nacellerotor,M_nacelle,M_rotor,M_spartower,z_hub,xCG_rotor,xCG_nacelle,zCG_nacelle,zCM_spartower);
    [M] = mass_matrix(m_tot,CG,Ibxx,Ibyy,Ibzz,Ib13);
    % Added Mass
    [A] = addedmass_matrix();
    % Linearization Point
    x1 = 0; % surge [m]
    x2 = 0; % sway [m]
    x3 = 0; % heave [m]
    x4 = 0.00; % roll [rad]
    x5 = 0.00; % pitch [rad]
    x6 = 0.00; % yaw [rad]
    yawyaw = 7.10*10^8;
    %[Cmoor] = mooring_matrix(Xmoor,hmoor,amoor,wmoor,Zm,phimoor,moorxy,x1,x2,x3,x4,x5,x6,TH,TV,yawyaw);
    [Cmoor] = mooring_matrix2([x1;x2;x3;x4;x5;x6],hmoor,h,yawyaw);
    [Chydro] = hydrostatics_matrix(rho_H2O,g,D_spar,m_tot,CG,z_buoy);
    [C] = Chydro + Cmoor;
    %C(3,5) = 0;
    [B] = zeros(6);
end
%% Natural Frequencies
for iwf = 1:1
    [V,omega] = eig(C*(M+A)^-1);
    omega = diag(omega).^0.5/(2*pi);
    wref = [0.01221, 0.03052, 0.02441, 0.09155];
end
%% Add Controller to System, 7 DOF
for ic = 1:1
    M(7,7) = Irg;
    A(7,7) = 0;     % accounted for in M
    B(7,7) = -kp3*pQpTh0 - pQpOm0;
    C(7,7) = -ki3*pQpTh0;
end
%% System Unforced Response
% Currently doesn't incorporate water forcing/presence
if freeResponse == 1
    for id = 1:1
        Y_decay = ode4(@dqdtsparbuoy,tode4,q0,M+A,B,C,1);
        % time series and frequency domain plotting
        figure
        PSD_single(tode4,Y_decay(:,1),fHighCut,wref(1),"No Forcing",["Surge [m]","PSD [m^2/Hz]"]);
        figure
        PSD_single(tode4,Y_decay(:,2),fHighCut,wref(1),"No Forcing",["Sway [m]","PSD [m^2/Hz]"]);
        figure
        PSD_single(tode4,Y_decay(:,3),fHighCut,wref(2),"No Forcing",["Heave [m]","PSD [m^2/Hz]"]);
        figure
        PSD_single(tode4,Y_decay(:,4),fHighCut,wref(3),"No Forcing",["Roll [deg]","PSD [deg^2/Hz]"]);
        figure
        PSD_single(tode4,Y_decay(:,5)*180/pi,fHighCut,wref(3),"No Forcing",["Pitch [deg]","PSD [deg^2/Hz]"]);
        figure
        PSD_single(tode4,Y_decay(:,6),fHighCut,wref(4),"No Forcing",["Yaw [deg]","PSD [deg^2/Hz]"]);
    end
end
%% System Response to Wave Forcing
global u udot v vdot
if Waves == 1
    for iwaves = 1:1
        if Jonswap == 0
            % no waves, just water presence
            u = zeros(length(z_spar),length(t));
            udot = u;
            v = zeros(length(z_spar),length(t));
            vdot = v;
        elseif Jonswap == 1
            % JONSWAP wave spectrum for forcing
            [fvec,amp,~] = jonswap(Hs,Tp,df,fHighCut,gammaJS);
            [w,wdot] = IrregularWaveKinematics(fvec,amp);
            u = w*cosd(WaveAngle);
            udot = wdot*cosd(WaveAngle);
            v = w*sind(WaveAngle);
            vdot = wdot*sind(WaveAngle);
        elseif Jonswap == 2
            % linear wave forcing
            [w,wdot]=LinearWaveKinematics();
            u = w*cosd(WaveAngle);
            udot = wdot*cosd(WaveAngle);
            v = w*sind(WaveAngle);
            vdot = wdot*sind(WaveAngle);
        end
        % forced response
        Y_waves = ode4(@dqdtsparbuoy,tode4,q0,M+A,B,C,2);
        % time series and frequency domain plotting
        if Jonswap == 1 || Jonswap == 2
            figure
            PSD_single(t,u(1,:),fHighCut*1.5,1,"Wave Frequency",["MSL Velocity [m/s]","PSD [(m/s)^2/Hz]"]);
        end
        figure
        PSD_single(tode4,Y_waves(:,1),fHighCut,wref(1),"Wave Forcing only",["Surge [m]","PSD [m^2/Hz]"]);
        figure
        PSD_single(tode4,Y_waves(:,2),fHighCut,wref(1),"Wave Forcing only",["Sway [m]","PSD [m^2/Hz]"]);
        figure
        PSD_single(tode4,Y_waves(:,3),fHighCut,wref(2),"Wave Forcing only",["Heave [m]","PSD [m^2/Hz]"]);
        figure
        PSD_single(tode4,Y_waves(:,4),fHighCut,wref(3),"Wave Forcing only",["Roll [deg]","PSD [deg^2/Hz]"]);
        figure
        PSD_single(tode4,Y_waves(:,5),fHighCut,wref(3),"Wave Forcing only",["Pitch [deg]","PSD [deg^2/Hz]"]);
        figure
        PSD_single(tode4,Y_waves(:,6),fHighCut,wref(4),"Wave Forcing only",["Yaw [deg]","PSD [deg^2/Hz]"]);
    end
end
%% System Response to Wind Forcing
if Wind == 1
    for iwind = 1:1
        if Kaimal == 1
            [~,V_hub,~] = Kaimal_Timeseries(TI,TL,fHighCut);
        elseif Kaimal == 0
            V_hub = ones(size(t))*V_10;
        end
        % forced response
        Y_wind = ode4(@dqdtsparbuoy,tode4,q0,M+A,B,C,3);
        % change angles to degrees
        Y_wind(:,4:7) = Y_wind(:,4:7)*180/pi;
        % time series and frequency domain plotting
        if Kaimal == 1
            figure
            PSD_single(t(20001:end),V_hub(20001:end),fHighCut*0.5,1,"Hub Wind Velocity",["Velocity [m/s]","PSD [(m/s)^2/Hz]"]);
        end
        figure
        PSD_single(tode4(10001:end),Y_wind(10001:end,1),fHighCut,wref(1),"Wind Forcing only",["Surge [m]","PSD [m^2/Hz]"]);
        figure
        PSD_single(tode4(10001:end),Y_wind(10001:end,2),fHighCut,wref(1),"Wind Forcing only",["Sway [m]","PSD [m^2/Hz]"]);
        figure
        PSD_single(tode4(10001:end),Y_wind(10001:end,3),fHighCut,wref(2),"Wind Forcing only",["Heave [m]","PSD [m^2/Hz]"]);
        figure
        PSD_single(tode4(10001:end),Y_wind(10001:end,4),fHighCut,wref(3),"Wind Forcing only",["Roll [deg]","PSD [deg^2/Hz]"]);
        figure
        PSD_single(tode4(10001:end),Y_wind(10001:end,5),fHighCut,wref(3),"Wind Forcing only",["Pitch [deg]","PSD [deg^2/Hz]"]);
        figure
        PSD_single(tode4(10001:end),Y_wind(10001:end,6),fHighCut,wref(4),"Wind Forcing only",["Yaw [deg]","PSD [deg^2/Hz]"]);
    end
end
%% System Response to Wave and Wind Forcing (With Controller)
if WindWaves == 1    
    for iwindwaves = 1:1
    % If waves not activated
    if Waves == 0
        if Jonswap == 0
            % no waves, just water presence
            u = zeros(length(z_spar),length(t));
            udot = u;
            v = zeros(length(z_spar),length(t));
            vdot = v;
        elseif Jonswap == 1
            % JONSWAP wave spectrum for forcing
            [fvec,amp,~] = jonswap(Hs,Tp,df,fHighCut,gammaJS);
            [w,wdot] = IrregularWaveKinematics(fvec,amp);
            u = w*cosd(WaveAngle);
            udot = wdot*cosd(WaveAngle);
            v = w*sind(WaveAngle);
            vdot = wdot*sind(WaveAngle);
        elseif Jonswap == 2
            % linear wave forcing
            [w,wdot]=LinearWaveKinematics();
            u = w*cosd(WaveAngle);
            udot = wdot*cosd(WaveAngle);
            v = w*sind(WaveAngle);
            vdot = wdot*sind(WaveAngle);
        end
    end
    % create wind state, if not already activated
    if Wind == 0
        if Kaimal == 1
            [~,V_hub,~] = Kaimal_Timeseries(TI,TL,fHighCut);
        elseif Kaimal == 0
            V_hub = ones(size(t))*V_10;
        end
    end
    % forced response
    Y_windwaves = ode4(@dqdtsparbuoy,tode4,q0,M+A,B,C,4);
    % time series and frequency domain plotting
    figure
        PSD_single(tode4,Y_windwaves(:,1),fHighCut,wref(1),"Wind & Wave Forcing",["Surge [m]","PSD [m^2/Hz]"]);
        figure
        PSD_single(tode4,Y_windwaves(:,2),fHighCut,wref(1),"Wind & Wave Forcing",["Sway [m]","PSD [m^2/Hz]"]);
        figure
        PSD_single(tode4,Y_windwaves(:,3),fHighCut,wref(2),"Wind & Wave Forcing",["Heave [m]","PSD [m^2/Hz]"]);
        figure
        PSD_single(tode4,Y_windwaves(:,4),fHighCut,wref(3),"Wind & Wave Forcing",["Roll [deg]","PSD [deg^2/Hz]"]);
        figure
        PSD_single(tode4,Y_windwaves(:,5),fHighCut,wref(3),"Wind & Wave Forcing",["Pitch [deg]","PSD [deg^2/Hz]"]);
        figure
        PSD_single(tode4,Y_windwaves(:,6),fHighCut,wref(4),"Wind & Wave Forcing",["Yaw [deg]","PSD [deg^2/Hz]"]);    
    end
end