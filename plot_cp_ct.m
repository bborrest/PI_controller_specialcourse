clear; clc;

load IEA_15MW_HWIND_Land_Based_hs2.mat
operating = readmatrix('operation_equinor.txt');

figure
contourf(turbine.pitchList,turbine.tsrList,max(turbine.Cp,0))
xlabel('Pitch [deg]')
ylabel('Tip-Speed Ratio [-]')

figure
contourf(turbine.pitchList,turbine.tsrList,max(turbine.Ct,0))
xlabel('Pitch [deg]')
ylabel('Tip-Speed Ratio [-]')

CPsurf = turbine.Cp;
CTsurf = turbine.Ct;
theta = turbine.pitchList;
tsr = turbine.tsrList;

% dCPdTh = (CPsurf(:,1:24) - CPsurf(:,2:25))./(theta(1:24) - theta(2:25));
% dCPdTSR = (CPsurf(1:55,:) - CPsurf(2:56,:))./(tsr(1:55) - tsr(2:56))';
% dCTdTh = (CTsurf(:,1:24) - CTsurf(:,2:25))./(theta(1:24) - theta(2:25));
% dCTdTSR = (CTsurf(1:55,:) - CTsurf(2:56,:))./(tsr(1:55) - tsr(2:56))';

del_tsr = tsr(2) - tsr(1)
del_theta = theta(2) - theta(1)
[testCPTh,testCPTSR] = gradient(CPsurf,del_tsr,del_theta)