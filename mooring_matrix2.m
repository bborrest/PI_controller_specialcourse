function [C_M] = mooring_matrix2(X,hmoor,water_depth,K66)
%% Mooring stiffness calculations
% This function calculates the mooring stiffness matrix based on catenary
% equations for a stiff catenary line resting on the seabed. For the 
% yaw-yaw case, the parameter is adjusted to match natural frequencies. The
% inertial coordinate system is set to the mean sea level at the axial
% center of the spar-buoy with zero displacement.
%
% This code can be adapted for the number of mooring lines, n, attached by
% a single fair-lead point to the spar. The first mooring line is attached
% at assumed angle beta = 0, which is set to be poiting upwind. This offset
% can be adjusted manually.
% Equations and work taken from "Stiffness of Slack Mooring Lines" by 
% Mohammed Khair Al-Solihat & Meyer Nahon in Ships and Offshore Structures, 
% 2016.
%% Inputs
% X, [x1, x2, x3, x4, x5, x6], the linearized displacement coordinates,
% with x4, x5, x6 in radians!!
% K66, manually adjusted yaw-yaw stiffness parameter for split fairleads
%% Mooring data
global rho_H2O
ld = 39.1;                                      % delta lines projected length (43.9) (39.1)
hd = 3.1104;                                      % delta lines projected height of connection adjusted, positive down (2.48) (3.1104)
Zm = (hmoor - hd) - water_depth;                % mooring fairlead attachment Z-coordinate in body coordinates [m]
L = 565 + ld;                                   % line length
Rc = 600;                                       % anchor radius [m]
Rs = 9.3;                                       % sparbuoy radius at fairlead connection
wmoor = 561.25*(1 - (rho_H2O*pi*0.08^2)/561.25)*9.81;   % N/m
%% Coordinate Geometry
% fix below for H,V,l being 3x1
n = 3;          % number of mooring lines
beta_1 = 0;     % offset angle for the first mooring cable (in radians!!)
% mooring line angles (radians!!)
beta = zeros(n,1);
for i = 1:n
    beta(i) = beta_1 + (2*pi/n)*(i-1);
end
phi = X(4);     % euler angles
theta = X(5);
psi = X(6);
% rotation matrix
R(1,1) = cos(psi)*cos(theta);
R(1,2) = cos(psi)*sin(theta)*sin(phi) - sin(psi)*cos(phi);
R(1,3) = cos(psi)*sin(theta)*cos(phi) + sin(psi)*sin(phi);
R(2,1) = sin(psi)*cos(theta);
R(2,2) = sin(psi)*sin(theta)*sin(phi) + cos(psi)*cos(phi);
R(2,3) = sin(psi)*sin(theta)*cos(phi) - cos(psi)*sin(phi);
R(3,1) = -sin(theta);
R(3,2) = cos(theta)*sin(phi);
R(3,3) = cos(theta)*cos(phi);
r_po = [Rs*cos(beta)';Rs*sin(beta)';Zm*ones(size(beta'))]; % position fairlead location in body system (X-row, Y-row, Z-row)
rm = X(1:3) + R*r_po; % rm vector defined in the inertial coordinate system
% anchor coordinates (global coordinates)
XA = Rc*cos(beta)'; % x-locations of the anchors
YA = Rc*sin(beta)'; % y-locations of the anchors
ZA = -water_depth*ones(size(beta')); % z-locations of the anchors
% fairlead coordinates (global coordinates)
Xp = rm(1,:);
Yp = rm(2,:);
Zp = rm(3,:);
% body frame center coordinates in inertial frame
rX = X(1);
rY = X(2);
rZ = X(3);
% l and h
l = (sqrt((Xp - XA).^2 + (Yp - YA).^2))';
h = (Zp - ZA)';
%% Forces!!
H = zeros(size(beta));
V = zeros(size(H));
for i = 1:n
    [H(i),V(i),~] = inelastic_catenary_line(h(i),l(i),L,wmoor);
end
% mooring geometry check
betafairlead = atan(V./H);
%% Kinematics Differentials
% R (rotation matrix)- rotational partial derivatives
dR_dphi(1,1) = 0;
dR_dphi(1,2) = cos(psi)*sin(theta)*cos(phi) + sin(psi)*sin(phi);
dR_dphi(1,3) = - cos(psi)*sin(theta)*sin(phi) + sin(psi)*cos(phi);
dR_dphi(2,1) = 0;
dR_dphi(2,2) = sin(psi)*sin(theta)*cos(phi) - cos(psi)*sin(phi);
dR_dphi(2,3) = -sin(psi)*sin(theta)*sin(phi) - cos(psi)*cos(phi);
dR_dphi(3,1) = 0;
dR_dphi(3,2) = cos(theta)*cos(phi);
dR_dphi(3,3) = -cos(theta)*sin(phi);
%
dR_dtheta(1,1) = -cos(psi)*sin(theta);
dR_dtheta(1,2) = cos(psi)*cos(theta)*sin(phi);
dR_dtheta(1,3) = cos(psi)*cos(theta)*cos(phi);
dR_dtheta(2,1) = -sin(psi)*sin(theta);
dR_dtheta(2,2) = sin(psi)*cos(theta)*sin(phi);
dR_dtheta(2,3) = sin(psi)*cos(theta)*cos(phi);
dR_dtheta(3,1) = -cos(theta);
dR_dtheta(3,2) = -sin(theta)*sin(phi);
dR_dtheta(3,3) = -sin(theta)*cos(phi);
%
dR_dpsi(1,1) = -sin(psi)*cos(theta);
dR_dpsi(1,2) = -sin(psi)*sin(theta)*sin(phi) - cos(psi)*cos(phi);
dR_dpsi(1,3) = -sin(psi)*sin(theta)*cos(phi) + cos(psi)*sin(phi);
dR_dpsi(2,1) = cos(psi)*cos(theta);
dR_dpsi(2,2) = cos(psi)*sin(theta)*sin(phi) - sin(psi)*cos(phi);
dR_dpsi(2,3) = cos(psi)*sin(theta)*cos(phi) + sin(psi)*sin(phi);
dR_dpsi(3,1) = 0;
dR_dpsi(3,2) = 0;
dR_dpsi(3,3) = 0;
% rm (fairlead attachment location)- position partial derivatives
drm_dphi = dR_dphi*r_po;
dXp_dphi = drm_dphi(1,:)';
dYp_dphi = drm_dphi(2,:)';
dZp_dphi = drm_dphi(3,:)';
%
drm_dtheta = dR_dtheta*r_po;
dXp_dtheta = drm_dtheta(1,:)';
dYp_dtheta = drm_dtheta(2,:)';
dZp_dtheta = drm_dtheta(3,:)';
%
drm_dpsi = dR_dpsi*r_po;
dXp_dpsi = drm_dpsi(1,:)';
dYp_dpsi = drm_dpsi(2,:)';
dZp_dpsi = drm_dpsi(3,:)';
% l (horizontal cable projection)- partial derivatives
dl_drx = cos(beta);
dl_dry = sin(beta);
dl_drz = 0;
dl_dphi = cos(beta).*dXp_dphi + sin(beta).*dYp_dphi;
dl_dtheta = cos(beta).*dXp_dtheta + sin(beta).*dYp_dtheta;
dl_dpsi = cos(beta).*dXp_dpsi + sin(beta).*dYp_dpsi;
% beta- partial derivatives
dbeta_drx = -sin(beta)./l;
dbeta_dry = cos(beta)./l;
dbeta_drz = 0;
dbeta_dphi = (1./l).*(cos(beta).*dYp_dphi - sin(beta).*dXp_dphi);
dbeta_dtheta = (1./l).*(cos(beta).*dYp_dtheta - sin(beta).*dXp_dtheta);
dbeta_dpsi = (1./l).*(cos(beta).*dYp_dpsi - sin(beta).*dXp_dpsi);
% h (vertical cable projection)- partial derivatives
dh_drx = 0;
dh_dry = 0;
dh_drz = 1;
dh_dphi = dZp_dphi;
dh_dtheta = dZp_dtheta;
dh_dpsi = dZp_dpsi;
%% Two-Dimensional Stiffness Matrix of a single mooring line
K_P11 = 1./((1/wmoor)*(-V./sqrt(H.^2 + V.^2) + asinh(V./H)));
K_P12 = 1./((1/wmoor)*(H./sqrt(H.^2 + V.^2) - 1));
K_P21 = K_P12;
K_P22 = 1./((V/wmoor).*(1./sqrt(H.^2 + V.^2)));
%% Exact generic mooring stiffness matrix coefficients, for n mooring cables
% First row, K1j
K11 = K_P11.*(cos(beta)).^2 + (H./l).*(sin(beta)).^2;
K12 = sin(beta).*cos(beta).*(K_P11 - H./l);
K13 = cos(beta).*K_P12;
K14 = cos(beta).*(K_P11.*dl_dphi + K_P12.*dh_dphi) - H.*sin(beta).*dbeta_dphi;
K15 = cos(beta).*(K_P11.*dl_dtheta + K_P12.*dh_dtheta) - H.*sin(beta).*dbeta_dtheta;
K16 = 0;
C1j = [sum(K11), sum(K12), sum(K13), sum(K14), sum(K15), sum(K16)];
% Second row, K2j
K21 = K12;
K22 = (sin(beta)).^2.*K_P11 + (cos(beta)).^2.*(H./l);
K23 = sin(beta).*K_P12;
K24 = sin(beta).*(K_P11.*dl_dphi + K_P12.*dh_dphi) + H.*cos(beta).*dbeta_dphi;
K25 = sin(beta).*(K_P11.*dl_dtheta + K_P12.*dh_dtheta) + H.*cos(beta).*dbeta_dtheta;
K26 = 0;
C2j = [sum(K21), sum(K22), sum(K23), sum(K24), sum(K25), sum(K26)];
% Third row, K3j
K31 = K13;
K32 = K23;
K33 = K_P22;
K34 = K_P12.*dl_dphi + K_P22.*dh_dphi;
K35 = K_P21.*dl_dtheta + K_P22.*dh_dtheta;  % why is this value nonzero?
K36 = 0;
C3j = [sum(K31), sum(K32), sum(K33), sum(K34), sum(K35), sum(K36)];
% Fourth row, K4j
K41 = (Yp' - rY).*K31 - (Zp' - rZ).*K21;
K42 = (Yp' - rY).*K32 - (Zp' - rZ).*K22;
K43 = (Yp' - rY).*K33 - (Zp' - rZ).*K23;
K44 = (Yp' - rY).*K34 - (Zp' - rZ).*K24 + V.*dYp_dphi - H.*sin(beta).*dZp_dphi;
K45 = (Yp' - rY).*K35 - (Zp' - rZ).*K25 + V.*dYp_dtheta - H.*sin(beta).*dZp_dtheta;
K46 = 0;
C4j = [sum(K41), sum(K42), sum(K43), sum(K44), sum(K45), sum(K46)];
% Fifth row, K5j
K51 = (Zp' - rZ).*K11 - (Xp' - rX).*K31;
K52 = (Zp' - rZ).*K12 - (Xp' - rX).*K32;
K53 = (Zp' - rZ).*K13 - (Xp' - rX).*K33;
K54 = (Zp' - rZ).*K14 - (Xp' - rX).*K34 + H.*cos(beta).*dZp_dphi - V.*dXp_dphi;
K55 = (Zp' - rZ).*K15 - (Xp' - rX).*K35 + H.*cos(beta).*dZp_dtheta - V.*dXp_dtheta;
K56 = 0;
C5j = [sum(K51), sum(K52), sum(K53), sum(K54), sum(K55), sum(K56)];
% K66
K61 = 0;
K62 = 0;
K63 = 0;
K64 = 0;
K65 = 0;
K66 = K66;
C6j = [sum(K61), sum(K62), sum(K63), sum(K64), sum(K65), sum(K66)];
% Total Matrix
C_M = [C1j;C2j;C3j;C4j;C5j;C6j];