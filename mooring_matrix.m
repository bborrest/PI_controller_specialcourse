function [C_M] = mooring_matrix(X,h,a,w,Zm,phi,moorxy,x1,x2,x3,x4,x5,x6,TH,TV,yawyaw)
%% Mooring stiffness calculations
% This function calculates the mooring stiffness based on catenary
% equations. For the yaw-yaw case, the parameter is adjusted to match
% natural frequencies.
%% Inputs
% X,    (3x1) the 0-point position sea-bed distance for each mooring line[m]*** leads the dimension of TH
% L,    (1x1) the total line length of each line [m]
% h,    (1x1) the 0-point position vertical fair-lead distance [m]
% a,    (3x1) slope of TH/w *** follows the dimension of X
% w,    (1x1) "wet" weight per length of iron chain [kg/m]
% Zm,   (1x1) the vertical distance between the rotational center and fair-lead connection [m]
% phi,  (3x1) angle between mooring line and the coordinate system [rad]
% moorxy, (3x2) the X-Y coordinates of each mooring anchor location relative to the spar [m]
% x1,   (1x1) the current x-location of the spar-buoy relative to equilibrium [m]
% x2,   (1x1) the current y-location of the spar-buoy relative to equilibrium [m]
% x4,   (1x1) the current roll angle [rad]
% x5,   (1x1) the current pitch angle [rad]
% TH,   (3x1) the horizontal tension of each mooring line [N]*** follows the dimension of X
% TV,   (3x1) the vertical tension of each mooring line [N]*** follows the dimension of X
% yawyaw, (1x1) the K66 term yaw stiffness
%% Coordinate Geometry
% these terms make some simplifying assumptions about geometry: they treat
% each of the fairlead connection coordinates as separate points but
% behaving in the same way, so there is no coordinate transformation for
% 100% accurate updates.
% *** later update these to use rigid body coordinate transformation
% why are these terms always positive??
R = 9.3;        % spar radius
Xf = x1*cos(phi) + x5*Zm*cos(phi) + R*cos(phi + x3);   % fairlead connections x-coordinates (3x1)
Yf = x2*sin(phi) - x4*Zm*sin(phi) + R*sin(phi + x3);   % fairlead connections y-coordinates (3x1)
Zf = [1;1;1]*Zm*(3 - cos(x4) - cos(x5)) + x3;   % fairlead connections z-coordinates (3x1)
%% Base partial derivatives
dTHda = w;  %(1x1)
dTVda = w*h./sqrt(h^2 + 2*h*a); %(3x1)
dXda = acosh(1 + h./a) - 2./sqrt(1 + 2*a/h);    %(3x1)
dadX = 1./dXda; %(3x1)
dXdh = 1./(sqrt(h./a).*sqrt(2 + h./a)) - (a + h)./(h*sqrt(1 + 2*a/h));   %(3x1)
dhdX = 1./dXdh; %(3x1)
%% Common partial derivatives
dTHdX = dTHda*dadX;	%(3x1)
dTHdh = dTHda*dadX.*dXdh;	%(3x1)***WRONG SIGN
dTHdh = - dTHdh;
dTVdX = dTVda.*dadX;	%(3x1)
dTVdh = w*(h + a)./sqrt(h^2 + 2*h*a);   %(3x1)
%% Geometry partial derivatives -> check moorxy terms in Xf, Yf, Zf need to be added
% x1
dXdx1 = cos(phi); %(3x1)
dcosdx1 = -1./X - cos(phi)./X.^2.*(moorxy(:,1) - x1);   %(3x1)
dsindx1 = -cos(phi)./X.^2.*(moorxy(:,2) - x2);  %(3x1)
dZmdx1 = dhdX.*dXdx1;   %(3x1)
dphidx1 = -sin(phi)./X; % (3x1)
% x2
dXdx2 = sin(phi);	%(3x1)
dcosdx2 = -sin(phi)./X.^2.*(moorxy(:,1) - x1);	%(3x1)
dsindx2 = -1./X - sin(phi)./X.^2.*(moorxy(:,2) - x2);	%(3x1)
dZmdx2 = dhdX.*dXdx2;   %(3x1)
dphidx2 = cos(phi)./X;  % (3x1)
% x3
dcosdx3 = -1./(X.^2).*dXdh.*(moorxy(:,1) - x1);	%(3x1)
dsindx3 = -1./(X.^2).*dXdh.*(moorxy(:,2) - x2); %(3x1)
dZmdx3 = -1;	%(1x1)
% x4
dx2dx4 = -Zm*cos(x4); %(1x1)
dhdx4 = -Zm*sin(x4);  %(1x1)
dXfdx4 = 0;	%(1x1)
dYfdx4 = -Zm*sin(phi);  %(3x1)
dZfdx4 = Zm*sin(x4);    %(1x1)
% x5
dx1dx5 = Zm*cos(x5);	%(1x1)
dhdx5 = -Zm*sin(x5);	%(1x1)
dXfdx5 = Zm*cos(phi);   %(3x1)
dYfdx5 = 0; %(1x1)
dZfdx5 = Zm*sin(x5);    %(1x1)
%% x1 derivatives
% Critical Term
dF1dx1 = sum(dTHdX.*dXdx1.*cos(phi) - TH.*sin(phi).*(-sin(phi)./X)); %(1x1)
% ^^^
dF2dx1 = sum(dTHdX.*dXdx1.*sin(phi) + TH.*cos(phi).*(-sin(phi)./X)); %(1x1)
dF3dx1 = sum(dTVdX.*dXdx1); %(1x1)
dT1dx1 = sum(-Zf.*(dTHdX.*dXdx1.*sin(phi) + TH.*cos(phi).*(-sin(phi)./X)) + Yf.*(dTVdX.*dXdx1)); %(1x1)
% Critical Term
dT2dx1 = sum(Zf.*(dTHdX.*dXdx1.*cos(phi) - TH.*sin(phi).*(-sin(phi)./X)) - Xf.*(dTVdX.*dXdx1)); %(1x1)
% ^^^
dT3dx1 = 0;
M1j = [dF1dx1;dF2dx1;dF3dx1;dT1dx1;dT2dx1;dT3dx1]; %(6x1)
%% x2 derivatives
dF1dx2 = sum(dTHdX.*dXdx2.*cos(phi) - TH.*sin(phi).*(cos(phi)./X)); %(1x1)
% Critical Term
dF2dx2 = sum(dTHdX.*dXdx2.*sin(phi) + TH.*cos(phi).*(cos(phi)./X)); %(1x1)
% ^^^
dF3dx2 = sum(dTVdX.*dXdx2); %(1x1)
% Critical Term
dT1dx2 = sum(-Zf.*(dTHdX.*dXdx2.*sin(phi) + TH.*cos(phi).*(cos(phi)./X)) + Yf.*(dTVdX.*dXdx2)); %(1x1)
% ^^^
dT2dx2 = sum(Zf.*(dTHdX.*dXdx2.*cos(phi) - TH.*sin(phi).*(cos(phi)./X)) - Xf.*(dTVdX.*dXdx2)); %(1x1)
dT3dx2 = 0;
M2j = [dF1dx2;dF2dx2;dF3dx2;dT1dx2;dT2dx2;dT3dx2]; %(6x1)
%% x3 derivatives
dF1dx3 = sum(dTHdh.*cos(phi)); %(1x1)
dF2dx3 = sum(dTHdh.*sin(phi)); %(1x1)
% Critical Term
dF3dx3 = sum(dTVdh); %(1x1)
% ^^^
dT1dx3 = sum(-Zf.*(dTHdh.*sin(phi)) + Yf.*(dTVdh)); %(1x1)
dT2dx3 = sum(Zf.*(dTHdh.*cos(phi)) - Xf.*(dTVdh)); %(1x1)
dT3dx3 = 0;
M3j = [dF1dx3;dF2dx3;dF3dx3;dT1dx3;dT2dx3;dT3dx3]; %(6x1)
%% x4 derivatives
dF1dx4 = sum((dTHdX.*dXdx2*dx2dx4 + dTHdh.*dhdx4).*cos(phi) - TH.*sin(phi).*(cos(phi)./X)*dx2dx4); %(1x1)
% Critical Term
dF2dx4 = sum((dTHdX.*dXdx2*dx2dx4 + dTHdh.*dhdx4).*sin(phi) + TH.*cos(phi).*(cos(phi)./X)*dx2dx4); %(1x1)
% ^^^
dF3dx4 = sum(dTVdX.*dXdx2*dx2dx4 + dTVdh*dhdx4); %(1x1)
% Critical Term
dT1dx4 = sum(-Zf.*((dTHdX.*dXdx2*dx2dx4 + dTHdh.*dhdx4).*sin(phi) + TH.*cos(phi).*(cos(phi)./X)*dx2dx4) - dZfdx4.*TH.*sin(phi) + Yf.*(dTVdX.*dXdx2*dx2dx4 + dTVdh*dhdx4) + dYfdx4.*TV.*sin(phi)); %(1x1)
% ^^^
dT2dx4 = sum(Zf.*((dTHdX.*dXdx2*dx2dx4 + dTHdh.*dhdx4).*sin(phi) + TH.*cos(phi).*(cos(phi)./X)*dx2dx4) + dZfdx4.*TH.*cos(phi) - Xf.*(dTVdX.*dXdx2*dx2dx4 + dTVdh*dhdx4) - dXfdx4.*TV.*cos(phi)); %(1x1) %(1x1)
dT3dx4 = 0;
M4j = [dF1dx4;dF2dx4;dF3dx4;dT1dx4;dT2dx4;dT3dx4]; %(6x1)
%% x5 derivatives
% Critical Term
dF1dx5 = sum((dTHdX.*dXdx1*dx1dx5 + dTHdh.*dhdx5).*cos(phi) - TH.*sin(phi).*(-sin(phi)./X)*dx1dx5); %(1x1)
% ^^^
dF2dx5 = sum((dTHdX.*dXdx1*dx1dx5 + dTHdh.*dhdx5).*sin(phi) + TH.*cos(phi).*(-sin(phi)./X)*dx1dx5); %(1x1)
dF3dx5 = sum(dTVdX.*dXdx1*dx1dx5 + dTVdh.*dhdx5); %(1x1)
dT1dx5 = sum(-Zf.*((dTHdX.*dXdx1*dx1dx5 + dTHdh.*dhdx5).*sin(phi) + TH.*cos(phi).*(-sin(phi)./X)*dx1dx5) - dZfdx5.*TH.*sin(phi) + Yf.*(dTVdX.*dXdx1*dx1dx5 + dTVdh.*dhdx5) + dYfdx5.*TV.*sin(phi)); %(1x1)
% Critical Term
dT2dx5 = sum(Zf.*((dTHdX.*dXdx1*dx1dx5 + dTHdh.*dhdx5).*cos(phi) - TH.*sin(phi).*(-sin(phi)./X)*dx1dx5) + dZfdx5.*TH.*cos(phi) - Xf.*(dTVdX.*dXdx1*dx1dx5 + dTVdh.*dhdx5) - dXfdx5.*TV.*cos(phi)); %(1x1)
% ^^^
dT3dx5 = 0;
M5j = [dF1dx5;dF2dx5;dF3dx5;dT1dx5;dT2dx5;dT3dx5]; %(6x1)
%% x6 derivatives/values
M6j = [0;0;0;0;0;yawyaw];
%% Cosntructing mooring stiffness matrix
C_M = [M1j,M2j,M3j,M4j,M5j,M6j];