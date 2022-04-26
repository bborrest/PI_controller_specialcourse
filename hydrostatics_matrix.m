function [C_H] = hydrostatics_matrix(rho_H20,g,D_s,m_tot,CG,zb)
%% This function calculates the hydrostatics restoring
% See notes for calculations
I11 = 5875.2;
I22 = I11;
A = pi*(D_s/2)^2;
factor = 1;
%
H1j = [0,0,0,0,0,0];
H2j = [0,0,0,0,0,0];
H3j = [0,0,rho_H20*g*A,0,0,0];
H4j = [0,0,0,(rho_H20*g*I22 + m_tot*g*(zb-CG(3))),0,-m_tot*g*(0-CG(1))];
H5j = [0,0,0,0,(rho_H20*g*I11 + m_tot*g*(zb-CG(3))),-m_tot*g*(0-CG(2))];
H6j = [0,0,0,0,0,0];

C_H = [H1j;H2j;H3j;H4j;H5j;H6j];