function [M] = mass_matrix(m_tot,CG,Ibxx,Ibyy,Ibzz,Ib13)
%% This function creates the mass matrix, assuming Ib12=Ib21=0 and Ib23=Ib32=0
xg = CG(1);     % referenced to the water line
yg = CG(2);
zg = CG(3);
Ib31 = Ib13;
%
M1j = [m_tot,0,0,0,m_tot*zg,-m_tot*yg];
M2j = [0,m_tot,0,-m_tot*zg,0,m_tot*xg];
M3j = [0,0,m_tot,m_tot*yg,-m_tot*xg,0];
M4j = [0,-m_tot*zg,m_tot*yg,Ibyy,0,-Ib31];
M5j = [m_tot*zg,0,-m_tot*xg,0,Ibxx,0];
M6j = [-m_tot*yg,m_tot*xg,0,-Ib13,0,Ibzz];
%
M = [M1j;M2j;M3j;M4j;M5j;M6j];