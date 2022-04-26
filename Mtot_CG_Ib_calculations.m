function [m_tot,CG,Ibxx,Ibyy,Ibzz,Ib13] = Mtot_CG_Ib_calculations(M_nacellerotor,M_nacelle,M_rotor,M_spartower,z_hub,xCG_rotor,xCG_nacelle,zCG_nacelle,zCM_spartower)
%% Description
% This function calculates the body moment of inertias and total mass/cg of
% the sparbuoy and tower system
m_tot = M_nacellerotor + M_spartower;
CG_nr = [(xCG_rotor*M_rotor+xCG_nacelle*M_nacelle)/M_nacellerotor,0,(z_hub*M_rotor+zCG_nacelle*M_nacelle)/M_nacellerotor];
CG_st = [0,0,zCM_spartower];
%
xg = (M_nacellerotor*CG_nr(1) + M_spartower*CG_st(1))/(M_nacellerotor + M_spartower);
yg = (M_nacellerotor*CG_nr(2) + M_spartower*CG_st(2))/(M_nacellerotor + M_spartower);
zg = (M_nacellerotor*CG_nr(3) + M_spartower*CG_st(3))/(M_nacellerotor + M_spartower);
CG = [xg,yg,zg];
%
Ibxx_spartower = 1.5536*10^11;  % about CG
Ibyy_spartower = Ibxx_spartower;
Ibzz_spartower = 1.9025*10^9;
%
Ibxx_rotor = 3.493315913*10^8;    % about CG, calculated for 1 blade
Ibyy_rotor = Ibxx_rotor/2;  % approximating using perpendicular axis theorem
Ibzz_rotor = Ibyy_rotor;
%
Ibxx_nacelle = 12607277;
Ibyy_nacelle = 21433958;
Ibzz_nacelle = 18682468;
%
Ibxx = (Ibxx_spartower + M_spartower*(0-CG_st(3))^2) + (Ibxx_rotor + M_rotor*(0-z_hub)^2) + (Ibxx_nacelle + M_nacelle*(0-zCG_nacelle)^2);
%
Ibyy = (Ibyy_spartower + M_spartower*(0-CG_st(3))^2) + (Ibyy_rotor + M_rotor*(0-z_hub)^2) + (Ibyy_nacelle + M_nacelle*(0-zCG_nacelle)^2);
%
Ibzz = (Ibzz_spartower) + (Ibzz_rotor + M_rotor*(0-xCG_rotor)^2) + (Ibzz_nacelle + M_nacelle*(0-xCG_nacelle)^2);
%
Ib13 = -(M_rotor*xCG_rotor*z_hub) -(M_nacelle*xCG_nacelle*zCG_nacelle);