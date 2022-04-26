function [A] = addedmass_matrix()
%% This function just returns the added mass matrix
a11 = 4.15752*10^7;
a22 = a11;
a51 = -3.21128*10^9;
a15 = a51;
a42 = -a51;
a24 = a42;
a44 = 3.25571*10^11;
a55 = a44;
a33 = 1.27039*10^6;
%
A1j = [a11,0,0,0,a15,0];
A2j = [0,a22,0,a24,0,0];
A3j = [0,0,a33,0,0,0];
A4j = [0,a42,0,a44,0,0];
A5j = [a51,0,0,0,a55,0];
A6j = [0,0,0,0,0,0];
%
A = [A1j;A2j;A3j;A4j;A5j;A6j];