function [Fx,Fz,ls] = inelastic_catenary_line(h,X,l,w)

%% function [Fx,Fz,ls] = inelastic_catenary_line(h,X,l,w)
% Solution to the inelastic catenary line equations
% pp.260 of Faltinsen (1990), 'Sea loads on ships and offshore structures',
% Cambridge University Press.
%
% DTU Wind Energy, August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INPUTS %%%%
% h             Vertical distance between anchor and fairlead       [m]
% X             Horizontal distance between anchor and fairlead     [m]
% l             Un-stretched line length                            [m]
% w             Weight per unit length of line submerged in water [N/m]
%
%%%% OUTPUTS %%%%
% Fx            Horizontal line force at fairlead                   [N]
% Fz            Vertical line force at fairlead                     [N]
% ls            Suspended line length                               [m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if geometric limits are exceeded
if l<sqrt(h^2+X^2)
    error('Line is too short')
end

if l>(h+X)
    error('Line is too long')
end

% Initial guess of 'a' quantity, ref. Eq.(8.21) in Faltinsen (1990)
if l>20    	% Full-scale mooring line
    X0 = 1e3;                 
else 		% Lab-scale mooring line
    X0 = 5; 
end

% Solving catenary shape equation iteratively
fun = @(a) (l-h*(1+2*a/h)^(1/2)+a*acosh(1+h/a))-X;
options = optimset('display','off');
[a, fval, exitflag, output] = fzero(fun,X0,options);
while isnan(a)
    X0 = X0/2;
    [a, fval, exitflag, output] = fzero(fun,X0,options);
end
    

% Calculation of horizontal fairlead force, Fx
Fx = w*a;

% Calculation of suspending line length, ls
ls = sqrt(h^2+2*h*a);

% Calculation of vertical fairlead force, Fz
Fz = w*ls;

return;