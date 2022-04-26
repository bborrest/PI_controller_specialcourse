function k = wave_number(f,g,h)
% this function calculates the wave number from frequency (f) and depth (h)
% w = radian frequency
% 
% w=2*pi*f;

fun = @(k) (2*pi*f)^2 - g*k*tanh(k*h);

k=fzero(fun,[0 3]);
