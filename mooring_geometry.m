ld = 39.9;
betafairlead = 0.3702;
sin_beta = sin(betafairlead);
sin2_beta = sin_beta^2;
hd = roots([(1), (sin_beta*ld/2), (sin2_beta*(ld^2 - 50^2)/4)])