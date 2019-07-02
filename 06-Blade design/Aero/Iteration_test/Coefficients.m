function [Cl,Cd] = Coefficients(alpha)
%Coefficients:  We will interpolate the values of
%coefficients of lift and drag for a given alpha

%Data

CL = [0 -0.72 -0.65 -0.62 -0.6 -0.62 -0.65 -0.72 -0.92 -0.87 -0.765 -0.637 -0.488 -0.31 0 0.31 0.488 0.637 0.765 0.87 0.92 0.72 0.65 0.62 0.6 0.62 0.65 0.72 0];
CD = [0 -0.4 -0.345 -0.3 -0.26 -0.21 -0.17 -0.09 -0.049 -0.0335 -0.021 -0.015 -0.012 -0.011 0.013 0.011 0.012 0.015 0.021 0.0335 0.049 0.09 0.17 0.21 0.26 0.3 0.345 0.4 0];
ALPHA = [-180 -24 -22 -20 -18 -16 -14 -12 -11 -10 -8 -6 -4 -2 0 2 4 6 8 10 11 12 14 16 18 20 22 24 180];

Cl = interp1 (ALPHA, CL, alpha);
Cd = interp1 (ALPHA, CD, alpha);

end

