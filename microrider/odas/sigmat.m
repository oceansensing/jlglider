function sigma = sigmat(t,s)

% SIGMA = SIGMAT(T,S) computes sigma-t knowing t (deg C), and s (o/oo).
% From Millero and  Poisson,DSR 28(6a),625 (1981)
%
% checkvalues: s (psu)  t (deg C)  sigma +1000 (kg m^-3)
%                0         5           999.96675
%                0         25          997.04796
%                35        5          1027.67547
%                35        25         1023.34306
%
% Modified May 2, 1988 by RGL
%
% Converted to a matlab function from a fortran subroutine
% June 30, 1993 by TDM
% Function checked with their check values.

a0      =  0.824493;
a1      = -4.0899e-03;
a2      =  7.6438e-05;
a3      = -8.2467e-07;
a4      =  5.3875e-09;

b0      = -5.72466e-03;
b1      =  1.02270e-04;
b2      = -1.65460e-06;

c       =  4.8314e-04;

r0      =   -0.157406;
r1      =    6.793952e-02;
r2      =   -9.095290e-03;
r3      =    1.001685e-04;
r4      =   -1.120083e-06;
r5      =    6.536336e-09;

sigma   = r0 + r1*t + r2*t.^2 + r3*t.^3 + r4*t.^4 + r5*t.^5 ...
        +(a0 + a1*t + a2*t.^2 + a3*t.^3 + a4*t.^4).*s ...
        +(b0 + b1*t + b2*t.^2).*sqrt(s.^3) ...
        + c*s.^2;
