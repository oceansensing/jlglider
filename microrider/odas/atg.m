function atg = atg(s,t,p)

% ATG = ATG(S,T,P) computes the adiabatic temperature gradient deg C/dBar
%
%       pressure        p       dBar
%       temperature     t       Deg C
%       salinity        s       PSU
%       adiabatic       atg     deg C/dBar
%
% Converted to a matlab function from a fortran subroutine
% June 30, 1993 by TDM

s = s - 35.0;
atg = (((-2.1687e-16 * t + 1.8676e-14) .* t - 4.6206e-13) .* p ...
    + ((2.7759e-12 * t - 1.1351e-10) .* s + ((-5.4481e-14 * t ...
    + 8.733e-12) .* t - 6.7795e-10) .* t + 1.8741e-8)) .* p ...
    + (-4.2393e-8 * t+1.8932e-6) .* s ...
    + ((6.6228e-10 * t - 6.836e-8) .* t + 8.5258e-6) .* t + 3.5803e-5;
