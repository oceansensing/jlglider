%% salinity
% Compute salinity at a specified pressure, temperature and electrical conductivity
%%
% <latex>\index{Functions!salinity}</latex>
%
%%% Syntax
%
%   sal = salinity( P, T, C )
%
% * [P] Pressure [dbar]
% * [T] Temperature [C]
% * [C] Conductivity [mS/cm]
% * []
% * [sal] Salinity in practical salinity units [PSU]
%
%%% Description
% Compute salinity from vectors of pressure, temperature and electrical
% conductivity. 
%
% checkvalue: c = 42.914, s = 35 psu, t = 15.0 deg C , p = 0
%
% Based on Fortran routine developed at WHOI.

function sal = salinity(p,t,c)

% SAL = SALINITY(P,T,C) computes the salinity given the pressure (p),
% conductivity (c), and temperature (t).
%
% checkvalue: c = 42.914, s = 35 psu, t = 15.0 deg C , p = 0
%
% Converted to a matlab function from a WHOI fortran subroutine
% June 30, 1993 by TDM
% * 2015-11-01 RGL Updated documentation 
% * 2015-11-19 RGL Updated documentation.

aa1=2.07e-05; aa2=-6.37e-10; aa3=3.99989e-15;
bb1=3.426e-02; bb2=4.464e-04; bb3=.4215; bb4=-3.107e-03;
c0=6.766097e-01; c1=2.00564e-02; c2=1.104259e-04;
c3=-6.9698e-07; c4=1.0031e-09;
al0=.008; al1=-.1692; al2=25.3851; al3=14.0941;
al5=2.7081; bl0=.0005; bl1=-.0056; bl2=-.0066; bl3=-.0375;
bl4=.0636; bl5=-.0144; al4=-7.0261; k1 = .0162;

% The Lewis and Perkin Method

r = c/42.914;
r = r ./ ((c0+t.*(c1+t.*(c2+t.*(c3+t*c4)))) .* (1 + p.*(aa1+p.* ...
     (aa2+p.*aa3))./(1.+t.*(bb1+t*bb2)+bb3*r+bb4*t.*r)));
t = t - 15;
sal = al0 + (t*bl0)./(1 + k1*t) ...
    + (al1 + (t*bl1)./(1 + k1*t)).*sqrt(abs(r)) ...
    + (al2 + (t*bl2)./(1 + k1*t)).*r ...
    + (al3 + (t*bl3)./(1 + k1*t)).*sqrt(abs(r)).^3 ...
    + (al4 + (t*bl4)./(1 + k1*t)).*r.^2 ...
    + (al5 + (t*bl5)./(1 + k1*t)).*sqrt(abs(r)).^5;
