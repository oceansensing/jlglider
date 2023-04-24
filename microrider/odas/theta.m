function theta = theta(s,t,p,pr)

% THETA = THETA(S,T,P,PR) computes potential temperature at pr using
% Bryden 1973 polynomial for adiabatic lapse rate
% and Runge=Kutta 4 th order integration algorithm.
% Ref: Bryden,H.1973, DeepSea Research,20etc.
%
% Units:
%       p       pressur dBars
%       t       temp    deg C
%       s       salin   PSU
%       pr      ref press
%       theta   potent.temp deg C
%
% checkvalues: s (psu)  t (deg C)  p (dB)  pr (dB)  theta (deg C)
%               34.75      1.0      4500      0         0.640
%               34.75      1.0      4500     4000       0.944
%               34.95      2.5      3500      0         2.207
%               34.95      2.5      3500     4000       2.558
%
% Converted to a matlab function from a fortran function
% June 30, 1993 by TDM

h = pr-p;
xk = h.*atg(s,t,p);
t = t + 0.5*xk;
q = xk;
p = p + 0.5*h;
xk = h.*atg(s,t,p);
t = t + 0.29289322*(xk-q);
q = 0.58578644*xk + 0.121320344*q;
xk = h.*atg(s,t,p);
t = t + 1.707106781*(xk-q);
q = 3.414213562*xk - 4.121320344*q;
p = p + 0.5*h;
xk = h.*atg(s,t,p);
theta= t + (xk-2.0*q)/6.0;