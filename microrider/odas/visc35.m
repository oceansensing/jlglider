%% visc35
% The kinematic viscosity of seawater for _S_ = 35.
%%
% <latex>\index{Functions!visc35}</latex>
%
%%% Syntax
%
%   v = visc35( T )
%
% * [T] temperature in units of degrees Celsius. 
% * []
% * [v] viscosity in units of metres-squared per second. 
%
%%% Description
%
% Return an approximation of the kinematic viscosity, $\nu$, in units of
% $\si{\square\m\per\s}$, based on temperature, in units of
% $\si{\celsius}$. The kinematic viscosity is derived from a
% $3^{\mathrm{rd}}$-order polynomial fit of $\nu$ against $T$ for
% salinity 35. The error of the approximation is less than
% $\SI{1}{\percent}$ for $0\le T \le \SI{20}{\celsius}$ and  $30\le S \le
% \SI{40}{PSU}$ at atmospheric pressure. 
%
% (see also viscosity)

% *Version History*
%
% * 1999-06-01 (FW)
% * 2002-10-09 (FW) revised
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-10-24 (WID) documentation update for publishing
% * 2015-11-02 (RGL) Documentation corrections.

function v = visc35(t)

pol=[-1.131311019739306e-011
      1.199552027472192e-009
     -5.864346822839289e-008
      1.828297985908266e-006];
v = polyval(pol,t);    
     
