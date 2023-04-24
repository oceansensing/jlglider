%% visc00
% The kinematic viscosity of freshwater,  _S_ = 0.
%%
% <latex>\index{Functions!visc00}</latex>
%
%%% Syntax
%
%   v = visc00( T )
%
% * [T] temperature in degrees Celsius 
% * []
% * [v] kinematic viscosity, in metres-squared per second 
%
%%% Description
%
% Returns an approximation of the kinematic viscosity, $\nu$, in units of
% $\si{\square\m\per\s}$, based on temperature, in units of
% $\si{\celsius}$. The kinematic viscosity is derived from a
% $3^{\mathrm{rd}}$-order polynomial fit of $\nu$ against $T$ for
% salinity 0. The error of the approximation is less than
% $\SI{1}{\percent}$ for $0\le T \le \SI{20}{\celsius}$ at atmospheric
% pressure. 
%
% (see also viscosity)

% Version History:
%
% * 1999-06-01 (FW)
% * 2002-10-09 (FW) revised
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-10-24 (WID) documentation update for publishing
% * 2015-11-02 (RGL) Documentation corrections.

function v = visc00(t)

pol=[-1.8377764060e-011
      1.4447664328e-009
     -6.0996917629e-008
      1.7903678481e-006];
v = polyval(pol,t);    
     
