%% viscosity
% The kinematic and dynamic viscosity of sea-water
%%
% <latex>\index{Functions!viscosity}</latex>
%
%%% Syntax
%   [nu, mu] = viscosity( S, T, r )
%
% * [T] temperature in units of degrees Celsius
% * [S] salinity in practical salinity units
% * [r] density in units of kg per cubic-metre
% * []
% * [nu] the kinematic viscosity in units of metres-squared per second
% * [mu] the dynamic viscosity in units of kg per metre per second
%
%%% Description
%
% Calculates the kinematic viscosity, $\nu$ $[\si{\square\m\per\s}]$, and
% dynamic viscosity, $\mu$ $[\si{\kg\per\m\per\s}]$, of 
% sea-water following Millero (1974). The range of validity is $5 \le T \le
% \SI{25}{\celsius}$ and $0 \le S \le \SI{40}{PSU}$ at atmospheric pressure.
%
% Check values are: $\texttt{T = 25, S = 40, r(T,S): mu = 9.6541e-4}$.
%
% References:
%
% # Millero, J. F., 1974, The Sea, Vol 5, M. N. Hill, Ed, John Wiley, NY, p. 3.
% # Peters and Siedler, in Landolt-Bornstein New Series V/3a (Oceanography), pp 234.

% * 2015-11-02 RGL Documentation changes

function [nu, mu] = viscosity(s,t,r)
%
% Output parameters:
%    kinematic viscosity nu [m^2/s]
%    molecular viscosity mu [m^2/s]
%
% See also: visc0

% References: 
% 1) Millero, J. F., 1974, The Sea, Vol 5, M. N. Hill, Ed, John Wiley, NY, p. 3.
% 2) Peters and Siedler, in Landolt-Bornstein New Series V/3a (Oceanography), pp 234.
%

% Check value: t = 25, s = 40, r(t,s): mu = 9.6541e-4

% Fabian Wolk 1999/08/01
% (C) 2002 Rockland Oceanographic Services Inc.
% Author: Fabian Wolk
% Revision: 2002/07/05

error(nargchk(2,3,nargin));

% coefficients
a0=2.5116e-6; a1=1.2199e-6;
b0=1.4342e-6; b1=1.8267e-8;
c0=1.002e-3;
d0=1.1709; d1=1.827e-3; d2=89.93;

% molecular viscosity at s=0
mu0 = c0*10.^( (d0*(20-t) - d1*(t-20).^2) ./ (t+d2) );

% molecular visc at s
mu = mu0 .* (1 + (a0 + a1*t).*(r.*s).^0.5 + (b0 + b1*t).*r.*s );

% kinematic visc
nu = mu./r;
