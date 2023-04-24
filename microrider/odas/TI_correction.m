%% TI_correction
% Derive the corrected temperature of the fluid in a conductivity cell.
%%
% <latex>\index{Functions!TI\_correction}</latex>
%
%%% Syntax
%
%   T_cell = TI_correction( T, alpha, beta, fs )
%
% * [T]      Environmental temperature measured by a thermometer.
% * [alpha]  Amplitude coefficient for the thermal inertial of the
%            conductivity cell. 
% * [beta]   Inverse relaxation time for the thermal inertia of the
%            conductivity cell, in units of 1/s. 
% * [fs]     Sampling rate of the data, in samples per second.
% * []
% * [T_cell] Water temperature in the cell. This value should be used for
%            salinity calculations but not for density calculations.
%
%%% Description
% 
% Adjust the temperature measured with a thermometer to account
% for the temperature anomaly in a conductivity cell due to the thermal
% inertia of the cell. This function was originally written to process data
% collected with the Sea-Bird SBE3F and SBE4C sensors, but it can be used
% with any conductivity cell. This function is based on Lueck 
% (1990), Lueck and Picklo (1990) and Morrison et al (1994). Familiarity
% with these papers is nearly a pre-requisite for using this function. 
%
% The algorithm originally proposed by Lueck and Picklo (1990) tries to
% adjust the measured conductivity to account for the thermal inertial of the
% SBE4C cell. That is, it tries to adjust the measured conductivity to what
% it would be in the absence of thermal inertia. However, this is indirect
% and problematic because the temperature coefficient of conductivity:
%
% $$\left. \frac{\partial C}{\partial T}\right|_{P,S}$$ 
%
% is pressure-, salinity-, and temperature-dependent, and a single value
% can not be used for an entire profile. The algorithm proposed by Morrison
% et al (1994) is more direct and requires no further coefficients. They
% argue that the conductivity reported by the cell is the correct value for
% the fluid in the cell. The problem is that its temperature is different
% from that reported by the thermometer, and is wrong by the amount proposed
% by Lueck (1990). So a much simpler approach is to estimate the temperature
% in the cell and use it, with the measured conductivity and pressure, to
% calculate the salinity. The density is then calculated with this salinity,
% the real temperature in the water, and the pressure.
%
% The coefficients $\texttt{alpha}$ and $\texttt{beta}$ are pump-speed
% dependent. The original estimates by Lueck and Picklo (1990) have
% $\texttt{alpha}$ = 0.021, and $\texttt{beta}$ = 1/12 for an un-pumped
% system. A pumped system will have a larger $\texttt{beta}$. The user has
% to play with these coefficients to get the best correction but even poor
% values reduce significantly the errors in the calculated salinity. These
% errors are frequently called 'salinity spikes'.
%
% This function corrects only the thermal inertia. The short-term mismatch
% between the measured temperature and conductivity must be corrected
% before using this function. 
% 
% # Lueck, R.G., 1990, Thermal inertia of conductivity cells: Theory, J.
%   Atmos. Oceanic. Technol., 7, 741-768.
% # Lueck, R.G. and J.J. Picklo, 1990, Thermal inertia of conductivity cells:
%   Observations with a Sea-Bird cell, J. Atmos. Oceanic Technol., 7, 756-768.
% # Morrison, J., R. Anderson, N. Larson, E, D'Asaro and T. Boyd, 1994: The
%   correction for thermal-lag effects in Sea-Bird CTD data., J. Atmos. Oceanic
%   Technol., 11, 1151-1164.

% Version History
%
% * 2007-12-21 (RGL) original function
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-10-24 (WID) documentation update
% * 2015-11-02 (RGL) Documentation changes
% * 2015-11-19 (RGL) Documentation changes

function [T_cell]= TI_correction (T,alpha, beta, fs)

% dtime=1/0.6; % time step
% alpha=0.021; beta=1/12; % Lueck-Picklo coefficients
% a=2*alpha/(dtime*beta+2); % filter coefficients
% b=1-(2*a/alpha);

% develop filter coefficients
a = 2*alpha / (beta/fs + 2);   % coefficients for the filter function
b = 1 - 2*a/alpha;

z =   filtic([a -a], [1 b], 0, T(1));   % set up initial condition
T_cell = filter([a -a], [1 b], T, z);
T_cell = T - T_cell; 
