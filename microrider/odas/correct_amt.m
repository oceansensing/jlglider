%% correct_amt
% Adjust O2 data for temperature, salinity, and depth
%%
% <latex>\index{Functions!correct\_amt}</latex>
%
%%% Syntax
%  amt_O2 = correct_amt( O2, T, S, D, setupfilestr, name )
%
% * [O2] Disolved oxygen signal from AMT sensor.  Data must be of type 
%       'amt_o2' and converted using convert_odas before using this function.
% * [T] Temperature in degree celsius. Can be either empty, a scalar, or 
%       a vector.  If empty, it defaults to 10 degC.
% * [S] Salinity in units of PSU. Can be either empty, a scalar, or a
%       vector. If empty, it defaults to 0 PSU.
% * [D] Depth in units of metres. Can be either empty, a scalar, or a
%       vector.  If empty, it defaults to 0 m.
% * [setupfilestr] Configuration string or file name (with extension) used
%       for data acquisition. If data is recorded with a third-party
%       instrument then you must construct a configuation file according to
%       the specifications in the User Manual for your dissolved oxygen
%       measuring instrument.
% * [name] Section name, in the configuration-file or the
%       configuration-string, that contains the temperature coefficients
%       for the sensor. If not provided, the first section with type
%       'amt_o2' is selected.
% * []
% * [amt_O2] Oxygen data after being adjusted for temperature, salinity,
%       and depth.
%
%%% Description
% Apply the temperature-, salinity-, and depth-corrections to an AMT O2 signal
% to complete its conversion into physical units. The input is the AMT O2
% signal $\emph{after}$ conversion using the function convert_odas. The channel type
% must be 'amt_o2'.
%
% The temperature, salinity and depth information can each be supplied in
% one of three different forms. If you use an empty matrix, [ ], then the
% default value is applied to every data point of O2. If you specify a scalar
% (single) value, then this value is applied to every data point. If you
% supply a vector, then the vector is applied, point-for-point, to the
% oxygen signal, O2. If you supply a vector is must be the same size as O2. 
%
% Corrections are applied in the following manner.
%
% # Initial temperature adjustment by multiplying by ET.  The value of ET 
%   is provided by AMT.
%
% # Division by the factor (PN - Pw), where PN is a reference pressure 
%   (1013) and Pw is the water vapour saturtion pressure (dBar) calculated
%   from temperature.
%
% # Multiplication by the in situ O2 saturation concentration. The O2 
%   saturation concentration is calculated as a function of temperature and
%   salinity.
%
% # Adjust for depth my multiplying by "exp(-0.3775*Depth./T_K)" where
%   Depth is in meteres and T_K is in kelvin.
%
%%% Example
%
%    >> O2_tmp = convert_odas( AMT_O2, 'AMT_O2', setupfilestr )
%    >> AMT_O2 = correct_amt( O2_tmp, T1_fast, [], 15, setupfilestr )
%
% Before correct_amt is called, convert_odas is used to perform the initial
% stage of conversion into physical units.  correct_amt is then used to
% adjust for the specific temperature, salinity, and depth.  The above example
% demonstrates using a data vector for temperature (T1_fast), the default value for
% salinity ([ ]), and a constant value for depth (15).

% *Version History:*
%
% * 2015-05-07 RGL original version.
% * 2015-07-24 WID lots of formatting changes - basic algorithm unchanged.
% * 2012-10-27 (RGL) Minor documentation corrections.
% * 2015-10-28 (RGL) Documentation corrections.


function amt_O2_corrected = correct_amt( O2, T, S, D, setupfilestr, name )

if isempty(O2), error ('Empty oxygen signal'), end
if isempty(T), T = 10; end % default value
if isempty(S), S = 0;  end % default value
if isempty(D), D = 0;  end % default value

% Anyone who want this function to run fast will ensure that interpolation
% is not required.
T = interp1_if_required(T, O2);
S = interp1_if_required(S, O2);
D = interp1_if_required(D, O2);

% Work when setupfilestr is in a cell.
if iscell(setupfilestr), setupfilestr = setupfilestr{1}; end

% If setupfilestr is a string and points to a file, load and use as the
% configuration string.
if ischar(setupfilestr) && exist( setupfilestr, 'file' )
    setupfilestr = fileread( setupfilestr );
end

% If setupfilestr is a string, assume it is the configuration string.
% Otherwise it is a pre-compiled object that can be used directly.
if ischar(setupfilestr)
    obj = setupstr(setupfilestr);
else
    obj = setupfilestr;
end

if nargin == 5
    sec  = setupstr(obj, '', 'type', 'amt_o2');
    name = setupstr(obj, sec, 'name');
    if isempty(sec) || isempty(name)
        error('Unable to find channel of type amt_o2.');
    end
end

sec = setupstr(obj, name);
if isempty(sec)
    error('Unable to find channel with name: %s', char(name));
end
    
b0 = setupstr(obj, name, 'b0');   b0 = str2double(b0{1});
b1 = setupstr(obj, name, 'b1');   b1 = str2double(b1{1});
b2 = setupstr(obj, name, 'b2');   b2 = str2double(b2{1});
b3 = setupstr(obj, name, 'b3');   b3 = str2double(b3{1});

b = [ b3 b2 b1 b0 ];
if size(b,2) < 4
    error('Missing b-parameters in configuration string')
end

ET = polyval(b, T);
T_K = T + 273.15; % In degree kelvin

amt_O2_corrected = O2 .* ET; % the basic temperature correction

P_w = get_wvP(T_K); % water vapour pressure, mbar
P_N = 1013; % standrd pressure, mbar
amt_O2_corrected = amt_O2_corrected ./ (P_N - P_w);

% Calculate saturation concentration
C_T   = get_O2_sat_T(T_K);
C_Sal = get_O2_sat_S(T_K, S);
C_S = 1.4289*C_T .* C_Sal;
amt_O2_corrected = amt_O2_corrected .* C_S;

% Adjust for depth
amt_O2_corrected = amt_O2_corrected .* exp(-0.3775.* D ./ T_K);

end



function value = get_O2_sat_S(T,S)
%
% value = get_O2_sat_S(T,S)
%
% function to calculate the saturation concentration of O2, salinity-
% dependence. T in kelvin, S in grams/kilo-grams (PSU).
% In support of AMT disolved O2 sensor.
%
% 2015-05-07, RGL, original version.

P = [-0.033096  0.014259 -0.0017];
P = fliplr(P);
value = polyval(P, T/100);
value = value .* S;
value = exp(value);
end



function value = get_O2_sat_T(T)
%
%  value = get_O2_sat_T(T)
%
% function to calculate the saturation concentration of O2, T-dependence
% only. T in kelvin.
% In support of AMT disolved O2 sensor
%
% 2015-05-07, RGL, original version

value = -173.4292 + ...
         249.6339*(100 ./ T) + ...
         143.3483*log(T / 100) - ...
          21.8492*T / 100;
value = exp(value);
end



function value = get_wvP(T)
% 
% value = get_wvP(T)
% function to calculate the water vapour pressure, in mbar. Temperature, T,
% in % kelvin.
% In support of AMT disolved O2 sensor.
%
% 2015-05-07, RGL, orginal version.

P = [-2.16961e5 -3.8407e3 11.8571];
value = polyval(P, 1 ./ T);
value = 1e3*exp(value);
end

function X = interp1_if_required( slave, master )
% Perform an interpolation of 'slave' to ensure it has the same number of
% points as 'master'.  If an interpolation is not required, return slave
% unchanged.

X = slave;
if length(slave) == length(master)
    return
elseif length(slave) == 1
    slave = [slave; slave];
end

t_slave = (0:length(slave)-1)' .* (length(master)/length(slave));
t_master = (0:length(master)-1)';

X = interp1(t_slave, slave, t_master, 'pchip');
end
