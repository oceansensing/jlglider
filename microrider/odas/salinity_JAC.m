%% salinity_JAC
% Convert conductivity and temperature data collected with a JAC_CT into
% salinity.
%%
% <latex>\index{Functions!salinity\_JAC}</latex>
%
%%% Syntax
%   [S, C_match] = salinity_JAC (P, T, C, process_info, ...)
%
% * [P]  Presure in units of dbar.
% * [T]  Temperature in degrees celsius.
% * [C]  Electrical conductivity in units of milli-siemens per centimetre.
% * [process_info] Structure containing the processing parameters. A
%           template is generated and returned when salinity_JAC is called
%           with no input parameters. The parameters are described below. 
% * [...] Optional configuration parameters to supplement, or override,
%           those values included within process_info.
% * []
% * [S] The salinity in practical salinity units using based on WHOI
%           fortran code. The default process_info structure is returned
%           when salinity_JAC is called with no input parameters.
% * [C_match] The conductivity singal that was matched to the response of
%                the temperature sensor and used to compute the salinity,
%                S. 
%
%%% Description
% This function matches the JAC_C singal to the JAC_T signal so that the
% computed salinity has reduced errors due to the inherent miss-match of
% all conductivity and temperature sensors. The errors are usually in the
% form of transient anomalies that bias the profiles of salinity, and
% result in erroneous density estimates.
%
% The function first applies a lag to the JAC_C signal. It then applies a
% single-pole low-pass filter to complete the matching of the conductivity
% data to the temperature data. 
%
% If $\texttt{salinity\_JAC}$ is called without input arguments, it returns
% the default parameters used by $\texttt{salinity\_JAC}$. For example,
%  
%    >> process_info = salinity_JAC
% 
%    process_info = 
%                   f_TC: 0.7300
%                     fs: 64
%                    lag: 0.0234
%                  speed: 0.6200
%
% 
%%% Parameters 
% The parameters within $\texttt{process\_info}$ are as follows:
%
% * [f_CT]  The cut-off (half-power) frequency of the low-pass filter
%            applied to the conductivity data, in units of Hz. 
%            Default = 0.73
% * [fs]   The sampling rate, in units of samples per second. 
%            Default = 64.
% * [lag]   The lag, or delay, applied to the conductivity data in units of
%            seconds. Default = 0.0234 
% * [speed] The speed of the CT sensor through the water, in units of m/s.
%            Default = 0.62. 
%
% If you specify either the f_CT or the lag parameters, then they are used
%  as specified. Otherwise, the default values are used and they will be
%  scaled by the speed. According to:
%
%       lag  = default_lag  *     (default_speed / speed). 
%       f_CT = default_f_CT * sqrt(speed / default_speed).
%

% * 2016-05-27 RGL Original version
% * 2016-05-28 RGL added C_match to output
% * 2020-06-11 JMM Assigned default structure to output in the case of no
%                  input parameters

function [S, C_match] = salinity_JAC( P, T, C, varargin )

% Default parameter values. 
default_f_TC  = 0.73; % low-pass filter cut-off frequency.
default_fs    = 64; % samples per second
default_lag   = 0.0234; % in seconds
default_speed = 0.62; % m/s


if ~nargin,
    for d = whos('default_*')',
        param = regexp(d.name, 'default_(.*)', 'tokens');
        result.(char(param{1})) = eval(d.name);
    end
    S = result;
    return
end

p = inputParser;
p.KeepUnmatched = true;
p.CaseSensitive = true;

val_numeric     = @(x) isnumeric(x) && isscalar(x) || isempty(x);
val_vector      = @(x) isnumeric(x) && isvector(x);

addRequired(  p, 'P', val_vector);
addRequired(  p, 'T', val_vector);
addRequired(  p, 'C', val_vector);
addParamValue(p, 'f_TC',  default_f_TC,  val_numeric);
addParamValue(p, 'fs',    default_fs,    val_numeric);
addParamValue(p, 'lag',   default_lag,   val_numeric);
addParamValue(p, 'speed', default_speed, val_numeric);


% Parse the arguments.
parse(p, P, T, C, varargin{:});

% Perform last stages of input validation.
if p.Results.f_TC <= 0,
  error('Invalid filter cut-off frequency, f_TC <= 0 .');
end
if p.Results.fs <= 0,
  error('Invalid  sampling rate, fs <= 0 .');
end
if p.Results.speed <= 0,
  error('Invalid speed, speed <= 0 .');
end
if p.Results.lag < 0,
  error('Invalid lag, lag < 0 .');
end


f_TC   = p.Results.f_TC;
fs     = p.Results.fs;
speed  = p.Results.speed;
lag    = p.Results.lag;

if any((strcmp(p.UsingDefaults, 'lag')))
    lag = lag *default_speed / speed;
end
if any((strcmp(p.UsingDefaults, 'f_TC')))
    f_TC = f_TC * sqrt(speed / default_speed);
end

if lag == 0 % nothing needs to be done
    C_match = C;
else % handle the lagging properly
    lag = lag * fs; % Convert to units of samples
    
    b = zeros(ceil(lag) + 1, 1); % size of lag filter
    b(end - 1) = ceil(lag) - lag ;
    b(end)     = 1 - b(end - 1) ; % make sum = 1
    
% calculate initial condition using the slope of the start of JAC_C and
%   lag the data using the filter function.
    
    L_b = length(b); % 1 + number of previous inputs required for lag-filter
    x = 0:L_b;
    x = x';
    p = polyfit(x, C(1: L_b + 1), 1); % linear polyfit
    previous_inputs = polyval(p, -x(2:end));
    z = filtic(b, 1, [], previous_inputs);
    C_match = filter(b, 1, C, z);
end

% Now apply match filter
[b_match,a_match] = butter(1, f_TC/(fs/2)); % coefficients for matching filter

% estimate initial conditions for the match filter
delay = 1/(2*pi*f_TC); % the delay induced by the match filter

x = 0:round(delay*fs) + 1; % length of delay in samples
x = x';
p = polyfit(x, C_match(1:length(x)), 1); % linear fit
initial_input  = polyval(p, -1); 
initial_output = polyval(p, -x(end));
z = filtic(b_match, a_match, initial_output, initial_input);

C_temp = filter(b_match, a_match, C_match, z); % apply matching low-pass filter

S = salinity (P, T, C_temp);




