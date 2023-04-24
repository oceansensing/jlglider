%% odas_lag
% Lag a signal by an amount specified in seconds (and not samples)
%%
% <latex>\index{Functions!odas\_lag}</latex>
%
%%% Syntax
%   [x_lag] = odas_lag(x, lag, fs)
%
% * [x]     The signal to be lagged.
% * [lag]   The lag (or advance) in units of seconds.
% * [fs]    The sampling rate in units of samples per second.
% * []
% * [x_lag] The lagged version of input x
%
%%% Description
% This function applies a lag to the input signal, x. Thd lag can be
% negative for an advance.
%

% * 2016-06-23 RGL Original version

function x_lag = odas_lag( x, lag, fs )


p = inputParser;
p.KeepUnmatched = true;
p.CaseSensitive = true;

val_numeric     = @(z) isnumeric(z) && isscalar(z) || isempty(z);
val_vector      = @(z) isnumeric(z) && isvector(z);

addRequired(  p, 'x',   val_vector);
addRequired(  p, 'lag', val_numeric);
addRequired(  p, 'fs',  val_numeric);

% Perform last stages of input validation.
if fs <= 0,
  error('Invalid sampling rate, fs must be > 0.');
end

x= x(:);  % make it a column vector
flip_me = false;

if lag == 0 % nothing needs to be done
    x_lag = x;
    return
end

% handle the lagging properly
lag = lag * fs; % Convert to units of samples
if lag < 0
    flip_me = true;
    x = flipud(x);
end

b = zeros(ceil(lag) + 1, 1); % size of lag filter
b(end - 1) = ceil(lag) - lag ;
b(end)     = 1 - b(end - 1) ; % make sum = 1

% calculate initial condition using the slope of the start of x and
%   lag the data in x using the filter function.

L_b = length(b); % 1 + number of previous inputs required for lag-filter
y = 0:L_b;
y = y';
p = polyfit(y, x(1: L_b + 1), 1); % linear polyfit
previous_inputs = polyval(p, -y(2:end));
z = filtic(b, 1, [], previous_inputs);

x_lag = filter(b, 1, x, z);

if flip_me
    x_lag = flipud(x_lag);
end
