%% make_gradT_odas
% Generate the gradient of temperature using the raw pre-emphasized
% temperature signal (default) of the first-difference of temperature
% (optional for compatibility with the legacy method).
%%
% <latex>\index{Functions!make\_gradT\_odas}</latex>
%
%%% Syntax
%  gradT = make_gradT_odas(scalar_vector_with_pre_emphasis, scalar_info)
%
% * [scalar_vector_with_pre_emphasis] A vector containing the raw
%           pre-emphasized temperature signal to be converted into
%           temperature gradient in the direction of profiling.
% * [scalar_info] a structure of information used to make the conversion.
% * []
% * [gradT] The gradient of temperature in the direction of profiling, in
%           units of degrees C (or Kelvin) per meter.
%
%%% Description
% This function converts the raw pre-emphasized data
% ($\texttt{scalar\_vector\_with\_pre\_emphasis}$) into physical units of
% degrees C (or Kelvin) per meter. The information to do this conversion is
% in the fields within the structure named $\texttt{scalar\_info}$, which 
% must contain the following items:
%
% * [obj] The object containing the configuration file for the data. See
%       the function setupstr for more information.
% * [fs] The scalar of the sampling rate of the raw data.
% * [name_with_pre_emphasis] The string containing the name of the signal
%       containing the pre-emphasized tempertature. This name will be use
%       to find the required parameters in the object obj. It must exist.
% * [name_without_pre_emphasis] The string containing the name of the
%       signal of tempertature without pre-emphasis. This name is also used
%       to find paramters in the object obj. If all parameters can be found
%       using [name_with_pre_emphasis], then this string can be empty.
% * [speed] A scalar or vector of the speed of profiling. If it is a
%       vector, then it must be the same size as
%       scalar_vector_with_pre_emphasis. 
% * [method] The string containing the method to be used to compute the
%       gradient of temperature. The default is 'high_pass'. If it is
%       'first_difference' then this is the method that will be used. All
%       other strings will be converted to the default of 'high_pass'.

% * 2016-12-19 RGL, Original version.
% * 2016-12-22 RGL and WID, cleaned up the code and descriptions.
% * 2016-12-22 RGL, added catch for case of broken thermistor for which
%       Z>=1. This causes R_T/R_0 <=0, and a crash of the log-function. 
% * 2017-03-26 RGL, improved the handling of initial conditions for low-
%       and high-pass filtering, to derive the gradient of temperature.
% * 2017-07-26 RGL, Corrected syntax error for xmp_therm, strdouble ->
%       str2double.


function gradT = make_gradT_odas(scalar_vector_with_pre_emphasis, varargin)

%-----------------------------------------------------
% -- Default Values ----------------------------------
%-----------------------------------------------------
default_name_without_pre_emphasis  = ''; % no string for this name
default_method                     = 'high_pass'; 

p = inputParser;
p.KeepUnmatched = true;
p.CaseSensitive = true;

val_numeric     = @(x) isnumeric(x) && isscalar(x);
val_speed       = @(x) isnumeric(x) && isscalar(x)   || isempty(x);
val_vector      = @(x) (isnumeric(x) && isvector(x)) || isempty(x);
val_string      = @(x) ischar(x);
val_struct      = @(x) isstruct(x);

addRequired(  p, 'scalar_vector_with_pre_emphasis', val_vector);
addParamValue(p, 'speed',                           val_speed);
addParamValue(p, 'name_with_pre_emphasis',          val_string);
addParamValue(p, 'name_without_pre_emphasis', default_name_without_pre_emphasis, ...
                                                    val_string);
addParamValue(p, 'fs',                              val_numeric);
addParamValue(p, 'method', default_method,          val_string);
addParamValue(p, 'obj',                             val_struct)

%-----------------------------------------------------
% -- Parse Arguments ---------------------------------
%-----------------------------------------------------
parse(p, scalar_vector_with_pre_emphasis, varargin{:});

% extract the variables

T_dT                      = scalar_vector_with_pre_emphasis;
speed                     = p.Results.speed;
method                    = p.Results.method;
fs                        = p.Results.fs;
name_with_pre_emphasis    = p.Results.name_with_pre_emphasis;
name_without_pre_emphasis = p.Results.name_without_pre_emphasis;
obj                       = p.Results.obj;

%-----------------------------------------------------
% -- Check -------------------------------------------
%-----------------------------------------------------
if isempty (T_dT)
    error('The vector to be converted into temperature gradient is empty')
end
if isempty (name_with_pre_emphasis)
    error('The name of the temperature channel with pre-emphasis is empty')
end
if any(speed <= 0) 
    error('the speed is <= 0')
end
if ~(length(speed) == length(scalar_vector_with_pre_emphasis) || length(speed) == 1)
    error (['length of speed not equal to ' name_with_pre_emphasis ' or 1'])
end
if fs <= 0 
    error('the sampling rate is <= 0')
end

if ~(strcmpi(method,'high_pass') || strcmpi(method,'first_difference'))
    method = 'high_pass';
    warning(['The method = ' method ...
        ' is not excepted for computing the temperature gradient.' ...
        'Using high_pass'])
end

%-----------------------------------------------------
% -- Get Basic Parameters ----------------------------
%-----------------------------------------------------
diff_gain = setupstr( obj, name_with_pre_emphasis, 'diff_gain' );
if isempty(diff_gain)
    error(['There is no diff_gain parameter for the channel ' ...
        name_with_pre_emphasis])
else
    diff_gain = str2double(diff_gain);
end

% Make a cell array of the parameters that we need and the optional ones.
params_we_need = {...
    'a', 'b', 'adc_fs', 'adc_bits', 'g', 'e_b', 't_0', ...
    'beta', 'beta_1', 'beta_2', 'beta_3', 'zero', 'type', ...
    'coef0', 'coef1'};

beta = {[]}; % This is needed because beta is a Matlab built-in function.
type = {[]}; % This is needed because type is a Matlab built-in function.

% First check the section of the setup-file string for the channel with
% pre-emphasis.
name = name_with_pre_emphasis; % make it shorter
for index = 1 : length(params_we_need)
     junk = setupstr( obj, name, params_we_need(index));
     eval([params_we_need{index} ' = junk;'])
end

% Now go to the section without pre-emphsis to get the parameters that may
% still be missing.
if ~isempty (name_without_pre_emphasis)
    name = name_without_pre_emphasis; % make it shorter
    for index = 1 : length(params_we_need)
        if isempty (eval(params_we_need{index}))
            junk = setupstr( obj, name, params_we_need(index));
            eval([params_we_need{index} ' = junk;'])
        end
    end
end

if ~isempty(beta) && isempty(beta_1)
    beta_1 = beta;
end

if isempty (zero)
    zero = 0;
else
    zero = str2double(zero);
end
if isempty (beta_2)
    beta_2 = inf;
else
    beta_2 = str2double(beta_2);
end
if isempty (beta_3)
    beta_3 = inf;
else
    beta_3 = str2double(beta_3);
end

if strcmpi(type, 'xmp_therm')
    a        = 0; 
    b        = 1; 
    adc_fs   = 4.096;
    adc_bits = 16;
    g        = 1; 
    e_b      = 4.096;
    T_0      = str2double(coef0);
    beta_1   = str2double(coef1);
else
    a        = str2double(a);
    b        = str2double(b);
    adc_fs   = str2double(adc_fs);
    adc_bits = str2double(adc_bits);
    g        = str2double(g);
    e_b      = str2double(e_b);
    T_0      = str2double(t_0);
    beta_1   = str2double(beta_1);
end


%-----------------------------------------------------
% -- Deconvolve T_dT to remove pre-emphasis ----------
%-----------------------------------------------------
fc           = 1 / (2*pi*diff_gain);
% [b_lp, a_lp] = butter(1, fc / (fs/2));
% T            = filter(b_lp, a_lp, T_dT); % Now it is deconvolved- pre-emphasis is removed.

T = deconvolve (name_with_pre_emphasis, [], T_dT, fs, obj);

%-----------------------------------------------------
% -- High-pass T_dT to get pure time-derivative ------
%-----------------------------------------------------
[b_hp, a_hp] = butter(1, fc / (fs/2), 'high');

z_ic = filtic( b_hp, a_hp, 0, T_dT(1) ); 

gradT        = filter(b_hp, a_hp, T_dT, z_ic); % high passed
gradT        = gradT / diff_gain; % Now a derivative wrt time

% Now compute the absolute temperature
if strcmpi(type, 'therm') || strcmpi(type, 'xmp_therm')
    Z = ((T -a) / b) * (adc_fs / pow2(adc_bits)) * 2 / (g * e_b);
end
if strcmpi(type, 't_ms')
    Z = T * (adc_fs / pow2(adc_bits)) + zero;    
    Z = ((Z - a) / b) *2 / (g * e_b);
end

R = (1 - Z) ./ (1 + Z);
R(R<0.1) = 1; % in case of broken thermistor 
log_R = log(R);

T = 1/T_0 + log_R / beta_1;
if isfinite(beta_2)
    T = T + log_R.^2 / beta_2;
end
if isfinite(beta_3)
    T = T + log_R.^3 / beta_3;
end
T = 1./ T; % now in Kelvin

%-----------------------------------------------------
% -- Calculate the temperature gradient --------------
%-----------------------------------------------------
if strcmpi(method, 'high_pass')
    eta = (b / 2) * pow2(adc_bits) * g * e_b / adc_fs;
    scale_factor = 1;
    if isfinite(beta_2)
        scale_factor = 1 + 2 * (beta_1 / beta_2) * log_R;
    end
    scale_factor = scale_factor .* T.^2 .* (1 + R).^2 ./ (2 * eta * beta_1 * R);
    
    gradT = scale_factor .* gradT; % The time-derivative of temperature
end

if strcmpi(method, 'first_difference')
    gradT = fs*[diff(T) ; 0];
    gradT(end) = gradT(end-1);
end

gradT = gradT ./ speed; % That is it.
end

