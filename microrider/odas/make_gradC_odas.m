%% make_gradC_odas
% Generate the gradient of conductivity using the raw pre-emphasized
% conductivity signal (default), or the first-difference of conductivity
% (optional for comatibility with the legacy method).
%%
% <latex>\index{Functions!make\_gradC\_odas}</latex>
%
%%% Syntax
%  gradC = make_gradC_odas(scalar_vector_with_pre_emphasis, scalar_info)
%
% * [scalar_vector_with_pre_emphasis] A vector containing the raw
%           pre-emphasized conductivity signal to be converted into
%           conductivity gradient in the direction of profiling.
% * [scalar_info] a structure of information used to make the conversion.       
% * []
% * [gradC] The gradient of conductivity in the direction of profiling, in
%            units of milli-Siemens per centimeter per meter. 
%
%%% Description
% This function converts the raw pre-emphasized data
% ($\texttt{scalar\_vector\_with\_pre\_emphasis}$) into physical units of
% milli-Siemens per centimeter per meter. The information to do this
% conversion is in the fields within the structure named 
% $\texttt{scalar\_info}$, which must contain the following items:
%
% * [obj] The object containing the configuration file for the data. See
%       the function setupstr for more information.
% * [fs] The scalar of the sampling rate of the raw data.
% * [name_with_pre_emphasis] The string containing the name of the signal
%       containing the pre-emphasized conductivity. This name will be use
%       to find the required parameters in the object obj. It must exist.
% * [name_without_pre_emphasis] The string containing the name of the signal
%       of conductivity without pre-emphasis. This name is also used to find
%       the parameters in the object obj. If all parameters can be found
%       using [name_with_pre_emphasis], then this string can be empty.
% * [speed] A scalar or vector of the speed of profiling. If it is a
%       vector, then it must be the same size as
%       scalar_vector_with_pre_emphasis. 
% * [method] The string containing the method to be used to compute the
%       gradient of conductivity. The default is 'high_pass'. If it is
%       'first_difference', then this is the method that will be used. All
%       other strings will be converted to the default of 'high_pass'.

% * 2016-12-22 RGL, Original version.
% * 2017-03-26 RGL, improved the handling of initial conditions for low-
%       and high-pass filtering, to derive the gradient of temperature.


function gradC = make_gradC_odas(scalar_vector_with_pre_emphasis, varargin)

%-----------------------------------------------------
% -- Default values ----------------------------------
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

C_dC                      = scalar_vector_with_pre_emphasis;
speed                     = p.Results.speed;
method                    = p.Results.method;
fs                        = p.Results.fs;
name_with_pre_emphasis    = p.Results.name_with_pre_emphasis;
name_without_pre_emphasis = p.Results.name_without_pre_emphasis;
obj                       = p.Results.obj;

%-----------------------------------------------------
% -- Some final checking -----------------------------
%-----------------------------------------------------
if isempty (C_dC)
    error('The vector to be converted into conductivity gradient is empty')
end
if isempty (name_with_pre_emphasis)
    error('The name of the conductivity channel with pre-emphasis is empty')
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
        ' is not excepted for computing the conductivity gradient.' ...
        'Using high_pass'])
end

%-----------------------------------------------------
% -- Get Basic Parameters that are needed ------------
%-----------------------------------------------------
diff_gain = setupstr( obj, name_with_pre_emphasis, 'diff_gain' );
if isempty(diff_gain)
    error(['There is no diff_gain parameter for the channel ' ...
        name_with_pre_emphasis])
else
    diff_gain = str2double(diff_gain);
end

% Make a cell array of the parameters that we need and the optional ones.
params_we_need = {'a', 'b', 'adc_fs', 'adc_bits', 'K'};


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


a        = str2double(a);
b        = str2double(b);
adc_fs   = str2double(adc_fs);
adc_bits = str2double(adc_bits);
K        = str2double(K);

%-----------------------------------------------------
% -- Deconvolve C_dC to remove pre-emphasis ----------
%-----------------------------------------------------
%% 
fc           = 1 / (2*pi*diff_gain);
% [b_lp, a_lp] = butter(1, fc / (fs/2));
% C            = filter(b_lp, a_lp, C_dC); % Now it is deconvolved- pre-emphasis is removed.

C = deconvolve (name_with_pre_emphasis, [], C_dC, fs, obj);

%-----------------------------------------------------
% -- Deconvolve C_dC to get pure time-derivative -----
%-----------------------------------------------------
[b_hp, a_hp] = butter(1, fc / (fs/2), 'high');

z_ic = filtic( b_hp, a_hp, 0, C_dC(1) ); 

gradC        = filter(b_hp, a_hp, C_dC, z_ic); % high passed
gradC        = gradC / diff_gain; % Now a derivative wrt time

%-----------------------------------------------------
% -- Calculate conductivity gradient -----------------
%-----------------------------------------------------
if strcmpi(method, 'high_pass')
    gradC = adc_fs * gradC / (pow2(adc_bits) * b * K);
end

if strcmpi(method, 'first_difference')
    % First compute the conductivity
    C = ((adc_fs * C / pow2(adc_bits)) - a) / (b*K); % Siemens/m
    gradC = fs*[diff(C) ; 0];
    gradC(end) = gradC(end-1);
end

gradC = 10*gradC ./ speed; % That is it. Factor of 10 to get units of mS/cm

end

