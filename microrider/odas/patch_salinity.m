%% patch_salinity
% Clean, lag and match-filter conductivity data from a JAC-CT sensor
%%
% <latex>\index{Functions!patch\_salinity}</latex>
%
%%% Syntax
%   str = patch_salinity( filename, ... )
%
% * [filename] Name of the mat-file containing JAC-CT data.
% * [...] Structure or list of key/value pairs containing input parameters
%       used to control the patching procedure. 
% * []
% * [result] Structure containing pressure, temperature, conductivity,
%       corrected conductivity, salinity, fall-rate, and time data for use
%       with further processing or visualization. When patch_salinity is 
%       called with no input parameters, default parameter values are 
%       returned within a structure.
%
%%% Description
% Read a mat-file containing JAC conductivity and temperature
% data and use it to generate a clean and corrected conductivity signal
% suitable for computing salinity.  The mat-file is then updated with the
% corrected conductivity and salinity values.
%
% A structure containing default parameters values is returned when 
% patch_salinity is called with no input arguments. This structure can be 
% examined, edited, and used as an input to the patch_salinity function.
% Structure fields that are not recognized are silently ignored.
%
%%% Optional Inputs
%
% * [lag] The lag of the conductivity signal relative to the temperature
%     signal, in units of samples. Default = 1.
% * [f_CT] The half-power response frequency [Hz] of the thermometer
%     relative to the conductivity cell. Default = 0.6.
% * [P_min] Only use data with a pressure greater than P_min, in units of
%     dbar. Default = 1. 
% * [speed_min] Only use data that has a magnitude of the rate-of-change of
%     pressure greater than speed_min, in units of dbar per second. Default
%     = 0.50.
% * [MF_len] Length, in samples, of the segments used by the median filter
%     to remove spikes in the conductivity signal due to plankton and other
%     detritus. See the function median_filter for more detailed
%     information on parameter values. Default = 64. 
% * [MF_threshold] The threshold for detection of conductivity spikes.
%     Default = [ ].
% * [MF_st_dev] Flag that indicates the use of the standard deviation to set
%     the threhold. Default = 'st_dev'.
% * [MF_k] The scale factor to apply to the standard deviation to determine
%     the threshold. Default = 4.
% * [MF_extra_points] The number of extra points around a spike to also
%     replace with a linear interpolation. Default = 0.
%
% The function goes through the folowing steps. First, it uses the median
% filter to remove conductivity spikes. It then lags the conductivity
% signal by $\texttt{lag}$ samples. Next, it low-pass filters the
% conductivity signal using a first-order Butterworth filter with a cut-off
% frrequency of $\texttt{f\_CT}$ $\si{\hertz}$ to match this signal to the
% temperature signal. It then computes the salinity using the processed
% conductivity and the original temperature and pressure. Finally, it
% appends to the data file the low-pass filtered conductivity signal,
% $\texttt{JAC\_C\_LP}$, the salinity, $\texttt{JAC\_S}$, and the indices,
% $\texttt{JAC\_C\_bad\_points}$, to the data removed by the median filter.

% Version History
%
% * 2015-03-31 (RGL) Original version.
% * 2015-04-13 (WID) Configured documentation for publishing.
% * 2015-10-30 (RGL) Changed documentation
% * 2015-11-18 (RGL) Corrected the documentation.

function result = patch_salinity( filename, varargin )

%
% Default values for optional fields
default_lag              = 1; % lag on C in units of samples
default_f_CT             = 0.60; % respose of C relarive to T, in Hz
default_P_min            = 35; % min depth in dbar
default_speed_min        = 0.50; % in dbar/s
default_MF_len           = 64; % Length of median filter
default_MF_threshold     = []; % Threshold
default_MF_st_dev        = 'st_dev';% Filter based on standard deviation
default_MF_k             = 4;% Threshold is MF_k times std
default_MF_extra_points  = 0;% extra points to remove around a spike

if ~nargin,
    for d = whos('default_*')',
        param = regexp(d.name, 'default_(.*)', 'tokens');
        result.(char(param{1})) = eval(d.name);
    end
    return
end

p = inputParser;
p.KeepUnmatched = true;
p.CaseSensitive = true;

val_numeric  = @(x) (isnumeric(x) && isscalar(x))        || isempty(x);
val_positive = @(x) (isnumeric(x) && isscalar(x) && x>0) || isempty(x);
val_string   = @(x) ischar(x)                            || isempty(x);

addRequired(  p, 'file_name',                              val_string);
addParamValue(p, 'lag',            default_lag,            val_numeric);
addParamValue(p, 'f_CT',           default_f_CT,           val_positive);
addParamValue(p, 'P_min',          default_P_min,          val_positive);
addParamValue(p, 'speed_min',      default_speed_min,      val_positive);
addParamValue(p, 'MF_len',         default_MF_len,         val_positive);
addParamValue(p, 'MF_threshold',   default_MF_threshold,   val_positive);
addParamValue(p, 'MF_st_dev',      default_MF_st_dev,      val_string);
addParamValue(p, 'MF_k',           default_MF_k,           val_positive);
addParamValue(p, 'MF_extra_points',default_MF_extra_points,val_numeric);

% Parse the arguments.
parse(p, filename, varargin{:});

names = fieldnames(p.Results);
for name = names'
    eval([char(name) ' = p.Results.' char(name) ';']);
end

% Now save these parameter values so that they can be placed into the
CT_info = p.Results;

% Perform last stages of input validation.
% First the median filter parameters
if ~isempty(MF_len) && MF_len < 8,
  MF_len = default_MF_len;
  warning(['You must specify MF_len >= 8. Using MF_len = ' num2str(MF_len)]);
end
if ~isempty(MF_k) && MF_k < 2.5,
  MF_k = default_MF_k;
  warning(['You must specify MF_k >= 2.5. Using MF_k = ' num2str(MF_k)]);
end
if ~isempty(MF_st_dev) && ~strcmpi(MF_st_dev,'st_dev')
  error('MF_st_dev must be set to the string st_dev or be empty');
end
result = [];
% Now check the speed information file, if any, and extract the structure dd.
% end of input argument checking.
% 
%%
%____________________________________________________________
% File opening etc.
% Check file name, check for a .mat file, and open the .p file if no *.mat
% file
% _________________________________________________________________

% Try to open the RSI binary data file. It must have the extention '.p' or '.P'.
% Also, if a mat-file of this name already exists, no conversion is done.

[file_name_path, file_name_base, file_name_extension]=...
    file_with_ext(filename, {'','.mat'} );

% Check if a mat-file exists. If so, there carry on.
if ~exist([file_name_path file_name_base '.mat'],'file') % Then the mat-file does not exists
	error(['Cannot find the mat-file ' [file_name_path file_name_base]])
end

load(filename)

m = find (P_slow >= P_min & W_slow >= speed_min);

% Apply median filter to get rid of conductivity spikes from plankton
[C_clean, JAC_C_bad_points] = median_filter(...
    JAC_C(m), MF_threshold, MF_len, MF_extra_points, MF_st_dev, MF_k );

[b,a] = butter(1, f_CT/(fs_slow/2)); % coefficients for matching filter
JAC_C_LP = JAC_C; % pre-assign space
JAC_C_LP(m) = C_clean; % fill in the despiked data section

JAC_C_LP(1:lag)     = JAC_C(1)*ones(lag,1); 
JAC_C_LP(1+lag:end) = JAC_C_LP(1 : end - lag); % apply lag
JAC_C_LP = filter(b, a, JAC_C_LP);

JAC_S = salinity(P_slow, JAC_T, JAC_C_LP); % lag and filter corrected salinity


result.JAC_S    = JAC_S;
result.JAC_C_LP = JAC_C_LP;
result.JAC_T    = JAC_T;
result.JAC_C    = JAC_C;
result.P_slow   = P_slow;
result.t_slow   = t_slow;
result.JAC_C_bad_points   = m(JAC_C_bad_points);

save (filename, 'JAC_C_LP', 'JAC_S', 'JAC_C_bad_points','-append')



