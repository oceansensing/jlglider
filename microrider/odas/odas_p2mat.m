%% odas_p2mat
% Convert a RSI raw binary data file into a Matlab mat-file
%%
% <latex>\index{Functions!odas\_p2mat}</latex>
%
%%% Syntax
%   result = odas_p2mat( fname, ... )
%
% * [fname] Name of the raw binary RSI data file to process (extension optional).
% * [...] Optional parameters provided as a structure and-or list of
%       parameter name-value pairs. See below for a listing of optional
%       parameters.
% * []
% * [result] Return structure with different contents depending on the
%       context in which the function is called. When odas_p2mat is called
%       with no input arguments, the result structure contains a set of
%       default input arguments and no processing is performed. When
%       odas_p2mat has fname specified, the structure contains the results
%       traditionally saved within the resulting mat-file and no mat-file
%       is generated. The mat-file is only generated if result is not
%       explicitly used.
%
%%% Description
% Loading, converting, and processing of raw data files ($\texttt{.p}$) is
% performed by the $\texttt{odas\_p2mat}$ function. This function works
% with all RSI data files and automates much of the work involved in
% processing data files.
%
% This function reads and processes a raw binary $\texttt{p}$-file to
% generate data vectors converted into physical units. Data vectors are
% either saved within a $\texttt{mat}$-file named from the
% $\texttt{p}$-file or returned within the $\texttt{result}$ structure.
%
% Optional input parameters control the data-conversion process. Each
% parameter has a default value that is used if the parameter is not
% explictly specified. To simplify using these parameters, a structure
% containing all optional parameters, set to their default values, is
% returned when this function is called without input arguments. For
% example, 
%
%    >> default_parameters = odas_p2mat()
%
%    default_parameters =
%     MF_extra_points: 0
%                MF_k: 4
%            MF_k_mag: 1.7000
%              MF_len: 256
%           MF_st_dev: 'st_dev'
%        MF_threshold: []
%                 aoa: []
%      constant_speed: []
%       constant_temp: []
%        gradC_method: 'high_pass'
%        gradT_method: 'high_pass'
%          hotel_file: []
%        speed_cutout: 0.0500
%           speed_tau: []
%         time_offset: 0
%             vehicle: ''
%
% You can then change the parameters (which are the fields within the
% returned structure $\texttt{default\_parameters}$) to values suited to
% your particular processing requirements. 
%
%%% Vehicle Specification
% An instrument is identified by the $\texttt{vehicle}$ parameter. 
%
% * [vehicle] String identifying the vehicle that carries the RSI
% instrument. Typically empty [ ].
%
% The vehicle should be identified in the $\texttt{[instrument\_info}]$
% section of the configuration-file. The function searches the
% configuration string and tries to identify the vehicle. If it cannot find
% the vehicle, it will read this parameter from the input arguments. If the
% input parameter is empty, it will assume that the vehicle is
% $\texttt{vmp}$. You must specify this parameter only if it is not in the
% configuration-string (and the data is not from a VMP), or if you 
% want to explicitly override the value in the configuration string. This is
% bad practice but useful for testing purposes. Ultimately, you should
% correct the configuration string using the functions
% $\texttt{extract\_setupstr}$ and $\texttt{patch\_setupstr}$.
%
% The recognized $\texttt{vehicle}$-types are listed in the
% $\texttt{default\_vehicle\_attributes.ini}$ file included in the ODAS
% Matlab Library. 
%
%
%%% Profiling Speed
% Parameters used to determine the speed of profiling speed are:
%
% * [constant_speed] Speed used to generate gradients and convert
%       shear probe data into physical units. Leave empty if a constant
%       speed is not desired. Must be empty or positive. Default = [ ].
%       Units are m/s.
% * [hotel_file] Name of the mat-file containing speed information from a
%       source outside of the raw RSI data file. Default = [ ]
% * [speed_tau] Time-scale [s] of smoothing filter applied to speed data.
%       Default is vehicle dependent. See the
%       default_vehicle_attributes.ini file. 
% * [speed_cutout] Minimum profiling speed, [m/s], used to convert data 
%       into physical units. Slower speeds are set to this value. Default =
%       0.05. 
%
% The $\texttt{constant\_speed}$ parameter takes precedence over all other
% methods of deriving the speed of profiling. Specifying a
% $\texttt{constant\_speed}$ is useful if you have unreliable speed data,
% for example, when profiling through large turbulent up- and down-drafts.
% Or, if you have a vehicle for which it is not possible to calculate the
% speed using the raw RSI data file, and you do not yet have a hotel-file. 
%
% If $\texttt{constant\_speed}$ is empty, the speed of profiling is
% calculated from data in the RSI data file, or optionally, from data in a
% hotel-file. The algorithm used to calculate speed is determined by the
% $\texttt{vehicle}$ parameter (see above). If the speed algorithm requires
% data vectors not available in the current data file, the required data
% vectors must be provided within a hotel-file. You can also use a
% hotel-file if you prefer an alternative method of estimating the speed.
% See the ``hotel''-scripts for more information.
%
% The purpose of $\texttt{speed\_cutout}$ is to avoid division by a small
% number. The raw shear probe signals and all gradients are generated by
% dividing them by the speed of profiling. The estimated speed is very
% small when an instrument is not actually profiling, such as a VMP while
% it is suspended near the surface before being released to fall freely.
% Shear probe and gradient data are meaningless, during such occasions, and
% could be extremely large if the speed is not set to a cut-out value.
%
%%% Median Filter
%
% * [MF_len] Length of data segment [samples] used to calculate the median
%       and standard deviation. Default = 256.
% * [MF_threshold] Data points are defined as bad when they differ by more
%       than this value from the median. Default = [ ].
% * [MF_st_dev] empty string, or 'std_dev'.
% * [MF_k] empty, or the scale factor for the MF_st_dev option. Default = 4.
% * [MF_k_mag] Similar to MF_k but applied to the magnetometer signal.
%       Default = 1.7
% * [MF_extra_points] Number of points, adjacent to bad points, to also
%      replace with interpolated values. Default = 0.
%
% A median filter is used remove flyers, and other erroneous data, from the
% ADV (acoustic doppler velocimeter) velocity data. See the function
% $\texttt{median\_filter}$ for a more detailed description of these parameters. The
% median filter can be disabled by setting either $\texttt{MF\_threshold}$
% or $\texttt{MF\_k}$ to $\texttt{NaN}$ or $\texttt{inf}$. 
%
%%% Other parameters
%
% * [aoa] Angle-of-attack for a glider [degrees]. It is the difference
%       between the glide-path angle and the pitch angle. Default = 3.
% * [time_offset] Time [s] added to the instrument time. The resulting
%       absolute time is offset from the time recorded within the data file
%       header. This facilitates time synchronization with other sources of
%       data and a correction of an erroneously configured instrument clock.
% * [constant_temp] Temperature used for calculating the kinematic 
%       viscosity of water when measured temperature is unreliable, or
%       unavailable. Default = [ ]. 
% * [gradT_method] A string containing the method used to calculate the
%       gradient of temperature. There are two methods -- 'high_pass' and
%       'first_difference'. The high_pass method applies a high_pass filter
%       to the pre-emphasized temperature signal to make it a
%       time-derivative for all frequencies. The first-difference method
%       uses the high-resolution temperature to estimate the gradient by
%       way of the first-difference operator. See RSI Technical Note 005
%       for more information. Default = 'high_pass'.
% * [gradC_method] Same as gradT_method but for micro-conductivity signals.
% 

% Version History
%
% * 2015-03-12 RGL Original version.
% - 2015-04-19 RGL, Added features to handle a hotel file for a Remus AUV.
%       Changed the designation "speed" file to "hotel" file for speed and 
%       other vehicular information. 
% - 2015-04-27 RGL, Added "filesep" to filenames so that we donot get
%       dot-file names.
% - 2015-07-03 WID, Modified algorithm for calling deconvolve and convert.
%       Explicit "lists" are no longer required.
% - 2015-07-05 WID, Complete rework.  Fixed numerous minor bugs and changed
%       how some algorithms worked.  Should now be more flexable.
% - 2015-10-30 RGL, Document corrections.
% - 2015-11-01 RGL Reverted back to orginal version of trim_preemphasis_name
% - 2015-11-12 WID Added type t_ms to list of scalar gradients
% - 2015-11-18 RGL Document corrections.
% - 2015-12-31 WID Oops, forgot about the xmp.
% - 2016-01-21 WID Fixed bug in trim_preemphasis_name() function.
% - 2016-06-05 WID Allow the file to be moved into a new folder without
%                  crashing when loading the MAT file for a second time.
% - 2016-08-24, RGL, corrected error with checking of parameters of median
%           filter for Vector velocity data.
% - 2016-10-07, RGL, The Vector velocity components [U V W] were not being
%           saved because we forgot to pre-pend a "d." to their names. 
% - 2016-12-19, RGL, Modified the method for calculating the gradient of
%           temperature. odas_p2mat now supports 'high_pass' and
%           'first_difference' methods. 
% - 2016-12-21, RGL, Modified the method for calculating the gradient of
%           conductivity. It is now similar to that of temperature
%           gradient.
% - 2017-01-12 WID Disabled the fast version of hotel file data vectors.
% - 2017-01-13 WID Re-enabled fast vector for speed from hotel file.
% - 2017-01-23 WID Re-enabled fast vector for P from hotel file.
% - 2017-01-25 WID Added support for JAC EMC source for speed.
% - 2017-03-23, RGL, spelling mistake in description -- fhigh_pass instead
%           of high_pass. 
% - 2018-03-19, JMM, Fixed a bug that would not regenerate .mat file if a
%           ql_info field changed to from an empty vector. 
% - 2018-03-19 WID Improved debug output - code simplification. Add type to
%           list of things to check for.
% - 2020-10-11, JMM. Convert single variables to double in extract_vector, 
%           so that interp1 function works. Also added a warning.
% ==============================================

function result = odas_p2mat( fname, varargin )

%
% Default values for optional fields
default_speed_cutout    = 0.05; % dBar/s or m/s
default_constant_speed  = [];
default_constant_temp   = [];
default_vehicle         = '';
def_value_vehicle       = 'vmp';
default_aoa             = [];
default_speed_tau       = [];
default_time_offset     = 0;
default_hotel_file      = [];
default_MF_len          = 256;
default_MF_threshold    = [];
default_MF_st_dev       = 'st_dev';
default_MF_k            = 4;
default_MF_k_mag        = 1.7;
default_MF_extra_points = 0;
default_gradT_method    = 'high_pass';
default_gradC_method    = 'high_pass';

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

addRequired(  p, 'fname',          val_string);
addParamValue(p, 'speed_cutout',   default_speed_cutout,   val_positive);
addParamValue(p, 'constant_speed', default_constant_speed, val_positive);
addParamValue(p, 'constant_temp',  default_constant_temp);
addParamValue(p, 'vehicle',        default_vehicle,        val_string);
addParamValue(p, 'aoa',            default_aoa,            val_numeric);
addParamValue(p, 'speed_tau',      default_speed_tau,      val_numeric);
addParamValue(p, 'time_offset',    default_time_offset,    val_numeric);
addParamValue(p, 'hotel_file',     default_hotel_file,     val_string);
addParamValue(p, 'MF_len',         default_MF_len,         val_positive);
addParamValue(p, 'MF_threshold',   default_MF_threshold,   val_positive);
addParamValue(p, 'MF_st_dev',      default_MF_st_dev,      val_string);
addParamValue(p, 'MF_k',           default_MF_k,           val_positive);
addParamValue(p, 'MF_k_mag',       default_MF_k_mag,       val_positive);
addParamValue(p, 'MF_extra_points',default_MF_extra_points,val_numeric);
addParamValue(p, 'gradT_method',   default_gradT_method,   val_string);
addParamValue(p, 'gradC_method',   default_gradC_method,   val_string);

%input_parameters = join_arguments(varargin); % save for record keeping

% Parse the arguments.
parse(p, fname, varargin{:});

input_parameters = p.Results;

% Index results from the "p" variable.  Makes access much easier.
p = p.Results;


% Perform last stages of input validation.
% First the median filter parameters
if p.MF_len < 8,
  error('Invalid parameter MF_len=%f.  Value must be positive and >= 8.', p.MF_len);
end
if p.MF_k < 2.5,
  error('Invalid parameter MF_k=%f.  Value must be positive and >= 2.5.', p.MF_k);
end
if ~isempty(p.MF_st_dev) && ~strcmpi(p.MF_st_dev, 'st_dev')
  error('MF_st_dev must be set to the string st_dev or be empty');
end

% end of input argument checking.
% 
%



% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%%%%%%%%   OPEN DATA FILE
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% Try to open the RSI binary data file. It must have the extention '.p' or
% '.P'. Also, if a mat-file of this name already exists, no conversion is
% done.
errormsg = sprintf('Unable to find input file: %s', p.fname);
[P,N,E,fname_full] = file_with_ext(p.fname, {'.p','.P',''}, errormsg);

% The expected mat-file name.  The original fname variable can include
% wildcards so we have to decompose the returned full name of the file when
% searching for the .mat file.
[mP,mN] = fileparts(fname_full);

% First check if a mat-file alrady exists. If so, there is no conversion or
% any other action.
[mmP,mmN,mmE,matname_full] = file_with_ext([mP filesep mN], {'.mat','.MAT'});
if ~isempty(matname_full)
    disp(['Returning data from existing MAT-file: ' matname_full]);
    result = load(matname_full);
    % If the MAT and P files were moved, the "fullPath" variable will be
    % out of date.  Must update with current path.
    if isfield(result, 'fullPath')
        result.fullPath = fname_full;
    end
    status = verify_mat_p_differ(fname_full, result, input_parameters);
    user_input = '';
    if ~isempty(status)
        while isempty(user_input) || ...
                ~(strcmpi('y',user_input(1)) || ...
                  strcmpi('n',user_input(1)))
            warning(['MAT file: "' matname_full '" exists but can not ' ...
                     'be loaded for reason: "' status '".']);
            user_input = input('Delete MAT file and continue? [y,n] ', 's');              
        end
        if strcmpi('y',user_input(1))
            delete(matname_full);
        else
            return
        end
    else
        return
    end
end
matname_full = [mP filesep mN '.mat'];

% So, we next check this file for bad buffers
bad_records = check_bad_buffers(fname_full);
if ~isempty(bad_records)
    disp(['The file ' fname_full ' has bad buffers.'])
    disp('Making backup copy ...')
    
    % Backup current data file - patching bad buffers will modify the file.
    % Never overwrite the original data file.
    original_name = [P filesep N '_original' E];
    if ~exist(original_name, 'file')
        copyfile(fname_full, original_name);
    end
    
    disp('Calling function patch_odas ...')
    [bad_records, fix_manually] = patch_odas(fname_full);
    if ~isempty(fix_manually)
        disp('Some records cannot be patched,')
        disp('Calling function fix_bad_buffers ...')
        fix_bad_buffers(fname_full);
    end
end

% Read the binary file but do not convert it into a mat-file.  Unlike
% previously, data is transfered via the "data" structure.  This is faster.
% The mat-file will be generated at the end of the function.
[vars, d] = read_odas(fname_full); 
if isempty(vars), error(['Could not read binary file ' fname_full]); end

% -------------------------------------------------------------------------
% Generate a configration object - for quicker references.
% -------------------------------------------------------------------------
if ~isfield(d, 'cfgobj')
    obj = setupstr(d.setupfilestr);
else
    obj = d.cfgobj;
end


% Find default vehicle parameters for this vehicle.
p.vehicle = strip(p.vehicle);
if isempty(p.vehicle)
    p.vehicle = strip(char(setupstr(obj, 'Instrument_info', 'vehicle')));
    if isempty(p.vehicle), p.vehicle = def_value_vehicle; end
end
% All vehicles are defined in lower case characters.
p.vehicle = lower(p.vehicle);

%%%%%% Default vehicle parameters.
%
% Load and save the default vehilce attributes for the specified vehicle.
vehicle_info = iniparser( 'default_vehicle_attributes.ini' );
if ~isfield(vehicle_info, p.vehicle)
    error('Specified vehicle: "%s" not a recognized vehicle.', p.vehicle);
end
d.vehicle_info = vehicle_info.(p.vehicle);



% -------------------------------------------------------------------------
% Load the hotel file
% -------------------------------------------------------------------------
%
% t_fast and t_slow are required so loading the hotel file must occur after
% the data file is loaded using read_odas.

% Read the record length from the configuration string - almost always 1.
recsize = str2double(setupstr(obj, 'root', 'recsize'));

% Calculate the start time in Matlab time format.
t_start = d.filetime;
t_start = addtodate(t_start, -recsize, 'second');
t_start = addtodate(t_start, p.time_offset, 'second');

% Open hotel file if specified and save into structure dd.
if ~isempty(p.hotel_file)
    [hP,hB,hE,hF] = file_with_ext(p.hotel_file, ...
                                 {'', '.mat', '.MAT'}, ...
                                 'Specified hotel file not found.' );
    hotelFile = load(hF);
    
    for name = fieldnames(hotelFile)'
        field = char(name);
        [fast, slow] = extract_vector(hotelFile, field, t_start, d.t_fast, d.t_slow);
        if ~isempty(fast)
            disp(['      vector: "' field '" from hotel file']);
            % Only extract the fast vector for channels "speed" and "P"
            if strcmpi(field, 'speed') || strcmpi(field, 'P')
                d.([field '_fast']) = fast;
            end
            d.([field '_slow']) = slow;
        end
    end
end


% -------------------------------------------------------------------------
% Generate year-day time vectors
% -------------------------------------------------------------------------
y = datevec(t_start);
t_start_YD = t_start - datenum(y(1),1,0);
d.t_fast_YD = t_start_YD + d.t_fast/(3600*24);
d.t_slow_YD = t_start_YD + d.t_slow/(3600*24);

% Values adjusted by time_offset - differ from header.
d.Year  = y(1);
d.Month = y(2);
d.Day   = y(3);
d.Hour  = y(4);
d.Minute= y(5);
d.Second= floor(y(6));
d.Milli = y(6) - d.Second;


% -------------------------------------------------------------------------
% Extract the address matrix from configuration string
% -------------------------------------------------------------------------

% Read the address matrix from the configuration string.
n_fast = str2double(setupstr(obj, 'root',   'no-fast' ));
n_slow = str2double(setupstr(obj, 'root',   'no-slow' ));
n_rows = str2double(setupstr(obj, 'matrix', 'num_rows'));
rows = setupstr(obj, 'matrix', 'row[0-9]*');
matrix = zeros(n_rows,n_fast+n_slow);
for i=1:length(rows)
    row = textscan(rows{i}, '%f');
    matrix(i,:) = row{1}';
end


% -------------------------------------------------------------------------
% Generate channel list from configuration string
% -------------------------------------------------------------------------

% Valid channels have the three properties - id, name, and type.
ch_types = setupstr(obj, '', 'type', '');
ch_names = setupstr(obj, '', 'name', '');
ch_ids   = setupstr(obj, '', 'id',   '');

% Test each section to see if it has the required properties.  Valid
% sections are added to the "channels" list.
channels = {};
for chtype = ch_types
    if ~find(strcmp(chtype, ch_names(:)), 1), continue; end
    if ~find(strcmp(chtype, ch_ids(:)), 1), continue; end
    % Channels must also be in the address matrix.
    keep = true;
    for id = str2num(char( setupstr(obj, char(chtype), 'id') ))
        if isempty(find(id == matrix, 1))
            warning(['Channel "' char(setupstr(obj, char(chtype), 'name')) ...
                     '" declared but not found within the address matrix.']);
            keep = false;
            break;
        end
    end
    % Note we use the setupstr function here to ensure the correct case for
    % the chtype variable.  "chtype" is always forced to lowercase.
    if keep,
        channels(end+1) = setupstr(obj, char(chtype), 'name'); %#ok<*AGROW>
    end
end



% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%%%%%%%%   DECONVOLUTION
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% Test if the channel should be deconvolved into a "hres" channel.  This
% involves looking for two channels, one named "X" and the other "X_dX".
% When found, remove the channels from the list of channels, perform the
% deconvolution, then add the channel "X_hres" to the list.  Note that
% channels that have "diff_gain" parameters, are not of type "shear", and
% do not follow the previous naming convention, are also applied to the
% deconvolve algorithm.

convert_list = {channels{:}; channels{:}};
for i = 1:size(channels,2)  
    % Only channels with diff_gain parameters are deconvolved.
    if isempty(setupstr(obj, channels{i}, 'diff_gain')), continue; end
    
    % Shear probes have a diff_gain but are not deconvolved.
    if ~isempty(setupstr(obj, channels{i}, 'type', 'shear|xmp_shear'))
        continue
    end
    
    % Look for a non-preemphasized version of the channel.  If not found,
    % default to empty ([]).  If found, the resulting channel should be
    % saved with the "_hres" extension added to the base channel name.
    non_preemph = [];
    j = i;
    [tok, mat] = regexp(channels{i}, '(\w+)_d\1', 'tokens', 'match');
    if ~isempty(tok)
        % Token matches a "X_dX" channel name - check to see of a channel
        % named "X" exists.  If not, set the "X" channel to the same index
        % as "X_dX" -- ensures the algorithm works.
        j = find(strcmp(tok{1}{1}, channels(:)), 1);
        if ~isempty(j)
            % Generate an empty data array with the name that this channel
            % would have been given if the non-preemphasised channel
            % existed.  This allows the deconvolve algorithm to work.
            non_preemph = d.(channels{j});
        else
            j = i;
        end
                
        % Purge the two channels from the convert list.  Will later replace
        % with the "hres" channel while maintaining the correct source for
        % calibration coefficients (convert_list{2,:}).
        convert_list{1,i} = '';
        convert_list{1,j} = [tok{1}{1} '_hres'];
    end
    
    % Count number of times the channel occurs within an address
    % matrix.  Set the rate to fs_fast * occurances / row_count.
    id = str2double( setupstr(obj, channels{i}, 'id') );
    rate = d.fs_fast * length(find(id==matrix(:))) / n_rows;
    
    d.(convert_list{1,j}) = deconvolve(channels{i}, ...
                                        non_preemph, ...
                                        d.(channels{i}), ...
                                        rate, obj);
end
% end of deconvolution



% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%%%%%%%%   CONVERSION
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% Convert to physical units. The channels to convert are saved within
% the "convert_list{1,:}" cell array.  The source of the calibration
% coefficients are stored within "convert_list{2,:}".
for item = convert_list
    if isempty(item{1}), continue; end
    d.(item{1}) = convert_odas(d.(item{1}), item{2}, obj);
end


% Convert all "_hres" channels into "_slow" and "_fast" versions.  This
% facilitates a more consistant naming scheme thereby allowing hotel files
% to be of more use.
for field = fieldnames(d)'
    [tok, mat] = regexp(field{1}, '(\w+)_hres', 'tokens', 'match');
    if isempty(tok), continue; end
    
    name_fast = [tok{1}{1} '_fast'];
    name_slow = [tok{1}{1} '_slow'];
    
    % If a hotel file was used, the variable might already be set. Use
    % hotel file vector when available.
    if isfield(d, name_fast), continue; end
    
    % Generate _slow / _fast vectors.
    d.(name_slow) = interp1_if_req(d.(field{1}), d.t_slow);
    d.(name_fast) = interp1_if_req(d.(field{1}), d.t_fast);
    
    % Remove "_hres" channel to prevent saving within .mat file.
    d = rmfield(d, field{1});
end



% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%%%%%%%%   TEMPERATURE
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% Temperature_fast used to calculate the kinematic viscosity of water.
done = false;   % End search for temperature_fast.

% Constant temperature declared.  Override any hotel file values.
if ~isempty(p.constant_temp)
    if isnumeric(p.constant_temp)
        d.temperature_fast =  p.constant_temp * ones(size(d.t_fast));
        p.temperature_source = 'Declared constant temperature';
    else
        % Temperature declared - a string value representing the name of
        % the temperature vector that should be used. Do not extract the
        % value yet - do it below when finding an available thermistor.
        [tok, mat] = regexp(p.constant_temp, '(\w+)(_fast|_slow)+', 'tokens', 'match');
        if ~isempty(tok), p.constant_temp = tok{1}{1}; end
    end
    done = true;
end

% Temperature value provided within hotel file.
if ~done && isfield(d, 'temperature_fast')
    p.temperature_source = 'Found within hotel file.';
    done = true;
end

% Temperature from an existing vector should be used.  If a name was
% declared within "p.constant_temp", use it.  Otherwise use the first
% available channel or type "therm".
if ~done || ~isnumeric(p.constant_temp)
    
    if ~isnumeric(p.constant_temp)
        % The vector name was declared within p.constant_temp
        therms = {p.constant_temp};
    else
        % Use the first available thermistor - alpha-numerically sorted.
        therms = setupstr(obj, '', 'type', 'therm|t_ms|xmp_therm');
    end
    
    % If no thermistor is found, do not set "p.temperature_source".
    for i = 1:length(therms)
        % Many thermometers could be available.  Sort based on name then
        % select the first available thermometer.  If an alternative
        % thermometer should be used, extract the converted channel data
        % and genetate 
        therms = sortrows(therms');
        name = setupstr(obj, therms{i}, 'name');
        name = trim_preemphasis_name(name);
        if ~isfield(d, [char(name) '_fast']), continue; end
        d.temperature_fast   = d.([char(name) '_fast']);
        p.temperature_source = [char(name) '_fast'];
        done = true;
    end
end

% Temperature not defined.  Use a default constant value.
if ~done
    d.temperature_fast   = 10 * ones(size(d.t_fast));
    p.temperature_source = 'Forced to constant temperature of 10C';        
end



% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%%%%%%%%   FALL RATE
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% Fall rate is typically used to calculate speed for verticle profilers.
% Only generate a fall rate vector if a pressure vector exists and the fall
% rate was not previously declared within an hotel file.

if ~isfield(d, 'W_slow') && isfield(d,'P_slow')
    
    % Smooth the resulting vector based on default fall rate for the
    % current vehicle unless a specific smooting facter is explicitly
    % defined.
    tau = d.vehicle_info.tau;
    if ~isempty(p.speed_tau), tau = p.speed_tau; end

    d.W_slow = gradient(d.P_slow, 1/d.fs_slow);    
    [b,a] = butter(1, (0.68/tau)/(d.fs_slow/2));
    d.W_slow = filtfilt(b, a, d.W_slow);

    % It might be faster to interploate the data but probably does not make
    % much of a difference.
    d.W_fast = gradient(d.P_fast, 1/d.fs_fast);
    [b,a] = butter(1, (0.68/tau)/(d.fs_fast/2));
    d.W_fast = filtfilt(b, a, d.W_fast);
end



% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%%%%%%%%   SPEED
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% Speed is required for converting shear probe signals into physical units
% and for calculating gradients.  Several different algorithms can be used
% to calculate speed.  The correct algorithm is determined by the vehile
% type.  Alternatively, a constant speed can be declared or a speed vector
% can be provided within an hotel file.

% Commonly used values simplify reading the algorithm.
done = false;                          % End search for speed.
type = d.vehicle_info.speed_algorithm; % Algorithm to use.
aoa  = d.vehicle_info.aoa;             % Angle of attack.
tau  = d.vehicle_info.tau;             % Smoothing factor

% Default vehicle values assigned if required.
if isempty(p.aoa),       p.aoa = aoa;       end    
if isempty(p.speed_tau), p.speed_tau = tau; end


% A declared constant speed takes precidence over any other speed values.
if ~isempty(p.constant_speed) % force a constant speed
    d.speed_fast =  p.constant_speed * ones(size(d.t_fast));
    d.speed_slow =  p.constant_speed * ones(size(d.t_slow));
    p.speed_source = 'Forced to constant speed by input parameter';
    done = true;
end 

% Speed exists!!  Must have been from an hotel file.  Hotel files take
% precidence over everything except a declared constant speed.
if ~done && isfield(d, 'speed_slow')
    p.speed_source = 'Hotel file';
    done = true;
end

% We do not have speed in an hotel file - report error.
if ~done && strcmpi(type, 'hotel')
    error('Missing speed.  Declare a constant speed or provide speed within an hotel file.');
end

% We must calculate speed based on the fall rate for a verticle profiler.
if ~done && strcmpi(type, 'pressure')
    d.speed_fast = abs(d.W_fast);
    d.speed_slow = abs(d.W_slow);
    p.speed_source = 'Rate of change of pressure';
    done = true;
end

% We must calculate speed based on the fall rate and glide angle.
if ~done && strcmpi(type, 'glide')
    if ~isfield(d,'W_slow')
        error(['Can not calculate glide speed without either P ' ...
               '(pressure) or W (fall rate).  Provide values within ' ...
               'an hotel file.']);
    end
    
    % Find vectors - use hotel file vectors if available.
    if isfield(d,'Incl_Y'),      Incl_Y = d.Incl_Y;      end
    if isfield(d,'Incl_Y_slow'), Incl_Y = d.Incl_Y_slow; end
    if isfield(d,'Ax'),      Ax = d.Ax;      end
    if isfield(d,'Ax_slow'), Ax = d.Ax_slow; end
    if isfield(d,'Ay'),      Ay = d.Ay;      end
    if isfield(d,'Ay_slow'), Ay = d.Ay_slow; end
    if isfield(d,'Az'),      Az = d.Az;      end
    if isfield(d,'Az_slow'), Az = d.Az_slow; end
    
    if exist('Incl_Y','var') % We have a new MR with inclinometers
        glide_angle = abs(Incl_Y) + p.aoa;
    elseif exist('Ax','var') && exist('Ay','var') && exist('Az','var')
        glide_angle = abs(asind(Ax/9.81)) + p.aoa;
        glide_angle = interp1_if_req(glide_angle, d.t_slow);
    else
        error(['Unable to determin inclination when calculating speed. '...
               'Provide either "Incl_Y" or "Ax,Ay,Az" vectors within an'...
               ' hotel file.']);
    end
    d.speed_slow = abs(d.W_slow) ./ sind(glide_angle);
    p.speed_source = 'Rate of change of pressure and glide angle';
    done = true;
end

% Speed determined from the three velocity vectors - possibley recored from
% a "vector" instrument. 
if ~done && strcmpi(type, 'vector')

    d.Vector_bad_points = [];

    % Use velocity values from hotel file in place of recorded values when
    % available.  Hotel file values have "_slow" appended to their name.
    if isfield(d, 'U'), U = d.U; end
    if isfield(d, 'U_slow'), U = d.U_slow; end
    
    if isfield(d, 'V'), V = d.V; end
    if isfield(d, 'V_slow'), V = d.V_slow; end

    if isfield(d, 'W'), W = d.W; end
    if isfield(d, 'W_slow'), W = d.W_slow; end
    
    if ~(exist('U','var') && exist('V','var') && exist('W','var'))
        error(['Velocity data (U,V,W) is missing.  Please generate an ' ...
               'hotel file with the required values.']);
    end
    
    % Use median filter to remove flyers.  Filter parameters must be
    % correctly defined.
    if isempty(p.MF_threshold) && isempty(p.MF_st_dev)
        error('You must specify MF_threshold when MF_st_dev is empty.');
    end
    
    if true
        % isfinite(p.MF_threshold) && isfinite(p.MF_k)
        [results, d.Vector_bad_points] = ...
            median_filter([U V W], p.MF_threshold, p.MF_len, ...
                          p.MF_extra_points, p.MF_st_dev, p.MF_k);
        d.U = results(:,1);
        d.V = results(:,2);
        d.W = results(:,3);
    end
    
    d.speed_slow = sqrt(U.^2 + V.^2 + W.^2);
    d.speed_source = 'Recorded U V W';
    done = true;
end

% Speed determined from a AEM1-G electromagnetic velocity sensor on a
% channel of type 'aem1g_a' or 'aem1g_d'. Using vehicle (such as auv_emc)
% with speed_algorithm = emc.
if ~done && strcmpi(type, 'emc')
    
    % Find channels using the assumed type.  Use the first channel we find
    % and ignore the rest.
    name = '';
    for t = {'aem1g_a', 'aem1g_d', 'jac_emc'}
        for ch = setupstr(obj, '', 'type', t{1})
            % Extract the name from the configuration string - returns 
            % correct case.
            chname = char(setupstr(obj, ch{1}, 'name'));
            if isfield(d, chname)
                name = chname;
                break;
            end
        end
        if ~isempty(name), break; end
    end
    
    if isempty(name)
        error('Unable to find channel of type "aem1g_a" or "aem1g_d".');
    end

    if isfield(d, [name '_slow']) && isfield(d, [name '_fast'])
        % Use data amended using a hotel file if available. Such data will
        % contain a "_fast" and "_slow" postfix to the name.
        d.speed_fast = d.([name '_fast']);
        d.speed_slow = d.([name '_slow']);
    else
        % Otherwise, use original data when hotel file data not provided.
        d.speed_fast = interp1_if_req(d.(name), d.t_fast);
        d.speed_slow = interp1_if_req(d.(name), d.t_slow);
    end
    
    d.speed_source = ['EMC channel: ' name];
    done = true;
end

% Last one - vehicle assumes a constant speed.  For example, a stand within
% a river.
if ~done && strcmpi(type, 'constant')
    speed = d.vehicle_info.constantSpeed;
    d.speed_fast = ones(size(d.t_fast)) .* speed;
    d.speed_slow = ones(size(d.t_slow)) .* speed;
    d.speed_source = 'Use of constant speed due to vehicle type.';
end


% Perform final smoothing and generate a fast vector for speed - if a slow
% vector was already generated.
if isfield(d,'speed_slow')
    % Generate the final vector if required.
    if ~isfield(d,'speed_fast')
        d.speed_fast = interp1_if_req(d.speed_slow, d.t_fast);
    end

    % Vehicle specific smoothing of the speed vector.
    [b,a] = butter(1, (0.68/p.speed_tau)/(d.fs_slow/2));
    d.speed_slow = filtfilt(b,a,d.speed_slow);
    [b,a] = butter(1, (0.68/p.speed_tau)/(d.fs_fast/2));
    d.speed_fast = filtfilt(b,a,d.speed_fast);

    % Fish out the speeds that are smaller then the cutout of speed_cutout
    % to avoid singularities.
    d.speed_slow(d.speed_slow < abs(p.speed_cutout)) = abs(p.speed_cutout);
    d.speed_fast(d.speed_fast < abs(p.speed_cutout)) = abs(p.speed_cutout);
end



% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%%%%%%%%   SHEAR CONVERSIONS
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% Make final conversion of shear probes in to physical units.  All channels
% of type "shear" will be converted.
for ch = setupstr(obj, '', 'type', 'shear|xmp_shear')
    name = setupstr(obj, char(ch), 'name');
    if find(strcmp(name{1}, channels))
        d.(name{1}) = d.(name{1}) ./ d.speed_fast.^2;
    end
end



% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%%%%%%%%   SCALAR GRADIENTS
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% Calculate gradients for temperature and micro-conductivity data vectors.
% Channels of type therm and ucond that have diff_gain parameters are
% selected for conversion.
%
% This section has been modified to use look only for micro-conductivity signals.
% for t = {'ucond'}
%     for ch = setupstr(obj, '', 'type', t{1})
%         if isempty(setupstr(obj, ch{1}, 'diff_gain')), continue; end
%         
%         name = setupstr(obj, ch{1}, 'name');
%         name = char(name);
% 
%         [tok, mat] = regexp(name, '(\w+)_d\1', 'tokens', 'match');
%         if ~isempty(tok), name = tok{1}{1}; end
%         
%         % Gradient source and destination vector names.
%         src  = [name '_fast'];
%         dest = ['grad' name];
% 
%         % Generate gradient signal/source names
%         if isfield(d, src)
%             d.(dest) = [ 0; diff(d.(src)) * d.fs_fast ] ./ d.speed_fast;
%             d.(dest)(1) = d.(dest)(2);
%         end
%     end
% end
% 
% This section for micro-conductivity signals
% We need the names of the channels with and without pre-emphasis in order
% to get all of the required coefficients.
for t = {'ucond'}
    for ch = setupstr(obj, '', 'type', t{1})

        % Skip channels that do not have a "diff_gain" parameter.
        if isempty(setupstr(obj, ch{1}, 'diff_gain')), continue; end

        % Extract the name from the configuration string - returns correct case.
        name = char(setupstr(obj, ch{1}, 'name'));
        name_with_pre_emphasis = name;
        
        % Find channel without pre-emphasis - if it exists.
        [tok, mat] = regexp(name, '(\w+)_d\1', 'tokens', 'match');
        if ~isempty(tok) %&& isfield(d, tok{1}{1})
            name_without_pre_emphasis = tok{1}{1};
        else
            name_without_pre_emphasis = [];
        end
        
        % Assign vector name for the result.
        dest = ['grad' name_with_pre_emphasis];
        if ~isempty(name_without_pre_emphasis)
            dest = ['grad' name_without_pre_emphasis];
        end
        
        % Assign input parameters for the gradient function.
        scalar_vector_with_pre_emphasis       = d.(name_with_pre_emphasis);
        scalar_info.name_without_pre_emphasis = name_without_pre_emphasis;
        scalar_info.name_with_pre_emphasis    = name_with_pre_emphasis;
        scalar_info.fs                        = d.fs_fast;
        scalar_info.speed                     = d.speed_fast;
        scalar_info.obj                       = obj;
        scalar_info.method                    = p.gradC_method;

        % Make the temperature gradient
        d.(dest) = make_gradC_odas(...
            scalar_vector_with_pre_emphasis, scalar_info);
    end
end


% This section for thermistor signals
% We need the names of the channels with and without pre-emphasis in order
% to get all of the required coefficients.
for t = {'therm','t_ms','xmp_therm'}
    for ch = setupstr(obj, '', 'type', t{1})

        % Skip channels that do not have a "diff_gain" parameter.
        if isempty(setupstr(obj, ch{1}, 'diff_gain')), continue; end

        % Extract the name from the configuration string - returns correct case.
        name = char(setupstr(obj, ch{1}, 'name'));
        name_with_pre_emphasis = name;
        
        % Find channel without pre-emphasis - if it exists.
        [tok, mat] = regexp(name, '(\w+)_d\1', 'tokens', 'match');
        if ~isempty(tok) %&& isfield(d, tok{1}{1})
            name_without_pre_emphasis = tok{1}{1};
        else
            name_without_pre_emphasis = [];
        end
        
        % Assign vector name for the result.
        dest = ['grad' name_with_pre_emphasis];
        if ~isempty(name_without_pre_emphasis)
            dest = ['grad' name_without_pre_emphasis];
        end
        
        % Assign input parameters for the gradient function.
        scalar_vector_with_pre_emphasis       = d.(name_with_pre_emphasis);
        scalar_info.name_without_pre_emphasis = name_without_pre_emphasis;
        scalar_info.name_with_pre_emphasis    = name_with_pre_emphasis;
        scalar_info.fs                        = d.fs_fast;
        scalar_info.speed                     = d.speed_fast;
        scalar_info.obj                       = obj;
        scalar_info.method                    = p.gradT_method;

        % Make the temperature gradient
        d.(dest) = make_gradT_odas(...
            scalar_vector_with_pre_emphasis, scalar_info);
    end
end



% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%%%%%%%%   MAGNETOMETER
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% Cleanup the flyers in the magnetometer data. MF_K_mag = 1.7 seems to work
% well for the magnetometer.
for ch = setupstr(obj, '', 'type', 'magn')
    name = setupstr(obj, char(ch), 'name');
    badname = [name{1} '_bad'];
    [d.(name{1}), d.(badname)] = ...
        median_filter(d.(name{1}), [], 8*round(d.fs_slow), 0, 'st_dev', p.MF_k_mag);
end



% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%%%%%%%%   SAVE / RETURN DATA
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% Save the original input parameters and p.  Many of the values will be the
% same but "p" will contain any modified values and debug / status strings
% that were generated during processing.
d.input_parameters = input_parameters;
d.params = p;

if nargout == 0
    disp(['Saving raw data into file: ' matname_full]);
    save(matname_full, '-struct', 'd'); % save all variables.
end

% return the variables in a structure
result = d;

end

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%%%%%%%%   THE END   %%%%%%%%.



% Assemble input arguments into a 1X1 stracture.  Contents of the input
% struture are joined with parameter name / value pairs.  The resulting
% structure can then be saved for future reference.
function args = join_arguments(input)
    if isempty(input)
        args = [];
        return
    end

    skip = 0;
    for k = 1:length(input)
        if skip, skip = 0; continue; end
    
        if isstruct(input{k})
            args = input{k};
            continue
        else
            if ~ischar(input{k})
                error('Invalid input argument');
            end
            args.(input{k}) = input{k+1};
            skip = 1;
        end
    end
end



% Interpolate the vector only if required.  This function assumes the time
% interval should stay the same and only the number of points should
% change.
function vector = interp1_if_req(vec, t_new)

    ratio = length(vec) / length(t_new);
    
    % Do nothing.
    if ratio == 1, vector = vec; return; end
    
    % Must modify the number of data points.  Generate an appropriate
    % source time vector from which we can interpolate.  This actually
    % works for anything.
    fs = t_new(2) - t_new(1);
    t = 0:fs/ratio:t_new(end)+fs-fs/ratio;
    
    % The resulting value of "t" should be the same length as "vec".
    vector = interp1(t, vec, t_new, 'spline', 'extrap');
end



function [fast, slow] = extract_vector(hotelFile, field, t_start, t_fast, t_slow)
% Extract fast and slow data vectors from an hotel file.  Resulting vectors
% are alligned to t_fast and t_slow time domains.  Extra data is clipped
% and missing data is interpolated.  Missing data at the vector extents is
% not interpolated but instead is held constant until the data vectors
% matches both t_fast and t_slow.
%
% The hotel data vectors can be in a few different time formats.  Both
% MATLAB time and Year-Day time are acceptable.  Values < 1000 are assumed
% to be in Year-Day while values > 1000 are assumed to be in MATLAB time.
% Time values starting at 0 are assumed to be relative to the start of data
% acquisition.

fast = [];
slow = [];

if isfield(hotelFile, field) && isstruct(hotelFile.(field)) && ...
        isfield(hotelFile.(field), 'data') && ...
        isfield(hotelFile.(field), 'time')
    
    % Data and time is available.  Convert time if required then ensure
    % date exists over the required absolute time interval. 
    absolute_fast = t_fast ./ (3600*24) + t_start;
    absolute_slow = t_slow ./ (3600*24) + t_start;
    
    % Get time and data from the hotel file structure.
    time = hotelFile.(field).time;
    data = hotelFile.(field).data;
    
    % Check that fields are doubles
    if isa(time,'single')
        warning('Hotel file: Time vector is a single. Converting to double.')
        time = double(time);
    end
    if isa(data,'single')
        warning('Hotel file: Data vector is a single. Converting to double.')
        data = double(data);
    end
    
    % Remove NAN and INF data / time values.
    bad_idx = union(union(find(isnan(data)), find(isnan(time))), ...
                    union(find(isinf(data)), find(isinf(time))));
    idx = setxor(1:length(time), bad_idx);
    time = time(idx);
    data = data(idx);
    
    % Allocate the required output arrays
    fast = zeros(length(t_fast),1);
    slow = zeros(length(t_slow),1);
    
    % Accept data/time vectors with only a single value - another form of
    % constant.  Add an additional value to allow interpolation to work.
    if length(time) == 1
        time(2) = time(1) + 1;
        data(2) = data(1);
    end
    
    % Convert the hotel file time to Matlab datenum if required.  For time
    % Year-Day, we do not have the hotel file year so we assume it is the
    % same as ours.
    if time(1) == 0
        % Time vector is in "Instrument" time - in units of seconds and
        % starting at the start of the profile.
        time = time ./ (3600*24) + t_start;
    elseif time(1) >= 1 && time(1) < 1000
        % Time Year-Day format.  Starting at the start of the calendar year
        % and in units of days.
        y = datevec(t_start);
        time = datenum(y(1), 1, time);
    elseif time(1) >= 1000 && time(1) < 1e8
        % Must be MATLAB time (in days) - do nothing.
    elseif time(1) >= 1e8
        % Time in seconds since 1970 - POSIX or UNIX time.
        time = datenum(datetime(time, 'ConvertFrom', 'posixtime'));
    else
        error('Invalid time standard.');
    end
    
    % The require absolute time domain required.  The hotel file likely has
    % a different domain.  If the hotel file is smaller we adjust the
    % domain over which we interpolate then manually patch in the
    % non-overlapping data.
    rsi_start = absolute_fast(1);
    rsi_last  = absolute_fast(end);
    
    hf_start = time(1);
    hf_last  = time(end);
    
    start = max(rsi_start, hf_start);
    last  = min(rsi_last,  hf_last );
    
    if start >= last
        warning(['Error loading data from hotel file field: ', field, ...
            '.  Data does not overlap the required time duration.']);
        return
    end
    
    % Start / stop indices for interpolation of hotel file data.
    hf_start_idx = max(1, find(time>=start,1,'first'));
    hf_last_idx  = min(length(time), find(time<=last,1,'last'));

    % Start / stop indices for interpolation of our data.
    % Adjust start and last to compensate for NANs.
    rsi_start_idx = max(1, ...
                        find(time(hf_start_idx)<=absolute_fast,1,'first'));
    rsi_last_idx  = min(length(absolute_fast), ...
                        find(time(hf_last_idx)>=absolute_fast, 1, 'last'));

    % Interpolate the fast channel then patch in any missing data.
    fast(rsi_start_idx:rsi_last_idx) = ...
                interp1( time(hf_start_idx:hf_last_idx), ...
                         data(hf_start_idx:hf_last_idx), ...
                         absolute_fast(rsi_start_idx:rsi_last_idx), ...
                         'spline');
    fast(1:rsi_start_idx) = fast(rsi_start_idx);
    fast(rsi_last_idx:end) = fast(rsi_last_idx);
       
    % Run the algorithm again, this time for a slow channel.
    rsi_start = absolute_slow(1);
    rsi_last  = absolute_slow(end);
    
    start = max(rsi_start, hf_start);
    last  = min(rsi_last,  hf_last );
    
    % Start / stop indices for interpolation of hotel file data.
    hf_start_idx = max(1, find(time>=start,1,'first'));
    hf_last_idx  = min(length(time), find(time<=last,1,'last'));
    
    % Start / stop indices for interpolation of our data.
    % Adjust start and last to compensate for NANs.
    rsi_start_idx = max(1, ...
                        find(time(hf_start_idx)<=absolute_slow,1,'first'));
    rsi_last_idx  = min(length(absolute_slow), ...
                        find(time(hf_last_idx)>=absolute_slow, 1, 'last'));

    % Interpolate the slow channel then patch in any missing data.
    slow(rsi_start_idx:rsi_last_idx) = ...
                interp1( time(hf_start_idx:hf_last_idx), ...
                         data(hf_start_idx:hf_last_idx), ...
                         absolute_slow(rsi_start_idx:rsi_last_idx), 'pchip');
    slow(1:rsi_start_idx) = slow(rsi_start_idx);
    slow(rsi_last_idx:end) = slow(rsi_last_idx);

end

end



function r = iniparser( filename )

[fid,msg] = fopen(filename, 'r');
if ~isempty(msg), r = msg; return; end

r = struct();

sec = 'root';

while true
    line = fgetl(fid);
    
    % The end of file
    if isnumeric(line), break; end
    
    % Skip comments
    if isempty(line) || line(1) == ';', continue; end
    
    % Assign a new current section
    [tok, m] = regexp(line, '^\s*\[\s*(.+?)\s*\]\s*$', 'tokens', 'match');
    if ~isempty(tok)
        sec = tok{1}{1};
        continue;
    end
    
    % Parameter name / value entry
    [tok, m] = regexp(line, '^\s*(.+?)\s*=\s*(.+?)\s*$', 'tokens', 'match');
    if ~isempty(tok)
        value = str2double(tok{1}{2});
        if isnan(value), value = tok{1}{2}; end
        r.(sec).(tok{1}{1}) = value;
        continue;
    end
    
    % Skip everyting else
end

fclose(fid);

end



% Verify that a MAT file corresponds with a P file when using the specified
% parameters.  Used to determin if the MAT file is valid for the specific P
% file or if it has to be deleted.
function result = verify_mat_p_differ(fname_full, mat, params)
    result = '';
    if ~isfield(mat, 'odas_version') || ...
            mat.odas_version ~= odas_version_info()
        result = 'ODAS Library Version Differs';
        return
    end
    
    fid = fopen_odas(fname_full, 'r');
    header = fread(fid, 64, 'uint16');
    cfgstr = fread(fid, header(12), '*char');
    fclose(fid);
    if ~strcmpi(mat.setupfilestr, cfgstr')
        result = 'Configuration Strings Differ';
        return
    end
    
    % Check each of the input parameters. Compare the fields from params
    % that are also found within mat.input_parameters. Matching field names
    % should then have their values compared.  One last catch - the fname
    % parameter should be excluded.
    for p = fieldnames(params)'
        % Check if input parameter should be skipped
        if ~isfield(mat.input_parameters, char(p)) || strcmp('fname', char(p))
            continue
        end
        
        value1 = params.(char(p));
        value2 = mat.input_parameters.(char(p));
            
        % Skip when both parameters are empty
        if isempty(value1) && isempty(value2)
            continue
        end
            
        % Error due to missing one parameter
        if isempty(value1) || isempty(value2)
            if isempty(value1)
                result = ['Input parameters differ: e.g. ' char(p) ': new = [], old = ' num2str(value2)];
            elseif isempty(value2)
                result = ['Input parameters differ: e.g. ' char(p) ': new = ' num2str(value1) ', old = []'];
            end
            
            return
        end
        
        % Both values are present, check if they differ.
        % Three checks, check type, numeric, and string.        
        if ~strcmp(class(value1), class(value2))
           result = sprintf('Input Parameters Differ in Type: e.g. %s: new = %s, old = %s', ...
                            char(p), class(value1), class(value2));
           return
        end
            
        if isnumeric(value1) && value1 ~= value2
           result = sprintf('Input Parameters Differ: e.g. %s: new = %s, old = %s', ...
                            char(p), num2str(value1), num2str(value2));
           return
        end

        if ischar(value1) && ~strcmpi(value1, value2)
           result = sprintf('Input Parameters Differ: e.g. %s: new = "%s", old = "%s"', ...
                            char(p), value1, value2);
           return
        end
    end
end

% Strings with the name 'X_dX' are trimmed to 'X'.
function name = trim_preemphasis_name( name )
%    [tok, mat] = regexp(name, '(\w+)_d\1', 'tokens', 'match');
%    if ~isempty(tok), name = tok{1}{1}; end
%%%%% Alternative version that does not use regular expressions.
    if iscell(name), name = char(name); end
    idx = strfind( name, '_d' );
    if isempty(idx), return; end
    name = char(name);
    prefix = name(1:idx(1)-1);
    postfix = name(idx(1)+2:end);
    if strcmp(prefix, postfix),
        name = prefix;
    end
end



