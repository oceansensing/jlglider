%% cal_FP07_in_situ
% Determine calibration coefficients for a thermistor probe using in situ data.
%%
% <latex>\index{Functions!cal\_FP07\_in\_situ}</latex>
%
%%% Syntax
%   [T_0, beta, Lag] = cal_FP07_in_situ( file_name, 
%                           T_ref_string, T_string, SN, cal_info,... )
%
% * [file_name] Name of mat-file containing data used to
%       calibrate a thermistor. Note: Data needs to have been previously 
%       extracted from .P file using either quick_look or odas_p2_mat.
% * [T_ref_string] Name of vector within the mat-file (e.g. 'JAC_T') 
%       that contains the reference temperature in degrees celsius. Usually 
%       from a SBE4F thermometer or a JAC CT. 
% * [T_string] Name of thermistor to calibrate, typically 'T1' or 'T2'.
% * [SN] Serial number of thermistor.
% * [cal_info] Structure (optional) containing configuration parameters. 
%       A template is generated and returned when cal_FP07_in_situ is called 
%       with no input parameters. The parameters are described below.
% * [...] Optional configuration parameters to supplement, or override, 
%       those values included within cal_info. Inputs are accepted as
%       string/value pairs.
% * []
% * [T_0] Value of parameter T_0, used in the Steinhart-Hart equation. When
%       called with no input parameters, a structure containing default
%       input parameter values is returned for reference.
% * [beta] beta coefficients, in ascending order, of the fit to the Steinhart-Hart
%       equation. i.e. beta_1, beta_2, beta_3. 
% * [Lag] Delay, in seconds, between the thermistor and the reference
%       thermometer. Typically a negative value because the reference
%       sensor is usually behind the thermistor being calibrated.
%
%%% Description
%
% This function can be used to calibrate a FP07 thermistor probe using 
% in-situ data. The reference temperature data is usually provided by a Sea-Bird 
% SBE3F or a JAC-CT. The reference temperature data must be contained within 
% the specified mat-file. If necessary the data can be incorportated into the 
% file using a hotel file at the time of data conversion. This is typically 
% necessary when analyzing data collected with a glider. 
%
% This function processes the data as follows:
% (1) Selects the portion of the data file that will be used for the 
% calibration based on the input parameters. The range will be plotted if 
% plot_range = true. 
% (2) It then detrends the thermistor data and scales it so that it is 
% approximately aligned with the reference thermometer. A plot of the
% scaled signals will only be shown if plot_scaled = true. 
% (3) It then computes the cross-correlation coefficient between the thermistor 
% and the reference thermometer and estimates the lag between these two 
% signals. A plot will only be generated if plot_xcorr = true.
% (4) Next it does a regression based on the Steinhart-Hart equation and estimates 
% the coefficients based on the chosen order of the regression. A plot of the 
% natural logarithm of the resistance ratio against the inverse of the 
% absolute temperature is shown if plot_regress = true.
% (5) The regression coefficients are then used to convert the thermistor
% data into physical units. The depth profile of the calibrated thermistor
% data and the reference signal, and their difference, will be plotted if
% plot_result = true. 
%
% If $\texttt{cal\_FP07\_in\_situ}$ is called without input arguments, it 
% returns the default parameters used by $\texttt{cal\_FP07\_in\_situ}$. 
% You can then customize this structure to your particular processing 
% requirements. For example,
%
%    >> cal_info = cal_FP07_in_situ
%    
%    cal_info = 
% 
%             make_figures: 1
%                    order: 2
%               plot_range: 1
%             plot_regress: 1
%              plot_result: 1
%              plot_scaled: 0
%               plot_xcorr: 0
%            profile_min_P: 1
%            profile_min_W: 0.2000
%     profile_min_duration: 20
%              profile_num: 1
%             vehicle_info: []
%
% The configuration parameters (fields) within the structure $\texttt{cal\_info}$ 
% that control the behaviour of $\texttt{cal\_FP07\_in\_situ}$,are listed below. 
% They are grouped for clarity, but are all part of the single structure 
% $\texttt{cal\_info}$.
%
%%% Parameters that control the calibration
%
% * [order] Fit order to the Steinhart-Hart equation. Value can be 1,
%       2, or 3. Default = 2. For small temperature ranges, order = 1 is 
%       recommended.
%
%%% Parameters that specify a profile
% * [profile_num] Index to the requested profile from the set of detected
%       profiles.  The first profile, 1, is the default value.
% * [profile_min_P] The minimum pressure of a profile, [dbar]. Default = 1.
% * [profile_min_W] The minimum vertical speed of profiling  
%       in [dbar/s]. Default = 0.2. Use a smaller value, e.g. 0.1 dBar/s, 
%       for gliders.
% * [profile_min_duration] The minimum duration in which the minimum
%      pressure and speed must be satisfied [s]. Default = 20.
% * [vehicle_info] A structure found in all .mat data files that is used
%       to determine the direction of profiling. Default = [].
%       If empty, the information will be loaded from the datafile. If a 
%       change is desired, the profiling direction needs to be specified, 
%       i.e. Set vehicle_info.profile_dir to 'down' for vmps, 'up' for rvmps
%       or 'glide' for gliders.
%
%
%%%% Parameters that toggle data visualization
%
% * [make_figures] The parameter that determines if figures are generated.
%     make_figures = false suppresses the generation of figures to speed up
%     the data processing. Default = true.
% * [plot_range] A logical parameter that determines if a figure of the 
%     selected range is generated. Default = true.
% * [plot_scaled] A logical parameter that determines if a figure of the 
%     detrended signals and lag are plotted. Default = false.
% * [plot_xcorr] A logical parameter that determines if a figure of the 
%      cross correlation is plotted. Default = false.
% * [plot_regress] A logical parameter that determines if the regression
%      data and fit are plotted. Default = true.
% * [plot_result] A logical parameter that determines if figures of the 
%       calibration results (i.e. comparison and difference) are shown. 
%       Default = true.


%___________________________
%
% Version History
%
% * 2013-12-05 (RGL) original version.
% * 2013-12-06 (RGL) added varargin to define a profile with default values.
% * 2015-04-10 (WID) revised documentation for publishing.
% * 2015-04-27 (RGL) modified to allow the specification of the fit order.
% * 2015-07-27 (WID) use texstr in place of fix_underscore.
% * 2015-07-29 (WID) return default values when called with no input
%                       parameters.
% * 2015-10-27 (RGL) Changed description section.
% * 2016-06-07 (RGL) Changed legend call, added clf to the start of figures.
% * 2016-11-10 (RGL) Changed the call to the deconvolve function so that it
%      includes the thermistor signal without pre-emphasis. odas_p2mat uses
%      the the thermistor without pre-emphasis (if it is available) for
%      conversion in to physical units.
%      Consequently, the thermistor signal without pre-emphasis MUST be
%      called with this in situ calibration function. Otherwise, there is a
%      substantial error due to the offset (order ~10 counts) that may be
%      present in the signal with pre-emphasis. This version also includes
%      a test for the existence of the signal without pre-emphasis so that
%      it does not bomb in its absence. I also beautified the legends and
%      labels.
% * 2016-11-10 (RGL) Added a low-pass filter to the thermistor data, in the
%      case of a JAC-T reference, in order to get a tighter regression. The
%      cut-off frequency is sped dependent and follows the recommendation
%      for calculating salinity. There is an improvement when the
%      thermistor is filtered.
% * 2016-11-15 (RGL) Added fc in case of Sea-Bird thermometer.
% * 2017-11-28 (RGL) Added ability to handle both type=therm abd type=t_ms.
% * 2017-11-30 (RGL) Some more display changes and warning will be
%      suppressed in case the channel without pre-emphasis does not exist.
%      Made the naming of parameters (fields) consistent with usage in
%      quick_look.
% * 2019-06-27 (JMM) Changes to function plotting and inputs. Regression
% algorithm is unaffected. Changes include:
%       - Added an initial subplot of the pressure record, highlighting the 
%           portion of the file being used
%       - Load in vehicle_info section from .mat file (datafile is required)
%       - Change default profile_min_W to 0.2m/s to match quick_look inputs
%       - Added plotting flags to be able to supress figures.
%       - Changed variable name from T2 to T_prof to be more descriptive.
%       - Added a warning message if temperature range is too small for high
%         order fit. 
% * 2019-06-03 (JMM) Updated comments and documentation. 


function [T_0,beta,Lag] = cal_FP07_in_situ(file_name,T_ref_string,T_string,SN,varargin)

%-----------------------------------------------------------------
% ----- Default parameters ---------------------------------------
%-----------------------------------------------------------------
default_vehicle_info         = []; % will trigger profile_dir = 'down'
default_profile_min_P        = 1; % in dBar
default_profile_min_W        = 0.2; % in dbar/s
default_profile_num          = 1; % process the first profile
default_profile_min_duration = 20; % minimum duration [s] for a segment to be considered a profile
default_order                = 2; % The order of the fit to the SS equation
default_make_figures         = true; % render figures for data visulization
default_plot_range           = true; % flag to plot the range used for profiles
default_plot_scaled          = false; % flag to plot detrended signals showing lag
default_plot_xcorr           = false; % flag to plot cross-correlation
default_plot_regress         = true; % flag to plot regression
default_plot_result          = true; % flag to plot result of calibration


if ~nargin
    for d = whos('default_*')'
        param = regexp(d.name, 'default_(.*)', 'tokens');
        result.(char(param{1})) = eval(d.name);
    end
    T_0 = result;
    return
end

%-----------------------------------------------------------------
% ----- Parse Inputs ---------------------------------------------
%-----------------------------------------------------------------
p = inputParser;
p.KeepUnmatched = true;
p.CaseSensitive = true;

val_numeric     = @(x) isnumeric(x) && isscalar(x);
val_string      = @(x) ischar(x);
val_logical     = @(x) islogical(x);


addRequired(  p, 'file_name',           val_string);
addRequired(  p, 'T_ref_string',        val_string);
addRequired(  p, 'T_string',            val_string);
addRequired(  p, 'SN',                  val_string);
addParamValue(p, 'profile_num',          default_profile_num,          val_numeric);
addParamValue(p, 'profile_min_P',        default_profile_min_P,        val_numeric);
addParamValue(p, 'profile_min_W',        default_profile_min_W,        val_numeric);
addParamValue(p, 'profile_min_duration', default_profile_min_duration, val_numeric);
addParamValue(p, 'order',                default_order,                val_numeric);
addParamValue(p, 'vehicle_info',         default_vehicle_info);

addParamValue(p, 'make_figures', default_make_figures, val_logical);
addParamValue(p, 'plot_range',   default_plot_range,   val_logical);
addParamValue(p, 'plot_scaled',  default_plot_scaled,  val_logical);
addParamValue(p, 'plot_xcorr',   default_plot_xcorr,   val_logical);
addParamValue(p, 'plot_regress', default_plot_regress, val_logical);
addParamValue(p, 'plot_result',  default_plot_result,  val_logical);

% -- Parse the arguments.
parse(p, file_name, T_ref_string, T_string, SN, varargin{:});

% -- Define simple variables to use throughout function
profile_min_P        = p.Results.profile_min_P;
profile_min_W        = p.Results.profile_min_W;
profile_num          = p.Results.profile_num;
profile_min_duration = p.Results.profile_min_duration;
order                = p.Results.order;
vehicle_info         = p.Results.vehicle_info;
make_figures         = p.Results.make_figures;
plot_range           = p.Results.plot_range;
plot_scaled          = p.Results.plot_scaled;
plot_xcorr           = p.Results.plot_xcorr;
plot_regress         = p.Results.plot_regress;
plot_result          = p.Results.plot_result;
% end of input argument checking.

fig_num = 0;

% -----------------------------------------------------------------------
% --- Load the data -----------------------------------------------------
% -----------------------------------------------------------------------

%%%%
% We need the name of the FP07 signals with and without pre-emphasis. The
% name without pre-emphasis is usually the section in the setup.cfg-file
% that contains the processing parameters. However, sometimes this signal
% does not exist and the information has to be gathered from the section
% for the signal with pre-emphasis. For eample, with a Sea-glider.

T_with_pre_emphasis_string = [T_string '_d' T_string]; % Should be T1_dT1 or T2_dT2
T_without_pre_emphasis_string = T_string; % Should be T1 or T2

warning off; % In case some variables do not exist
load(file_name, ...
    'setupfilestr', ...
    T_without_pre_emphasis_string, ...
    T_with_pre_emphasis_string, ...
    T_ref_string, ...
    't_fast', 't_slow', 'P_slow', 'W_slow', 'fs_fast', 'fs_slow')
warning on;

if isempty(vehicle_info) % Load vehicle_info from file unless specified
    load(file_name,'vehicle_info') 
else
    warning(['Overriding vehicle_info used to convert to physical units'])

end

% -----------------------------------------------------------------------
% --- Get relevant temperature data (and deconvolve if necessary) -------
% -----------------------------------------------------------------------

% We must identify the name of the section in the setup.cfg-file that
% contains the processing paramters.
section_name = T_without_pre_emphasis_string;

% In case that there is no thermistor signal without pre-emphasis
if exist(T_without_pre_emphasis_string,'var')
    section_name = T_without_pre_emphasis_string;
    eval(['T_without_pre_emphasis = ' T_without_pre_emphasis_string ';'])
else
    T_without_pre_emphasis_string = '[]';
    section_name = T_with_pre_emphasis_string;
    T_without_pre_emphasis = [];
end

% In case that there is no thermistor signal with pre-emphasis, we can only
% look at the signal without pre-emphasis and do not have to bother with
% deconvolution of the pre-emphasized signal.
if ~exist(T_with_pre_emphasis_string,'var')
    T_with_pre_emphasis_string = '[]';
end

% make sure that there is at least one signal with FP07 data.
if ...
        isempty(T_without_pre_emphasis_string) && ...
        isempty(T_with_pre_emphasis_string)
    error(['Could not find any signals with base name = ' T_string])
end

% Finally, make sure that the temperature reference signal exists.
if ~exist(T_ref_string,'var')
    error(['Could not find any temperature reference signals with name = ' T_ref_string])
end

% Assign temperature vectors
eval (['T_ref = ' T_ref_string ';']) % T_ref is the reference thermometer, usually SBT or JAC_T
eval(['T = ' T_without_pre_emphasis_string ';']) % T is the thermistor signal without pre-emphasis


% If we have a signal with pre-emphasis, then we will use it to form the signal T.
if ~isempty(T_with_pre_emphasis_string)
    eval(['T = ' T_with_pre_emphasis_string ';']) % T is the thermistor signal with pre-emphasis
    T = deconvolve(...
        T_with_pre_emphasis_string, ...
        T_without_pre_emphasis, ...
        T, ...
        fs_fast, setupfilestr);
    
    % sampling rate ratio (used to downsample)
    ratio = round(fs_fast / fs_slow);
    
    % down size to match T_ref
    T = reshape(T, ratio, []); 
    T = mean(T)';
end

% -----------------------------------------------------------------------
% --- Plot temperature and pressure data for entire datafile ------------
% -----------------------------------------------------------------------

% common plot title
title_string{1} = ['\rmFP07 \itin situ\rm calibration , SN-' SN];

if plot_range && make_figures
    fig_num = fig_num+1;
    figure(fig_num), clf
    
    % Pressure 
    ax(1) = subplot(211);
    plot(t_slow, P_slow,'linewidth',1);
    ylabel('P [dBar]')
    set(gca,'Ydir','rev')
    title(title_string)
    hold all
    
    % Temperature 
    ax(2) = subplot(212);
    [axyy,p1,p2] = plotyy(t_slow, T,t_slow,T_ref);grid on
    set(axyy(2),'ycolor',[0 0.6 0])
    set(p2,'color',[0 0.6 0])
    xlabel('\it t \rm [s]')
    ylabel(axyy(1),[T_string '[counts]'])
    ylabel(axyy(2),[texstr(T_ref_string) ' [\circC]'])
    hold all
    
    linkaxes([ax,axyy],'x')
end

% -----------------------------------------------------------------------
% --- Extract profile data  ------------------------------------
% -----------------------------------------------------------------------


% Figure out which section of the file to use for the in situ calibration.
if isempty(vehicle_info)
    profile_dir = 'down';
else
    profile_dir = vehicle_info.profile_dir;
end

% Get the profile indices
if exist('P_slow','var') && exist('W_slow','var')
    if strcmpi(profile_dir, 'up') || strcmpi(profile_dir, 'down')
        profile = get_profile(P_slow, W_slow, profile_min_P, ...
            profile_min_W, profile_dir, ...
            profile_min_duration, fs_slow);
    elseif strcmpi(profile_dir,'glide')
        profile_down = get_profile(P_slow, W_slow, profile_min_P, ...
            profile_min_W, 'down', ...
            profile_min_duration, fs_slow);
        profile_up   = get_profile(P_slow, W_slow, profile_min_P, ...
            profile_min_W, 'up',   ...
            profile_min_duration, fs_slow);
        % Sort columns in ascending order
        profile = sort([profile_down profile_up],2);
    end
end

% Give warning and quit if requested profile doesn't exist
profiles_in_this_file = size(profile,2);
if profile_num > profiles_in_this_file
    warning(['There are only ' num2str(profiles_in_this_file) ' profiles in this file'])
    
    %  add info to title of plot
    if plot_range && make_figures
        title(ax(1),{['\rm ',texstr(file_name)],...
            ['Requesting profile ',num2str(profile_num), ...
            ' but ' num2str(profiles_in_this_file) ' profiles detected in this file']})
        
        drawnow
    end
    
    T_0 = [];
    return
end

profile_start = profile(1,profile_num);
profile_end   = profile(2,profile_num);
m = (profile_start:profile_end)';

% Highlight selected range on previous plot
if plot_range && make_figures
    % pressure
    plot(ax(1),t_slow(m), P_slow(m),'r','linewidth',2); grid on
    legend(ax(1),...
        'P', ...
        ['P Profile'], 'location', 'northeast')
    title(ax(1),title_string)
    
    % temperature
    p3 = plot(ax(2),t_slow(m), T(m),'r','linewidth',2);grid on
    legend([p1,p3,p2], T_string,[T_string ' Profile'], texstr(T_ref_string),...
        'location', 'northeast')
end

%----------------------------------
% - Warning if range is too small - 
%----------------------------------
if max(T_ref(m))-min(T_ref(m))<=8 & order>1
    warning(['Temperature range is less than 8 degrees '...
        '-> Recommend using FIRST-ORDER calibration.',...
        ' Exit function (typically by CRTL+C) and re-run with additional inputs,'...
        ' - i.e. ''order'',1 - '...
        ' or modify cal_info structure.',...
        ' Then, in setup file, delete entire ''beta_2=__'' line and ',...
        ' update T_0 and beta_1 values before patching.']), pause
end

%--------------------------------------------------------------------
% -- Detrend and scale signals so T and T_ref have same range -------
%--------------------------------------------------------------------
if plot_scaled && make_figures
    fig_num = fig_num+1;
    figure(fig_num), clf
    junk_T      = detrend(T(m));
    junk_T_ref  = detrend(T_ref(m));
    range_T_ref = max(junk_T_ref) - min(junk_T_ref);
    range_T     = max(junk_T) - min(junk_T);
    junk_T      = junk_T * range_T_ref / range_T; % T should now span the same range as T_ref.
    
    title_string{2} = [...
        'Detrended ' texstr(T_ref_string) ' & scaled ' texstr(T_string)];
    
    plot(t_slow(m), [junk_T  junk_T_ref]);grid on
    title (title_string)
    xlabel('\it t \rm [s]')
    ylabel('[ ^{\circ}C ]')
    legend(texstr(T_string), texstr(T_ref_string))
end

%-------------------------------------------------------------------------
% -- Low-Pass filter thermistor data and compute cross correlation -------
%-------------------------------------------------------------------------

% Low-pass filtering the thermistor data to make it more compatible with the JAC-T
if strcmp(T_ref_string,'JAC_T')
    W_mean = abs(mean(W_slow(m)));
    fc = 0.73 * sqrt(W_mean / 0.62); % in Hz
    [b,a] = butter(1, fc / (fs_slow/2));
    T = filter(b, a, T);
else
    fc = fs_slow/3; % It is a Sea-Bird Thermometer
    [b,a] = butter(1, fc / (fs_slow/2));
    T = filter(b, a, T);
end


% Compute cross-correlation
max_lag = round(10*fs_slow); % estimate of the max lag required to find the actual lag.
[bb, aa] = butter(2,4/(fs_slow/2)); % 4 Hz smoother to suppress high-frequency noise
[correlation, lags] = xcorr(...
    filter(bb,aa,detrend(diff(T(m)))),...
    filter(bb,aa,detrend(diff(T_ref(m)))),max_lag,'coeff');
[max_corr, m_lag] = max(abs(correlation));
junk_m = m_lag; % needed for figure
m_lag = m_lag - max_lag - 1;
Lag    = m_lag / fs_slow; % in seconds and should be negative.

% Plot
if plot_xcorr && make_figures
    fig_num = fig_num+1;
    figure(fig_num), clf
    plot(lags/fs_slow, correlation, m_lag/fs_slow, correlation(junk_m), 'r*');grid on
    xlabel('Lag [ s ]')
    legend_string_1 = ['max X_{corr} = ' num2str(max_corr,2)];
    legend_string_2 = ['@ \tau = ' num2str(m_lag/fs_slow,2) ' s'];
    legend(legend_string_1, legend_string_2, 'location', 'northeast')
    title_string{2} = texstr([...
        'X-correlation of ' T_string ...
        ' and ' T_ref_string]);
    title(title_string)
end

%-----------------------------------------------------
% -- Do regression to get thermistor coefficients ----
%           (using Steinhart-Hart equation)
%-----------------------------------------------------

% Copy only the profile data
T_ref_prof = T_ref(m);
T_prof     = T(m);
P_prof     = P_slow(m);

% First align the T and T_ref signals using m_lag.
if m_lag >0, m_lag = 0; end % m_lag is expected to be negative because 
                            % reference sensor is physically 'behind' probes
P_prof     = P_prof(1:end+m_lag); 
T_prof     = T_prof(1:end+m_lag);
T_ref_prof = T_ref_prof(1-m_lag:end); % shift reference temperature
T_ref_regress = T_ref_prof + 273.15; % in kelvin
T_ref_regress = 1 ./ T_ref_regress;

% Now gather information about the electronics for this thermistor.
my_object = setupstr( setupfilestr );
therm_type =    (char(setupstr( my_object, section_name, 'type')));
E_B = str2double(char(setupstr( my_object, section_name, 'E_B')));
a   = str2double(char(setupstr( my_object, section_name, 'a'  )));
b   = str2double(char(setupstr( my_object, section_name, 'b'  )));
G   = str2double(char(setupstr( my_object, section_name, 'G'  )));
adc_fs   = str2double(char(setupstr( my_object, section_name, 'adc_fs'  )));
adc_bits = str2double(char(setupstr( my_object, section_name, 'adc_bits'  )));
try zero = str2double(char(setupstr( my_object, section_name, 'adc_zero'  )));catch, zero = 0; end

% Compute non-dimensional thermistor voltage
if strcmp(therm_type, 'therm')
    factor = (adc_fs / 2^adc_bits)*2 / (G*E_B);
    Z = factor*(T_prof - a)/b;
elseif strcmp(therm_type, 't_ms')
    Z = T_prof * (adc_fs/2^adc_bits) + zero;
    Z = ((Z - a)/b) *2 / (G*E_B);
end

% Compute resistance ratio for this thermistor.
RT_R0 = (1 - Z) ./ (1 + Z); 
RT_R0 = log(RT_R0);

% Generate the coefficients for this thermistor.
beta = zeros(1,order);
p = polyfit(RT_R0, T_ref_regress, order); 
pp = p; % save for later usage
p = 1 ./ p;
p = fliplr(p); % place in ascending order
T_0    = p(1);
for index = 2:order+1
    beta(index-1) = p(index);
end

% make a smooth line using coefficients (for comparison)
R = linspace(min(RT_R0), max(RT_R0), 1000);
R = R';
T_inverse_predicted = polyval(pp,R);

%-----------------------------------------------------
% --------- Plot regression --------------------------
%   The plot should be a nearly straight line.
%-----------------------------------------------------
if plot_regress && make_figures
    fig_num = fig_num+1;
    figure(fig_num), clf
    
    h = plot(...
        T_ref_regress,       RT_R0, '.', ...
        T_inverse_predicted, R,     'r');
    xlabel ('\itT \rm^{-1} [K^{-1}]')
    ylabel ('log_{\ite\rm} (\itR_T\rm / \itR_{\rm0} \rm)')
    set(h(1), 'markersize', 15)
    set(h(2), 'linewidth',   2)
    
    legend('Observed', 'Predicted','location', 'southeast')
    
    betatext = '';
    for ii = 1:order
        betatext = [betatext ', \beta_{\it',num2str(ii),'}\rm = ' num2str(beta(ii))];
    end
    title_string{2} = ['\itT_{\rm0}\rm = ' num2str(T_0) , betatext];
    title(title_string)
    x_limits = get(gca,'xlim');
    x_text = x_limits(1) + (x_limits(2) - x_limits(1))/25;
    y_limits = get(gca,'ylim');
    y_text = y_limits(2) - (y_limits(2) - y_limits(1))/10;
    if order == 1
        my_text = [...
            '$$\frac{1}{T} = \frac{1}{T_0} + ' ...
            '\frac{1}{\beta_1} \log_{\ e}\left(\frac{R_T}{R_0}\right) $$'];
    elseif order == 2
        my_text = [...
            '$$\frac{1}{T} = \frac{1}{T_0} + ' ...
            '\frac{1}{\beta_1} \log_{\ e}\left(\frac{R_T}{R_0}\right) + ' ...
            '\frac{1}{\beta_2} \log^{\ 2}_{\ e} \left(\frac{R_T}{R_0}\right)$$'];
    end
    
    text(x_text, y_text, my_text, 'interpreter','latex', 'fontsize', 16)
end

%------------------------------------------------------------------
% -- Use the computed co-efficients to estimate temperature -------
%------------------------------------------------------------------
% Estimate temperature values using computed calibration coeffients
T_calibrated = polyval(pp, RT_R0);
T_calibrated = 1 ./ T_calibrated;
T_calibrated = T_calibrated - 273.15;

if plot_result && make_figures
    fig_num = fig_num+1;
    figure(fig_num), clf
    fig_aspectratio(gcf,1.8);
    
    lagtext = [texstr(T_ref_string) ' Lag = ',num2str(Lag,'%4.2f'),' s'];
    
    
    ax(1) = subplot(1,2,1);
    plot(...
        T_ref_prof,  P_prof, ...
        T_calibrated, P_prof);grid on
    set(gca,'ydir','rev')
    xlabel('[ ^{\circ}C]')
    ylabel('\itP   \rm[dBar]')
    title({'Temperature Profile',lagtext},...
        'fontweight','normal')
    legend(...
        ['lagged ' texstr(T_ref_string)], ...
        texstr(T_string), ...
        'location','southeast')
    
    ax(2)=subplot(1,2,2);
    c = get(0,'defaultaxescolororder');
%     plot(T_calibrated(m) - T_ref(m-m_lag), P_slow(m),'color',c(4,:));grid on
    plot(T_calibrated - T_ref_prof, P_prof,'color',c(4,:));grid on
    set(gca,'ydir','rev')
    xlimits = get(gca,'xlim');
    xlim(max(abs(xlimits))*[-1 1])
    xlabel('[ ^{\circ}C]')
    legend(texstr([T_string  ' - lagged ' T_ref_string]),...
        'location','southeast')
    title({'\rm Temperature Difference',lagtext})
    linkaxes(ax,'y')
end
