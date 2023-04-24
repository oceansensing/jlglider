%% quick_look
% Visualize contents of a RSI data file, compute spectra for a selected
% range, and return a profile of the rate of dissipation of kinetic energy.
%%
% <latex>\index{Functions!quick\_look}</latex>
%
%%% Syntax
%   diss = quick_look (fname, P_start, P_end, ql_info, ...)
%
% * [fname] Name of the binary data file to process (extension optional).
% * [P_start] Starting point, in pressure, of the segment to be used for
%           displaying a spectrum of all shear, scalar-gradient, and
%           acceleration signals. Can be empty, if P_end is also empty, to
%           suppress the display of the spectra.
% * [P_end] End point, in pressure, for the segment used to show spectra. Can
%           be empty, if P_start is also empty.
% * [ql_info] Structure containing configuration parameters. A template
%           is generated and returned when quick_look is called with no
%           input parameters. The parameters are described below.
% * [...] Optional configuration parameters to supplement, or override,
%           those values included within ql_info.
% * []
% * [diss] Dissipation estimates for the specified profile. A default
%           ql_info structure is returned when quick_look is called with no
%           input parameters.
%
%%% Description
% This function generates a variety of figures that visualize the data in
% the file $\texttt{fname}$. The function works in three stages. In the
% first stage, data is converted into physical units using the
% $\texttt{odas\_p2mat}$ function, if a mat-file does not already exist. In
% the second stage, the data are plotted in several figures. In the third
% stage, the function computes a profile of the rate of dissipation of
% kinetic energy from the shear probe data.
%
% To convert the data into physical units, $\texttt{quick\_look}$ calls the
% $\texttt{odas\_p2mat}$ function. The function gets executed if a mat-file 
% with the same name as the data file does not exist, or, if the processing 
% parameters are not the same as those that were used to create the
% existing mat-file. You will be prompted for permission to overwrite the 
% mat-file. Therefore, input parameters used with $\texttt{quick\_look}$ 
% can include the parameters required to convert your data into physical 
% units, but these parameters are only required if the mat-file does not 
% exist.
%
% If $\texttt{quick\_look}$ is called without input arguments, it returns
% the default parameters used by $\texttt{quick\_look}$ and by
% $\texttt{odas\_p2mat}$, within a single structure. You can then customize
% this structure to your particular processing requirements. For example,
%
%    >> ql_info = quick_look
%
%    ql_info =
%                   HP_cut: 0.4000
%                   LP_cut: 30
%                     YD_0: 0
%                despike_A: [8 0.5000 0.0400]
%                despike_C: [10 1 0.0400]
%               despike_sh: [8 0.5000 0.0400]
%              diss_length: 8
%                     f_AA: 98
%                  f_limit: Inf
%               fft_length: 2
%                fit_2_isr: 1.5000e-05
%                fit_order: 3
%                  goodman: true
%             make_figures: 1
%                  op_area: 'open_ocean'
%                  overlap: 4
%             plot_battery: 0
%         plot_dissipation: 1
%          plot_kinematics: 1
%            plot_rawaccel: 0
%             plot_sensors: 1
%             plot_spectra: 1
%        plot_spectrograms: 0
%            profile_min_P: 1
%            profile_min_W: 0.2000
%     profile_min_duration: 20
%              profile_num: 1
%          MF_extra_points: 0
%                     MF_k: 4
%                 MF_k_mag: 1.7000
%                   MF_len: 256
%                MF_st_dev: 'st_dev'
%             MF_threshold: []
%                      aoa: []
%           constant_speed: []
%            constant_temp: []
%               hotel_file: []
%             speed_cutout: 0.0500
%                speed_tau: []
%              time_offset: 0
%                  vehicle: ''
%
% The configuration parameters (fields) within the structure
% $\texttt{ql\_info}$ that control the behaviour of $\texttt{quick\_look}$,
% are listed below. They are grouped for clarity, but are all part of the
% single structure $\texttt{ql\_info}$, including those required only by
% $\texttt{odas\_p2mat}$.
%
%%% Parameters that specify a profile
% A single data file may contain multiple profiles. For example, a vertical
% profiler that was raised and lowered multiple times while recording data
% into a single file. Profiles are detected within $\texttt{quick\_look}$
% using 5 specifications. These are:
%
% * [profile_num] Index to the requested profile from the set of detected
%       profiles.  The first profile, 1, is the default value.
% * [profile_min_P] The minimum pressure of a profile, [dbar]. Default = 1.
% * [profile_min_W] The minimum vertical speed of profiling, [dbar/s].
%      Default = 0.2. Use a smaller value, ~0.1 dBar/s, for gliders.
% * [profile_min_duration] The minimum duration in which the minimum
%      pressure and speed must be satisfied [s]. Default = 20.
%
% The $\texttt{direction}$ of profiling is the implict 5th specification
% and is determined by the $\texttt{vehicle}$ used for profiling. See the
% file $\texttt{default\_vehicle\_attributes.ini}$ for the direction of each
% recognized $\texttt{vehicle}$.
%
% <latex>\clearpage</latex>
%%% Parameters that control the calculation of dissipation rate
% The rate of dissipation of turbulent kinetic energy, $\epsilon$, is
% estimated using the function $\texttt{get\_diss\_odas}$ and follows the
% method described in RSI Technical Note 028. The controlling parameters
% are:
%
% * [diss_length] The time span, in seconds, over which to make each
%      estimate of the rate of dissipation. Default = 8.
% * [overlap] The overlap, in seconds of each dissipation estimate. Default
%      = 4.
% * [fft_length] The length, in seconds, of the fft-segments that will be
%      ensemble-averaged into a spectrum. Default = 2.
% * [fit_2_isr] The rate of dissipation, in W/kg, above which the function will
%      switch from the method of spectral-integration to the method of
%      fitting to the inertial subrange, to estimate the rate of
%      dissipation. Default = 1.5e-5.
% * [fit_order] The order of the polynomial to fit to the shear spectrum,
%      in log-log space, to estimate the wavenumber at which the spectrum
%      has a minimum. This is one of several constraints on the upper limit
%      of spectral integration. Default = 3.
% * [ f_limit] The upper frequency limit, in Hz, to be used for estimating
%      the rate of dissipation. Default = inf (no unconditional limit).
% * [f_AA] The cut-off frequency, in Hz, of the anti-aliasing filter in
%      your instrument. Default = 98. This value is instrument dependent
%      but is almost always 98, unless you have an instrument that has been
%      customized to samples at rates faster than 512 per second.
% * [goodman] A logical variable that controls whether or not the Goodman
%      noise removal algorithm is implemented.
%
% <latex>\clearpage</latex>
%%% Parameters that control the processing of microstructure signals
% The parameters that control the processing of the microstructure signals
% pertain to the high-pass filtering of the shear-probe signals, and the
% despiking of the shear-probe, acceleration and micro-conductivity signals.
%
% * [HP_cut] The cut-off frequency, in Hz, of the high-pass filter applied
%      to the shear-probe signals. Default = 0.4.
% * [despike_sh] The triplet of parameters for the despike
%      function applied to the shear-probe signals. The first value is the
%      threhold. The second value is the cut-off frequency, in Hz, of the
%      low-pass smoothng filter. The third value is the duration, in
%      seconds, of data to remove around a spike. Default = [8 0.5 0.04].
%      See the function despike for more information on the parameters. You
%      can suppress the despike function by specifying an infinite
%      threshold, for example [inf 0.5 0.07].
% * [despike_A] The triplet of parameters for the despike
%      function applied to the accelerometer signals. For data collected
%      with a glider, it may be necessary to supress despiking. The
%      intermittent vibrations from battery movement and fin actuators
%      creates short duration vibrations that are easily confused with
%      spikes, but such data is needed for coherent-noise removal. Default
%      = [8 0.5 0.04].
% * [despike_C] The triplet of parameters for the despike function applied
%      to the micro-conductivity signals. Default = [10 1.0 0.04].
%
% <latex>\clearpage</latex>
%%% Parameters for data visualization
% The parameters that control the data visualization include plotting flags
% and a low pass filter that is applied to the shear data to better see
% turbulent patches. Additional parameters control the axes labels and
% scales. 
%
% * [LP_cut] The cut-off frequency, in Hz, of the low-pass filter applied
%      to the microstructure profile signals, for graphical display only.
%      It does not affect the estimation of the rate of dissipation.
%      Default = 30.
% * [YD_0] The year-day subtracted from the time axis of figures. It is
%      currently not used. Default = 0.
% * [op_area] The operational area of your instrument. Recognized values are
%     'open_ocean' and 'tidal_ch'. It controls the
%     scale on certain figures. Default = 'open\_ocean'.
% * [make_figures] The parameter that determines if figures are generated.
%     make_figures = false suppresses the generation of figures to speed up
%     the data processing. Default = true.
% * [plot_battery] A logical parameter that determines if a figure of the 
%     battery voltage is generated. Default = false.
% * [plot_dissipation] A logical parameter that determines if a figure of the 
%     dissipation rate is generated. Default = true.
% * [plot_kinematics] A logical parameter that determines if figures of the 
%     instrument kinematics are generated. This includes the inclination
%     angles, the fall rate (i.e. dP/dt), the vehicle speed (if not a VMP)
%     and the accelerometers. Default = true.
% * [plot_rawaccel] A logical parameter that determines if the raw 
%     acceleration is plotted in addition to the despiked signal. 
%     Default = false.
% * [plot_sensors] A logical parameter that determines if figures of the 
%     oceanographic data are plotted. This includes microstructure data and
%     ancillary data (e.g. JAC-CT). Default = true.
% * [plot_spectra] A logical parameter that determines if the frequency and
%     wavenumber spectra are plotted for the specified input range. 
%     Default = true.
% * [plot_spectrograms] A logical parameter that determines if spectrograms
%     of the shear and accelerometer data are generated. Default = false.
%
%%% Parameters for the odas_p2mat function
% The parameters starting with $\texttt{MF\_extra\_points}$ are used by the
% $\texttt{odas\_p2mat}$ function, to convert your data into physical
% units. They are described in the section for that function.
%
%
%%% The diss structure output
% The output structure from $\texttt{quick\_look}$ depends slightly on the
% channels in your instrument. The typical fields are shown below, in groups.
% The first group is associated with the calculation of the profile of the
% rate of dissipation, $\epsilon$, and is described in the section for
% $\texttt{get\_diss\_odas}$.
%
%                         [e]  [2x217 double]
%                     [K_max]  [2x217 double]
%                    [method]  [2x217 double]
%                  [dof_spec]  15.2000
%                     [dof_e]  [2x217 double]
%                       [mad]  [2x217 double]
%                        [FM]  [2x217 double]
%              [Nasmyth_spec]  [513x2x217 double]
%                  [sh_clean]  [4-D double]
%                        [sh]  [4-D double]
%                        [AA]  [4-D double]
%                        [UA]  [4-D double]
%                         [F]  [513x217 double]
%                         [K]  [513x217 double]
%                     [speed]  [217x1 double]
%                        [nu]  [217x1 double]
%                         [P]  [217x1 double]
%                         [T]  [217x1 double]
%                         [t]  [217x1 double]
%                       [AOA]  []
%                 [Data_fast]  [18x217 double]
%                 [Data_slow]  [19x217 double]
%                      [f_AA]  88.2000
%                   [f_limit]  Inf
%                 [fit_order]  3
%               [diss_length]  4096
%                   [overlap]  2048
%                [fft_length]  1024
%
% The next group is associated with the despiking of the shear-probe,
% micro-conductivity and acceleration signals.
%
%                  [spikes_A]  {[4x1 double], []}
%              [pass_count_A]  [2x1 double]
%                [fraction_A]  [2x1 double]
%                 [spikes_sh]  {[247x1 double], [269x1 double]}
%             [pass_count_sh]  [2x1 double]
%               [fraction_sh]  [2x1 double]
%                  [spikes_C]  {}
%              [pass_count_C]  [0x1 double]
%                [fraction_C]  [0x1 double]
%
% The indices to the spikes located in the signals from the accelerometers,
% are given in $\texttt{spikes\_A}$. The number of passes of the despike
% function used to remove the spikes is in $\texttt{pass\_count\_A}$. The
% fraction of data removed by the despike function is in
% $\texttt{fraction\_A}$. Similarly for the shear probe and the
% micro-conductivity signals.
%
% The next group is associated with the scalar signals.
%
%            [scalar_spectra]  [1x1 struct]
%        [scalar_vector_list]  {'gradT1'  'gradT2'}
%               [scalar_info]  [1x1 struct]
%
% The structure
% $\texttt{scalar\_spectra}$ (a structure within a structure) is described in
% the section for the function $\texttt{get\_scalar\_spectra\_odas}$. The
% processing parameters are in $\texttt{scalar\_info}$. The names of the
% scalar vectors are in $\texttt{scalar\_vector\_list}$.
%
% The remaining fields are:
%
%                   [fs_fast]  511.9454
%                   [fs_slow]  63.9932
%                     [piezo]  1
%     [profiles_in_this_file]  1
%              [speed_source]  'Rate of change of pressure'
%                 [fast_list]  {1x18 cell}
%                 [slow_list]  {1x19 cell}
%                   [ql_info]  [1x1 struct]
%                [ql_info_in]  []
%
% $\texttt{fs\_fast}$ and $\texttt{fs\_slow}$ are the actual fast and slow
% sampling rates of the data. The number of profiles detected in this data
% file is given in $\texttt{profiles\_in\_this\_file}$. The source of the
% speed of profiling is identified in $\texttt{speed\_source}$.
%
% All column vectors that have a length matching the length of the time
% vector $\texttt{t\_fast}$ are combined into a single matrix and passed to
% the dissipation function for averaging over the interval of each
% dissipation estimate. Similarly for all column vectors that match
% $\texttt{t\_slow}$. Each row of $\texttt{Data\_fast}$ and
% $\texttt{Data\_slow}$ hold the values from a single vector. The names of
% the signals are identified in the cell arrays $\texttt{fast\_list}$ and
% $\texttt{slow\_list}$. In this example there are 19 fast vectors and 19
% slow vectors. Use the Matlab $\texttt{find}$ and $\texttt{strcmpi}$
% functions to identify the row of a particular signal.
%
% The structure that was input to this function is saved in field
% $\texttt{ql\_info\_in}$ and is empty in this example because the call
% used default values. The structure that was actually used to process
% the file, in this case a structure of default values, is given in the
% field $\texttt{ql\_info}$.

% * 2015-10-31 RGL Update documentation
% * 2015-11-04 RGL Changed position of legend on Fig 1 for gliders.
%      Supressed spectrogram outputs. Improved handling of pressure ranges
%      that are empty or do not exist. Changed CTD figure.
% * 2015-11-09 RGL Amended parameters associated with despike.
% * 2015-11-12 RGL Changed scaling on piezo-accelerometers for case of
%      Sea-glider. Changed profile_P_min to profile_min_P for consistency
%      with other profile parameters. Change profile_min_speed to
%      profile_min_W to make it more explictly the vertical component of
%      speed. Signicantly changed the way we handle accelerometers by
%      distinquishing between piezo and linear accelerometers, in a maner
%      that is backwards compatible with older data files that use
%      type=accel instead of type=piezo. Removed the ability to make pdf
%      fileas of the figures. Added option to supress making figures for
%      greater speed of execution.
% * 2015-11-18 RGL Force drawnow after every figure.
% * 2015-11-20 Changed recognisition of piezo-accelerometers to accomodate
%      changes to setupstr. Removed BP shear from figure for tidal_channel
%      data.
% * 2015-12-07 Check for existance of P_limit before using it. Required by
%      XMPs.
% * 2015-12-23 Fixed frequency spectra plot. Indices were erroneous so the
%      wrong data range was being plotted.  Some minor cleaning of the
%      code.
% * 2016-01-21 Removed third linear accelerometer plot from Figure 1.
% * 2016-04-19 RGL, added factor of 10 scaling of mico-conductivity
%           gradient for case of "tidal_ch". Corrected legend bug --
%           numbers are no longer supported for location of legend. Fixed
%           bug with despiking the piezo-accelerometers. Now the statistics
%           are only for the profile and not the entire file.
% * 2016-04-25 RGL, added factor of 10 scaling of thermistor temperture
%           gradient for case of "tidal_ch".
% * 2016-05-28 RGL, changed thermistor correction to double-pole from
%           single-pole to make it consistent with Vachon and Lueck, 1984.
% * 2016-07-22 RGL, added salinity_JAC for the case that we have a JAC_CT.
%           We still use salinity function in case of Sea-Bird CT.
% * 2016-08-31 RGL, added vehicle_info and params to the returned structure
%           so that the direction of profiling can be detected. show_spec
%           will now indicate the time since the start of the file instead
%           of preesure for the case of horizontal profiles.
% * 2016-09-06 RGL, The filffilt function is still causing problems with
%           the high-pass filtering of the shear-probe signals,
%           particularly when there is very strong shear near the end of a
%           profile. So, I have reverted to using the filter function
%           forwards and backwards. This was noticed in the Arnoldo data
%           collected in Florida.
% * 2016-11-08 RGL, Corrected the indexing problem associated with the
%           linear accelerometers. Changed the way that the linear
%           accelerometers are plotted.
% * 2016-11-14 RGL, Allowed xticklabel in case of linear accelerometers.
% * 2016-12-12 RGL, Trimmed all vectors to the range of a profile to
%           simplify processing. Fluorometer, Backscatter, Suna, etc can
%           now be fast or slow channels. Linear accelerometers can now be
%           fast or slow. However, piezo-accelerometers must be fast. Added
%           check in Goodman routine to make sure that the sample rate is
%           fast.
% * 2016-12-13 RGL, Added check for electromagnetic current meter reference
%           signal, and use it to clean the shear-probe and scalar signals.
%           Coherent noise removal can now be controlled by the usage of
%           the logical variable "goodman". The default is "true".
% * 2016-12-13 RGL, Using the who-function to simplify finding the
%           variables for plotting.
% * 2016-12-14 RGL, More streamlining of code wrt gnerating variable names.
% * 2016-12-15 RGL, More of the same. Increased efficiency of extracting
%           the accelerometer signals by using pre-assignment. Also, big
%           mistake in the piezo-accelerometer despiking -- we despiked the
%           AA_piezo matrix but not the AA matrix, which is later used for
%           coherent noise removal. So, we never used the despiked
%           accelerometer signals. Stupid.
% * 2016-12-22 RGL, added support of gradT and gradC methods. Fixed bug
%           with scalar spectra of selected depth range. The profile of
%           spectra was correct, but not the specified depth range.
% * 2017-06-09 RGL, Corrected naming of JAC EM current meter. Now using
%           EMC_Cur instrad of EM_Cur.
% * 2017-06-22 RGL, Corrected naming accelerometers in the spectrograms.
% * 2017-07-18 RGL, Check if EMC_Cur is fast or slow. Adjust accel_count
%           down by 1 if EMC_Cur exists.
% * 2017-07-26 RGL, Fixed major bug -- piezo-accelerometers were not
%           being despiked. The despiked matrix was placed into AA_piezo
%           instead of AA.
% * 2018-01-05 JMM, Minor bug fixes
% * 2018-01-16 JMM, Major changes to simplify and streamline quick_look for
%           version 4.3 release. Changes include:
%               - a flag to only generate the basic figures
%               - only plot W_2 in kinematics plot
%               - added U_EMC to kinematics plots when appropriate
%               - set axes limits when acceleration counts are really high
%               - separating plots into subpanels
%               - improved comments
% * 2018-03-14 JMM, Added a plot of the pressure for the entire file
%           highlighting the profile that is being analyzed. Also modified
%           the plots of the spectra to be in one subpanel. Also some minor
%           formatting changes and cleaning up of the code. 
% * 2018-03-19 JMM, Fixed a plotting bug that arose when using a constant
%           speed
% * 2018-05-11 JMM, Fixed call to salinity_jac to include sampling rate and
%           mean speed            
% * 2018-05-18 HG/JMM, Added functionality to support horizontal profilers. 
%           Also added plots for magnetometer and rotation rate added. Code
%           was copied from an earlier version created by Helen Gemmrich
% * 2018-09-27 JMM, Small fix to accelerometer plot for seagliders (no
%           longer sets axes limits)
% * 2019-06-05 JMM, Added if statements to handle MicroPods when there is
%           shear is the only microstructure sensor 
%           Removed blank subpanels in kinematics plots
% * 2019-07-03 JMM, Modified offset for JAC sensors so that it is correctly 
%           plotted as a pressure lag. Note: The time lag of 0.14 s that is
%           used is only an approximation. The plots that have CTD data
%           also now plot C_CTD_match, which is the conductivity signal
%           that was matched to the response of the temperature sensor. 
% * 2019-07-04 JMM, Minor formatting fixes of plots. 
% * 2019-08-16 JMM, No longer switch inclinometer channels for large
%           angles, unless the instrument is a vmp or rvmp
% * 2019-09-24 JMM, Bug fix: Suppress plotting of raw dP/dt and only show
%           low pass filtered signal.
% * 2019-10-07 JMM, Now plot both despiked and raw vibration sensor data by 
%           default.
% * 2019-10-08 JMM, Bug fix: When linear accelerormeters and 
%           inclinometers were not sampled at the same rate, a mis-match 
%           of vector sizes was causing the code to crash. 
%           Axes labels for inclinometers were also reversed.
% * 2020-09-09 JMM, Bug fix: Corrected the linewidths for the vibration
%           sensor plots to have consistency between the two panels.





function result = quick_look( fname, P_start, P_end, varargin )

%-----------------------------------------------------------------
% ----- Default parameters ---------------------------------------
%-----------------------------------------------------------------
% -- Deployment params
default_op_area       = 'open_ocean';   % Deployment area (Adjust some figure scales)
default_YD_0        = 0;                % Deployment day [year day]

% -- Flags to make figures
default_make_figures      = true;    % render figures for data visulization
default_plot_kinematics   = true;    % flag to plot kinematics
default_plot_rawaccel     = true;   % flag to plot unspiked accelerometer data
default_plot_sensors      = true;    % flag to plot raw data from sensors
default_plot_battery      = false;   % flag to plot battery voltage
default_plot_spectra      = true;    % flag to plot spectra
default_plot_dissipation  = true;    % flag to plot dissipation
default_plot_spectrograms = false;   % flag to plot spectrograms

% -- Values used to determine profiles
default_profile_num          = 1;      % profile number
default_profile_min_P        = 1;      % minimum pressure [dBar] to identify profile
default_profile_min_W        = 0.2;    % minimum speed [m/s] to identify pressure
default_profile_min_duration = 20;     % minimum duration [s] for a segment to be considered a profile

% -- Values to compute spectra and dissipation rates
default_fft_length  = 2;                            % length [s] of each fft segment
default_diss_length = 4* default_fft_length;        % length [s] for each dissipation estimate
default_overlap     = round(default_diss_length/2); % overlap [ratio] between dissipation estimates
default_HP_cut      = 0.4;                          % frequency [Hz] of high-pass filter for shear probe data
default_LP_cut      = 30;                           % frequency [Hz] of low-pass  filter for shear probe data (FOR PROFILE DISPLAY PURPOSES ONLY!)
default_fit_2_isr   = 1.5e-5;                       % dissipation rate [W/kg] above which fit-to-ISR method is used
default_fit_order   = 3;                            % order of the polynomial used to compute upper limit for integration
default_f_AA        = 98;                           % frequency [Hz] of anti-aliasing filter
default_f_limit     = inf;                          % frequency [Hz] limit for integration of spectra
default_goodman     = true;                         % flag to implement Goodman noise removal algorithm

% -- Despiking parameters  [thresh, smooth, and length] (in seconds)
default_despike_sh  = [ 8  0.5 0.04];               % for shear probes
default_despike_A   = [ 8  0.5 0.04];               % for piezo-accelerometers
default_despike_C   = [10  1.0 0.04];               % for micro-C

% -- Return defaults if no inputs to function are specified
if ~nargin,
    for d = whos('default_*')',
        param = regexp(d.name, 'default_(.*)', 'tokens');
        result.(char(param{1})) = eval(d.name);
    end
    % -- Add the arguments from odas_p2mat (Not required but useful as a reference)
    p2mat = odas_p2mat();
    for name = fieldnames(p2mat)'
        result.(name{1}) = p2mat.(name{1});
    end
    return
end

%-----------------------------------------------------------------
% ----- Parse Inputs ---------------------------------------------
%-----------------------------------------------------------------
p = inputParser;
p.KeepUnmatched = true;
p.CaseSensitive = true;

val_numeric     = @(x) isnumeric(x) && isscalar(x);
val_speed       = @(x) isnumeric(x) && isscalar(x)   || isempty(x);
val_vector      = @(x) (isnumeric(x) && isvector(x)) || isempty(x);
val_string      = @(x) ischar(x);
val_logical     = @(x) islogical(x);

addRequired(  p, 'fname',   val_string);
addRequired(  p, 'P_start', val_speed);
addRequired(  p, 'P_end',   val_speed);

addParamValue(p, 'YD_0',                 default_YD_0,                 val_numeric);
addParamValue(p, 'op_area',              default_op_area,              val_string);

addParamValue(p, 'make_figures',         default_make_figures,     val_logical);
addParamValue(p, 'plot_kinematics',      default_plot_kinematics,  val_logical);
addParamValue(p, 'plot_rawaccel',        default_plot_rawaccel,    val_logical);
addParamValue(p, 'plot_sensors',         default_plot_sensors,     val_logical);
addParamValue(p, 'plot_battery',         default_plot_battery,     val_logical);
addParamValue(p, 'plot_spectra',         default_plot_spectra,     val_logical);
addParamValue(p, 'plot_dissipation',     default_plot_dissipation, val_logical);
addParamValue(p, 'plot_spectrograms',    default_plot_spectrograms,val_logical);

addParamValue(p, 'profile_num',          default_profile_num,          val_numeric);
addParamValue(p, 'profile_min_P',        default_profile_min_P,        val_numeric);
addParamValue(p, 'profile_min_W',        default_profile_min_W,        val_numeric);
addParamValue(p, 'profile_min_duration', default_profile_min_duration, val_numeric);

addParamValue(p, 'fft_length',           default_fft_length,           val_numeric);
addParamValue(p, 'diss_length',          default_diss_length,          val_numeric);
addParamValue(p, 'overlap',              default_overlap,              val_numeric);
addParamValue(p, 'fit_order',            default_fit_order,            val_numeric);
addParamValue(p, 'HP_cut',               default_HP_cut,               val_numeric);
addParamValue(p, 'LP_cut',               default_LP_cut,               val_numeric);
addParamValue(p, 'fit_2_isr',            default_fit_2_isr,            val_numeric);
addParamValue(p, 'f_limit',              default_f_limit,              val_numeric);
addParamValue(p, 'f_AA',                 default_f_AA,                 val_numeric);
addParamValue(p, 'goodman',              default_goodman,              val_logical);

addParamValue(p, 'despike_sh',           default_despike_sh,           val_vector);
addParamValue(p, 'despike_A',            default_despike_A,            val_vector);
addParamValue(p, 'despike_C'  ,          default_despike_C,            val_vector);

% -- save the input for record keeping
ql_info_in = join_arguments(varargin);

% -- Parse the arguments.
parse(p, fname, P_start, P_end, varargin{:});

% -- Validate input
if p.Results.diss_length < 2*p.Results.fft_length
    error('Invalid size for diss_length - must be greater than 2 * fft_length.');
end
if p.Results.P_end < p.Results.P_start,
    error('Starting pressure must be less than end pressure.');
end

% -- create ql_info structure (uses both user specified values and defaults)
names = fieldnames(p.Results);
for name = names'
    eval(['ql_info.' char(name) ' = p.Results.' char(name) ';']);
end
names = fieldnames(p.Unmatched);
for name = names'
    eval(['ql_info.' char(name) ' = p.Unmatched.' char(name) ';']);
end

% -- Define simple variables to use throughout function
names = fieldnames(p.Results);
for nn = 1:length(names)
    if ~strcmp(names{nn},'fname') || ~strcmp(names{nn},'P_start') || ~strcmp(names{nn},'P_end')
        eval([names{nn} '= p.Results.' names{nn} ';'])
    end
end
names = fieldnames(p.Unmatched);
for nn = 1:length(names)
    if ~strcmp(names{nn},'fname') || ~strcmp(names{nn},'P_start') || ~strcmp(names{nn},'P_end')
        eval([names{nn} '= p.Unmatched.' names{nn} ';'])
    end
end

% -- cleanup
YD_0 = floor(YD_0); 
    
clear default*
set(groot, 'defaultFigureColor' ,'White');

% -----------------------------------------------------------------------
% --- Convert to physical units ----------------------------------------
% -----------------------------------------------------------------------

% -- Call odas_p2mat and create mat file (if needed)
d = odas_p2mat(fname, ql_info);
[P,N,E] = fileparts(d.fullPath);
File_Name = [N E];
if ~exist([P filesep N '.mat'], 'file')
    disp(['Saving into MAT-file: ' P filesep N '.mat']);
    save([P filesep N '.mat'], '-struct', 'd', '-v6');
else
    disp(['Loading from MAT-file: ' P filesep N '.mat']);
end

% -- Extract fields from data structure into simple variables
for field = fieldnames(d)'
    eval([char(field) ' = d.' char(field) ';']);
end

% -- Rename Current Meter Vectors
if exist('U_slow','var'), U = U_slow; end
if exist('V_slow','var'), V = V_slow; end
if exist('W_slow','var'), W = W_slow; end
if exist('U_EMC','var'),  U_EM = U_EMC; end % Handles misnaming of EM channel
if exist('EM_Cur', 'var'), EMC_Cur = EM_Cur; end % Handles misnaming of EM channel

% -- Extract variables from "params" structure (returned by odas_p2mat)
% -------
% Note: If a variable of the same name already exists, keep it.  This 
% ensures that parameters passed into quick_look are used in place of 
% those parameters used when generating the MAT file.
%--------
names = fieldnames(params);
for name = names'
    if ~exist(name{1}, 'var')
        eval([name{1} ' = params.' name{1} ';']);
    end
end

% -- Define additional parameters
model       = char(setupstr(cfgobj,'instrument_info', 'model'));
serial_num  = char(setupstr(cfgobj,'instrument_info', 'serial_num'));
profile_dir = char(setupstr(cfgobj,'instrument_info', 'profile_dir'));
if isempty(profile_dir), profile_dir = vehicle_info.profile_dir; end

fft_num = round(fft_length*fs_fast);% fft length in units of samples
tidal_channel = strcmpi(op_area,'tidal_ch'); % logical variable to control scale on figures

% -- Determine if the vehicle is an HMP
isHMP = false;
if (strcmpi(profile_dir, 'horizontal')), isHMP = true; end

%------------------------------------------------------------------------
% ----- Define Independent Variables (VMP/glider: pressure, HMP: time)---
%------------------------------------------------------------------------
if isHMP
    Indep_var_slow = t_slow - YD_0;
    Indep_var_fast = t_fast - YD_0;
    Indep_labelstr = ['time [s]'];
else
    Indep_var_slow = P_slow;
    Indep_var_fast = P_fast;
    Indep_labelstr = ' \itP\rm [dbar]';
end
%-----------------------------------------------------------------
% ----- Figure Properties  ---------------------------------------
%-----------------------------------------------------------------

fig_num = 0;

% -- define colors
fuj  = [0.75 0    0.5]; % a shade of fujia
gold = [0.75 0.65 0];
green = [0 0.9 0];
cyan  = [0 0.9 0.9];
brown = [0.54 0.27 0.07];
figure(1),clf
ColOrd = get(gca,'ColorOrder');
close




%-----------------------------------------------------------------
% ----- Extract a profile  ---------------------------------------
%-----------------------------------------------------------------
profile = [1 ; length(t_slow)]; % start using entire file

if exist('P_slow','var') && exist('W_slow','var')
    
    if strcmpi(profile_dir, 'up') || strcmpi(profile_dir, 'down')
        % get profile based on direction, duration, etc.
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
        % sort columns in ascending order
        profile = sort([profile_down profile_up],2);
    end
end

% - Quit if there are no profiles in this file
profiles_in_this_file = size(profile,2);
if profile_num > profiles_in_this_file
    warning(['There are only ' num2str(profiles_in_this_file) ' profiles in this file'])
    result = [];
    
    % Plot the pressure record
    fig_num = fig_num + 1;
    figure (fig_num); clf
    clear sp
    
    plot(t_slow,P_slow);
    
    % - properties and labels
    set(gca,'ydir','rev')
    xlabel('\it t \rm [s]')
    ylabel('\it P\rm [dbar]')
    
    % - set title
    title({['\rm ',texstr(fname)],...
        ['Requesting profile ',num2str(profile_num), ...
        ' but ' num2str(profiles_in_this_file) ' profiles detected in this file']})
    
    drawnow
    
    return
end

start_index_slow = profile(1, profile_num);
end_index_slow = profile(2, profile_num);

start_index_fast = 1 + round((fs_fast/fs_slow)*(start_index_slow - 1));
end_index_fast =     round((fs_fast/fs_slow)*(  end_index_slow    ));
n = (start_index_fast : end_index_fast)';
m = (start_index_slow : end_index_slow)';

%-----------------------------------------------------------------
% ----- Plot of the section of data being analyzed ---------------
%-----------------------------------------------------------------
if make_figures
    fig_num = fig_num + 1;
    figure (fig_num); clf
    clear sp
    
    % - get limits for future plots 
    plot(t_slow(m),P_slow(m))
    if ~isHMP
        y_limits = get(gca,'ylim'); 
    else
        y_limits = get(gca,'xlim');
    end
    
    % - plot
    p=plot(t_slow,P_slow,t_slow(m),P_slow(m));
    
    % - properties and labels
    set(p(1),'linewidth',1)
    set(p(2),'linewidth',2)
    set(gca,'ydir','rev')
    xlabel('\it t \rm [s]')
    ylabel('\it P\rm [dbar]')   
        
    % - basic title (time extracted from header of the starting profile)
    hd = header(ceil((start_index_slow - 1) / fs_slow + 1),:);
    line1 = sprintf('%s;  %d-%02d-%02d %02d:%02d:%02dZ', ...
                 File_Name, hd(4), hd(5), hd(6), hd(7), hd(8), hd(9));
    line2 = sprintf('Profile = %02d', profile_num);
    title_string = texstr({line1, line2});
    title_string{1} = ['\rm' title_string{1}]; 
    
    % - set title
    title(title_string)
    
    drawnow
end
    
%-----------------------------------------------------------------
% ----- Trim fast and slow data files to requested profile -------
%-----------------------------------------------------------------
% Develop the lists of fast and slow vectors in this data file. Only 
% recognizes column vectors of the right length.
% -----------------------------------------------------------------------

% - All variables
S = whos; 

% - Get sizes and indices
junk = zeros(size(S,1),2);
for k = 1:size(S,1)
    junk(k,:) = S(k).size; % Get the size of every variable
end
index_fast = find((junk(:,1) == size(t_fast,1)) & junk(:,2) == 1); % index_fast points to the fast channels.
index_slow = find((junk(:,1) == size(t_slow,1)) & junk(:,2) == 1); % index_slow points to the slow channels.

% - Trim vectors to the profile segment
for k = 1:size(index_fast,1)
    fast_name = S(index_fast(k)).name;
    eval([fast_name '=' S(index_fast(k)).name '(n);'] ); % trim to n
end
for k = 1:size(index_slow,1)
    slow_name = S(index_slow(k)).name;
    eval([slow_name '=' S(index_slow(k)).name '(m);'] ); % trim to m
end

clear junk

%-----------------------------------------------------------------
% ----- Calculate profile-specific estimate of the fall rate -----
%-----------------------------------------------------------------
% The variable 'W_slow' is heavily biased by the zero speed after the 
% profile. It was calculated using a 4th order zero-phase filter that 
% "sees" the stop after the profile. Here W_slow is only estimated 
% from the data for the requested profile.
% -------------------------------------------------------------------
fc = 0.5;
[b,a]=butter(4,fc/(fs_slow/2));
dP_dt = gradient(P_slow,1/fs_slow);
dP_dt_LP = filtfilt(b,a,dP_dt);

%-----------------------------------------------------------------
% ----- Figure: Kinematics (i.e. vehicle behavior) ---------------
%-----------------------------------------------------------------
% Speed, instrument roll, rate of change of pressure  and XY
% inclination (gliders only) plotted as subplots on the same figure. 
% Note: Won't execute if inclinometers are not present
% ------------------------------------------------------------------

if make_figures && plot_kinematics && exist('Incl_X','var') && exist('Incl_Y','var')
    
    %-- inclinometer delay
    fs_ADIS = 482; % Internal sampling rate of the ADIS Inclinometer
    tau_N = 256; % The number of samples in the running average internal to the ADIS
    ADIS_delay = (tau_N/2) / fs_ADIS; % time delay of output from ADIS (delta t)
    ADIS_P_delay =  dP_dt_LP * ADIS_delay; % pressure delay (delta P = W * delta t)
    
    %-- filter inclinometer data
    LP_ADIS = 2; % low-pass filter applied to ADIS Inclinometer
    [b,a] =   butter(1,LP_ADIS/(fs_slow/2));% low-pass to show only the gravity signal
    X_data = [Incl_X Incl_Y];
    X_data = filtfilt(b,a,X_data);
    
    %-- check if inclinometer channels are reveresed
    if abs(mean(Incl_X)) > 50 && (strcmp(vehicle,'vmp') || strcmp(vehicle,'rvmp'))
        X_data = fliplr(X_data);
    end   
    
    %-- plot
    fig_num = fig_num + 1;
    figure (fig_num); clf
    clear sp
    
    MM = 2; % default is two subplots
%     if ~isempty(constant_speed) || strcmpi(profile_dir,'glide')
    if strcmpi(profile_dir,'glide')
        MM = MM + 1;
        fig_aspectratio(gcf,2);
    end
    
    %-- subplot of inclinometer data
    sp(1) = subplot(1,MM,1);
    
    if ~strcmpi(profile_dir,'glide') && ~strcmpi(profile_dir,'horizontal') % VMP - plot only theta X
        h = plot(X_data(:,1), Indep_var_slow + ADIS_P_delay,'linewidth',2); grid on
        xlabel('\theta_X [  ^{\circ} ]')
    else
        if strcmpi(profile_dir,'glide')
            h = plot(X_data, Indep_var_slow + ADIS_P_delay,'linewidth',2); 
        elseif strcmpi(profile_dir,'horizontal') % don't need delay
            h = plot(X_data, Indep_var_slow,'linewidth',2); 
        end
        grid on
        xlabel('[  ^{\circ} ]')
        lgnd = legend('\theta_X', '\theta_Y','location','northeast');
        lgnd_pos = get(lgnd,'Position');
        set(lgnd,'Position',[1.1*lgnd_pos(1) lgnd_pos(2:4)])
    end
    set(gca, 'ydir', 'rev')
    ylabel(Indep_labelstr)
    if MM == 2;  title(title_string); end
    
    %-- subplot of fall rate
    if ~strcmpi(profile_dir,'horizontal')
        sp(2) = subplot(1,MM,2);
        h = plot(dP_dt_LP, Indep_var_slow, 'linewidth', 2 ); grid on
%         hold all
%         h = plot(dP_dt, Indep_var_slow, 'linewidth', 2 ); grid on

        set(gca,'ylim',y_limits)
        set(gca, 'ydir', 'rev')
        xlabel('d\itP\rm / d\itt \rm[ dBar s^{-1}]')
    end
    if MM == 3;  title(title_string); end
    
    %-- subplot of estimated speed
    if (MM == 3 && strcmpi(profile_dir,'glide')) || strcmpi(profile_dir,'horizontal')
        sp(MM) = subplot(1,MM,MM);
        h = plot(speed_slow, Indep_var_slow ); grid on
        set(h(1), 'linewidth', 2)
        set(gca,'ylim',y_limits)
        set(gca, 'ydir', 'rev')
        xlabel('estimated U [ m s^{-1}]')
        if strcmp(speed_source,'Rate of change of pressure and glide angle')
            ll{1} = 'dP/dt, aoa';
        elseif strcmp(speed_source,'Hotel file')
            ll{1} = 'Hotel File';
        elseif strcmp(speed_source,'Forced to constant speed by input parameter')
            ll{1} = 'Constant Speed';
        elseif strcmp(speed_source,'Recorded U V W')
            ll{1} = 'Recorded U V W';
        else
            ll{1} = '';
        end
        legend(h,ll)
    end
    %-- add EM current meter, if available
%     if MM == 3 && exist('U_EM','var') 
    if exist('U_EM','var') 
        hold all
        h(2) = plot(U_EM,Indep_var_slow,'linewidth',1);
        ll{2} =texstr('U_EM');
        uistack(h(1),'top');
        legend(h,ll)
    end
    
    
    %-- show difference between fall rate and speed (if const speed is used)
    if ~isempty(constant_speed) && (strcmp(profile_dir,'up') || strcmp(profile_dir,'down'))
        sp(3) = subplot(1,MM,3);
        h = plot(constant_speed-dP_dt_LP, Indep_var_slow ); grid on
        set(h(1), 'linewidth', 2,   'color','b')
        set(gca,'ylim',y_limits)
        set(gca, 'ydir', 'rev')
        xlabel('dP/dt - W_0 [ m s^{-1}]')
        ll{1} = 'Constant Speed';
        legend(h,ll)
    end
    linkaxes(sp,'y')
%     if MM == 3; legend(h,ll); end
    
    drawnow    
end


%-----------------------------------------------------------------
% -----  Figure: Kinematics for an XMP ---------------------------
%-----------------------------------------------------------------
%  Plot the pitch and fall rate of an XMP - a separate plot is needed
%  because inclinometers are not used on this instrument
% -----------------------------------------------------------------------
if exist('Pitch','var') && make_figures && plot_kinematics
        
    %-- low pass filter pitch
    f_low_pass = 1;
    [b,a] = butter(2, f_low_pass /(fs_slow/2));
    X_data = Pitch;
    X_data = filtfilt(b,a,X_data);
    
    %-- plot
    fig_num = fig_num + 1;
    figure (fig_num); clf
    clear sp
    
    %-- subplot of pitch
    sp(1)=subplot(1,2,1);
    h = plot(X_data, Indep_var_slow); grid on
    set(h(1), 'linewidth', 2,   'color','b')
    set(gca, 'ydir', 'rev')
    ylabel(Indep_labelstr)
    xlabel('Pitch [  ^{\circ} ]')
    title(title_string)
    
    %-- subplot of fall rate
    sp(2) = subplot(1,2,2);
    h = plot(dP_dt_LP, Indep_var_slow ); grid on
    set(h(1), 'linewidth', 2,   'color','b')
    set(gca, 'ydir', 'rev')
    set(gca,'yticklabel',[])
    xlabel('d\itP\rm / d\itt \rm[ dBar s^{-1}]')
    
    linkaxes(sp,'y')
    
    drawnow
end

%--------------------------------------------------------------------------
% ----- Accelerometers ----------------------------------------------------
%--------------------------------------------------------------------------
% Identify if we have the true linear (DC) response accelerometers, or if
% we have piezo-accelerometers. A lot of data has been collected using type
% = accel for piezo acceleromters, with coefficients of 0 and 1. On
% 2015-11-12 we introduced the type=piezo to make it easier to handle both
% types of accelerometers, especially with respect to scaling the spectra of
% acceleration. This used to be a nightmare of if-statements. Here we
% figure out what type of accelerometers we have (including the possibility
% of both types, such as with the Nemo system), and place the data into
% matrices with appropriate names. In addition, we place AA_piezo into AA
% for coherent-noise removal. If the are no piezo acceleromters, then we
% place the linear accelerometers, AA_linear, into the matrix AA.
% -----------------------------------------------------------------------

%-- identify accelerometers
AA_names_linear = {};
AA_names_piezo  = {};
for ch = setupstr(cfgobj, '', 'type', 'accel|piezo')
    name = setupstr(cfgobj, ch{1}, 'name');
    type = setupstr(cfgobj, ch{1}, 'type');
    
    if strcmp(type, 'accel')
        AA_names_linear{end+1} = name{1};
    else
        AA_names_piezo{end+1} = name{1};
    end
end
AA_names_linear = sort(AA_names_linear);
AA_names_piezo  = sort(AA_names_piezo);

% --Build the matrix of linear accelerometers
L1 = 0;
L2 = 0;
if ~isempty(AA_names_linear)
    eval(['L1 = length(' AA_names_linear{1} ');'])
    L2 = length(AA_names_linear);
end
AA_linear = zeros(L1, L2);
for index = 1 : length (AA_names_linear)
    AA_linear(:, index) = eval(AA_names_linear{index}) ;
end

%-- Determine whether linear acclerometers are fast or slow channels
AA_linear_is_fast = [];
if size(AA_linear,1) == length(t_fast)
    AA_linear_is_fast = true;
    fs_AA = fs_fast;
end
if size(AA_linear,1) == length(t_slow)
    AA_linear_is_fast = false;
    fs_AA = fs_slow;
end

% --Build the matrix of piezo accelerometers
L1 = 0;
L2 = 0;
if ~isempty(AA_names_piezo)
    eval(['L1 = length(' AA_names_piezo{1} ');'])
    L2 = length(AA_names_piezo);
end
AA_piezo = zeros(L1, L2);
for index = 1 : length (AA_names_piezo)
    AA_piezo(:, index) = eval(AA_names_piezo{index}) ;
end

piezo = ~isempty(AA_piezo);

% -- Create AA matrix to be used for noise removal
%     (using piezo accelerometers, if available)
AA = [];
if piezo
    AA = AA_piezo;
elseif ~isempty(AA_linear)
    AA = AA_linear;
end

%------------------------------------------------------------------------
% ----- Figure: Kinematics with linear accelerometers -------------------
%------------------------------------------------------------------------
% Plot inclination and fall rate for instuments without inclinometers
% (estimate inclination from accelerometer data)
% -----------------------------------------------------------------------
if ~isempty(AA_linear) && make_figures && plot_kinematics && ~isempty(AA_linear_is_fast)
    
    %-- low pass filter and estimate inclination
    f_low_pass = 1; % in Hz
    [b,a] = butter(2, f_low_pass /(fs_AA/2));
    X_data = filtfilt(b,a,AA_linear(:,1:min(2,end))); % because we want rotation around y-axis
    X_data = asind(X_data / 9.81); % inclination [degrees]
    
    %-- appropriate independent variable
    if AA_linear_is_fast
        y_data = Indep_var_fast;
    else
        y_data = Indep_var_slow;
    end
    
    %-- plot
    fig_num = fig_num + 1;
    figure (fig_num); clf
    fig_aspectratio(gcf,2);
    clear sp
    
    %%-- subplot acceleration
    sp(1) = subplot(1,3,1);
    h = plot(AA_linear, y_data,'linewidth',2); grid on
    set(gca,'ylim',y_limits)
    set(gca, 'ydir', 'rev')
    xlabel('[m s^{-2}]')
    legend(AA_names_linear,'location','southeast');  
    ylabel(Indep_labelstr)
    
    %-- subplot inclination
    sp(2) = subplot(1,3,2);
    h = plot(X_data, y_data); 
    grid on
    set(gca,'ylim',y_limits)
    set(h(1), 'linewidth', 2)
    if size(X_data, 2) > 1
        set(h(2), 'linewidth', 2)
    end
    set(gca, 'ydir', 'rev')
    xlabel('[  ^{\circ} ]')
    if ~exist('Incl_X','var') && ~exist('Incl_Y','var') % No inclinometers
        legend('\theta_Y = sin^{-1}(Ax/g)', '\theta_X = sin^{-1}(Ay/g)',...
            'location','southeast')
    else % add on inclinometers
        LP_ADIS = 2; % low-pass filter applied to ADIS Inclinometer
        [b,a] =   butter(1,LP_ADIS/(fs_slow/2));% low-pass to show only the gravity signal
        X_data = [-Incl_Y Incl_X];
        X_data = filtfilt(b,a,X_data);
        hold on
        plot(X_data,Indep_var_slow) % Assumes inclinometers are sampled on slow channels
        legend('\theta_Y = sin^{-1}(Ax/g)', '\theta_X = sin^{-1}(Ay/g)',...
            ' -\theta_Y', '\theta_X','location','southeast')
    end
    title(title_string)
    
    %-- subplot dP/dt
    sp(3) = subplot(1,3,3);
    if ~isHMP
        h = plot(dP_dt_LP, Indep_var_slow ); 
        xlabel('d\itP\rm / d\itt \rm[ dBar s^{-1}]')
    else
        h = plot(speed_slow, Indep_var_slow ); 
        xlabel('estimated \itU \rm[ m s^{-1}]')
    end
    grid on
    set(gca,'ylim',y_limits)
    set(h(1), 'linewidth', 2)
    set(gca, 'ydir', 'rev')
    %set(gca,'yticklabel',[])

    linkaxes(sp,'y')
    drawnow

end

%-----------------------------------------------------------------
% ----- Process Data: Despike piezo-accelerometer sensors --------
%-----------------------------------------------------------------

%-- pre-assign despiking variables (even if no piezo-accelerometer)
spikes_A     = {};  
pass_count_A = zeros(size(AA_piezo,2),1);
fraction_A   = zeros(size(AA_piezo,2),1);

%-- despike the piezo-accelerometer signals
 piezo_accel_num = size(AA,2);
if  ~isempty(AA_piezo) && despike_A(1) ~= inf
   
    AA_raw = AA;
    for probe = 1:piezo_accel_num
        [AA(:,probe), spikes_A{probe}, ...
            pass_count_A(probe), fraction_A(probe)]  = ...
            despike(AA_raw(:,probe),  despike_A(1), ...
            despike_A(2), fs_fast, round(despike_A(3)*fs_fast));
    end
else 
    AA_raw = AA;
end

%------------------------------------------------------------------------
% ----- Figure: Vibration sensors ---------------------------------------
%------------------------------------------------------------------------
if ~isempty(AA_piezo) && make_figures && (plot_kinematics || plot_rawaccel)
    
    %-- low-pass filter parameters
    fc = 1; %[Hz]
    [b,a] =   butter(1,fc/(fs_fast/2));
        
    %-- plot
    fig_num = fig_num + 1;
    figure (fig_num); clf
    clear sp
    
    %-- number of subplots
    MM = 1; % default
    if plot_rawaccel == 1
        fig_aspectratio(gcf,2);
        MM = 2; 
    end
    
    %-- legend/title info
    legend_string = AA_names_piezo;
    for nn = 1:piezo_accel_num
        legend_string = [legend_string [AA_names_piezo{nn} ' LP']];
    end
    
    despike_title = [];
    if ~isempty(spikes_A)
        for nn = 1:piezo_accel_num
            despike_title = [despike_title ' ' AA_names_piezo{nn} ,...
                ': ' num2str(length(spikes_A{nn})) ' '];
        end
        despike_title = ['No. of Spikes Removed  = ' despike_title];
    else
        despike_title = ['No. of Spikes Removed  = Ax: 0 Ay: 0'];
    end
    
    if ~isHMP
        xlimits = [-1000 1000];
        if mean(mean(abs(AA)))>1000 % seaglider
            xlimits = xlimits+mean(mean((AA)));
        end
    else
        xlimits = max(max(AA))*[-1 1];
    end

    
    %-- subplot: raw acceleration
    if MM == 2
        sp(1) = subplot(1,MM,1);
        plot(AA_raw, Indep_var_fast)
        grid on
        hold all
        plot(filtfilt(b, a, AA_raw), Indep_var_fast,'linewidth',2);
        title([title_string ['No despiking']])
        xlim(xlimits)
    end
    
    %-- subplot: depiked acceleration 
    sp(MM) = subplot(1,MM,MM);
    plot(AA, Indep_var_fast,'linewidth',1)
    grid on
    hold all
    plot(filtfilt(b, a, AA), Indep_var_fast,'linewidth',2)
    
    title([title_string despike_title])
    
    if exist('y_limits', 'var')
        set(sp,'ylim',y_limits)
    end
    xlim(xlimits)
    set(sp, 'ydir', 'rev')
    xlabel(sp(1),'[counts]')
    xlabel(sp(end),'[counts]')
    ylabel(sp(1), Indep_labelstr)
    legend(sp(1),legend_string,'location','southeast');
end


%-----------------------------------------------------------------
% ----- Figure: Magnetometers ------------------------------------
%-----------------------------------------------------------------
if make_figures && plot_kinematics && (exist('Mx','var') || exist('My','var') || exist('Mz','var'))
    fig_num = fig_num + 1;
    figure (fig_num); clf
    plot_data = [];
    legend_string = [];
    if exist('Mx', 'var')
        plot_data = [plot_data Mx];
        legend_string{end+1} = 'M_x';
    end
    if exist('My', 'var')
        plot_data = [plot_data My];
        legend_string{end+1} = 'M_y';
    end
    if exist('Mz','var')
        plot_data = [plot_data Mz];
        legend_string{end+1} = 'M_z';
    end

    h = plot(plot_data, Indep_var_slow);
    set(gca, 'ydir', 'rev')
    xlabel('[\muT]')
    ylabel(Indep_labelstr)
    legend(legend_string, 'location', 'NorthEast');
    title(title_string)
    drawnow
end

%-----------------------------------------------------------------
% ----- Figure: Rotation Rate ------------------------------------
%-----------------------------------------------------------------
if make_figures && plot_kinematics && (exist('Rx','var') || exist('Ry','var') || exist('Rz','var'))
  fig_num = fig_num + 1;
    figure (fig_num); clf
    plot_data = [];
    legend_string = [];
    if exist('Rx', 'var')
        plot_data = [plot_data Rx];
        legend_string{end+1} =  'R_x';
    end
    if exist('Ry', 'var')
        plot_data = [plot_data Ry];
        legend_string{end+1} =  'R_y';
    end
    if exist('Rz', 'var')
        plot_data = [plot_data Rz];
        legend_string{end+1} =  'R_z';
    end

    h = plot(plot_data, Indep_var_slow);
    
    set(gca,'ydir', 'rev')
    xlabel('[\circ s^{-1}]')
    ylabel(Indep_labelstr)
    legend(legend_string, 'location', 'NorthEast');
    title(title_string)
    drawnow
end




%-----------------------------------------------------------------
% ----- Figure: Battery Voltage ----------------------------------
%-----------------------------------------------------------------
if exist('V_Bat','var') && plot_battery && make_figures
    fig_num = fig_num + 1;
    figure(fig_num); clf
    [b,a]=butter(1,1/(fs_slow/2));
    h = plot([V_Bat filtfilt(b,a,V_Bat)], Indep_var_slow);grid on
    set(h(2),'linewidth',3,'color','m')
    set(gca, 'ydir', 'rev')
    ylabel(Indep_labelstr)
    xlabel('\itV_{Bat}\rm [volts]')
    title(title_string)
    legend('V_{Bat}','V_{Bat}-LP','location','southeast')
    set(gca,'ylim',y_limits)
    drawnow
end


%-----------------------------------------------------------------
% ----- Process data: High accuracy CTD (if available)  ----------
%-----------------------------------------------------------------
% - time offsets (approximate)
offset = 0;
sbt_offset = 0.35;
JAC_offset = 0.14;

% - get variables
T_name = {'SBT1', 'SBT', 'sbt', 'JAC_T'}; % assume that at most only one is available
C_name = {'SBC1', 'SBC', 'sbc', 'JAC_C'};
all_variables = who('*','var'); % every variable
index_to_T = ismember(T_name, all_variables);
index_to_C = ismember(C_name, all_variables);
T_list = T_name(index_to_T);
C_list = C_name(index_to_C);

% - set flags
we_have_T = false;
we_have_C = false;
we_have_JAC_T = false;
we_have_JAC_C = false;
we_have_SBE_T = false;
we_have_SBE_C = false;

% - populate temperature vectors
T_CTD = [];
if ~isempty(T_list)
    we_have_T = true;
    eval(['T_CTD = ' T_list{1} ';'])
    junk = T_list{1};
    if strcmpi(junk(1),'s') % SeaBird
        offset = sbt_offset;
        we_have_SBE_T = true;
    elseif strcmpi(junk(1),'j') % JAC
        offset = JAC_offset;
        we_have_JAC_T = true;
    end
end

% - populate conductivity vectors, legend strings and offsets
C_CTD = [];
if ~isempty(C_list)
    we_have_C = true;
    eval(['C_CTD = ' C_list{1} ';'])
    junk = C_list{1};
    if strcmpi(junk(1),'s')
        offset = sbt_offset;
        we_have_SBE_C = true;
    elseif strcmpi(junk(1),'j')
        offset = JAC_offset;
        we_have_JAC_C = true;
    end
end

% Compute pressure offset (if vertical profiler)
if ~isHMP
    offset = -offset*mean(dP_dt_LP);
end

% - compute salinity and density
we_have_S = we_have_C & we_have_T;
if we_have_S
    if we_have_JAC_T 
        % Get mean speed 
        if isempty(constant_speed)
            mean_speed = mean(speed_slow);
        else
            mean_speed = constant_speed;
        end
        % compute salinity
        [S_CTD, C_CTD_match] = salinity_JAC (P_slow, T_CTD, C_CTD,...
                          'fs',fs_slow,'speed',mean_speed);
    else 
        S_CTD = salinity     (P_slow, T_CTD, C_CTD);
        C_CTD_match = C_CTD;
    end
    potential_T = theta(S_CTD, T_CTD, P_slow, 0); %pot_T relative to P=0;
    potential_density = sigma_p(potential_T, S_CTD, 0); % relative to P=0
end

%-----------------------------------------------------------------
% ----- Figure: High accuracy CTD (if available)  ----------------
%-----------------------------------------------------------------
if make_figures && plot_sensors 

    % - set up plot
    windows = 0;
    if we_have_T; windows = windows + 1; end
    if we_have_C; windows = windows + 1; end
    if we_have_S; windows = windows + 2; end
       
    % - plot
    if windows > 0
        fig_num = fig_num +1;
        figure (fig_num); clf
        fig_aspectratio(gcf,2);
        clear sp legend_string
        
        plot_data = [T_CTD C_CTD_match S_CTD potential_density];
        x_label = {'[\circC]', '[mS cm^{-1}]', '[PSU]', '[kg m^{-3}]'};
        legend_string = {'\itT_{\rm??}','\itC_{\rm??}','S','\sigma_0'}; % defaults
        
        for index = 1:windows
            sp(index) = subplot(1,windows,index);
            
            % - plot data
            plot(plot_data(:,index), Indep_var_slow + offset);grid on
            
            % - labels and limits
            if index == 1, ylabel(Indep_labelstr), end
            set(gca, 'ydir', 'rev')
            xlabel(x_label{index})
            
            % - legend
            if index == 1
                if we_have_JAC_T; legend_string{index} = '\itT_{\rmJAC}'; end
                if we_have_SBE_T; legend_string{index} = '\itT_{\rmSB}';  end 
            elseif index == 2
                if we_have_JAC_C; legend_string{index} = '\itC_{\rmJAC}'; end
                if we_have_SBE_C; legend_string{index} = '\itC_{\rmSB}';  end
            end

            legend(legend_string{index},'location','southwest')
            
        end
        title(sp(1),title_string)
        
        linkaxes(sp,'y')
        
    end
    drawnow
end

clear junk

%--------------------------------------------------------------------
% ----- Figure: Compare Thermistors and MicroC data with CTD data ---
%--------------------------------------------------------------------
Thermistor_name = 'T*_fast';
Thermistor_list = who(Thermistor_name,'var');
microconductivity_name = 'C*_fast';
microconductivity_list = who(microconductivity_name,'var');

T_count = length(Thermistor_list);
C_count = length(microconductivity_list);
CTD_count = length(T_CTD)+length(C_CTD);


T_data  = [];
T_names = {};
for kk = 1: length(Thermistor_list)
    T_data = [T_data eval(Thermistor_list{kk})];
    T_names = [T_names Thermistor_list{kk}(1:end-5)];
end
C_data  = [];
C_names = {};
for kk = 1: length(microconductivity_list)
    C_data = [C_data eval(microconductivity_list{kk})];
    C_names = [C_names microconductivity_list{kk}(1:end-5)];
end


if   make_figures && plot_sensors  &&  (T_count+C_count)>0 && CTD_count>0
    
    % - Themistor plot set up
    T_legend_string = {};
    T_plot_data     = [];
    T_slow_data      = [];
    for index = 1 : T_count
        T_legend_string = [T_legend_string texstr(Thermistor_list{index})];
        eval(['T_plot_data = [T_plot_data ' Thermistor_list{index} '];'])
    end
    
    % - CTD temperature plot set up    
    if we_have_T
        T_slow_data   = T_CTD;
        if we_have_JAC_T
            T_legend_string = [T_legend_string '\itT_{\rmJAC}'];
        elseif we_have_SBE_T
            T_legend_string = [T_legend_string '\itT_{\rmSB}'];
        end
    else
        Indep_data_slow = [];
    end
        
    % - MicroC plot set up
    C_legend_string = {};
    C_plot_data     = [];
    C_slow_data      = [];
    for index = 1 : C_count
        C_legend_string = [C_legend_string texstr(microconductivity_list{index})];
        eval(['C_plot_data = [C_plot_data ' microconductivity_list{index} '];'])
    end
    
    if we_have_C
        C_slow_data   = C_CTD_match;
        if we_have_JAC_C
            C_legend_string = [C_legend_string '\itC_{\rmJAC}'];
        elseif we_have_SBE_C
            C_legend_string = [C_legend_string '\itC_{\rmSB}'];
        end
    else
        Indep_data_slow = [];
    end
    
    % - plot
    fig_num = fig_num +1;
    figure (fig_num); clf
    clear sp
    if T_count>0 && C_count>0
        % have both thermistor and micro C
        sp(1) = subplot(1,2,1);
        h = plot(T_plot_data, Indep_var_fast, T_slow_data, Indep_var_slow + offset);grid on
        %set(h(end),'linewidth',2)
        xlabel('[  ^{\circ}C ]')
        legend(T_legend_string,  'location','SouthEast')

        sp(2) = subplot(1,2,2);
        h = plot(C_plot_data, Indep_var_fast, C_slow_data, Indep_var_slow + offset);grid on
        %set(h(end),'linewidth',2)
        xlabel('[ mS cm^{-1} ]')
        legend(C_legend_string,  'location','SouthEast')
    elseif T_count>0
        % only have thermistor
        sp(1) = subplot(1,1,1);
        h = plot(T_plot_data, Indep_var_fast, T_slow_data, Indep_var_slow + offset);grid on
        %set(h(end),'linewidth',2)
        xlabel('[  ^{\circ}C ]')
        legend(T_legend_string,  'location','SouthEast')
    elseif C_count>0
        % only have microC
        sp(1) = subplot(1,1,1);
        h = plot(plot_data, Indep_var_fast, slow_data, Indep_var_slow + offset);grid on
        %set(h(end),'linewidth',2)
        xlabel('[ mS cm^{-1} ]')
        legend(C_legend_string,  'location','SouthEast')
    end
    
    if ~exist('y_limits', 'var')
        y_limits = get(sp(1),'ylim');
    end
    set(sp,'ylim',y_limits)
    set(sp, 'ydir', 'rev')
    ylabel(sp(1),Indep_labelstr)
    
    title(sp(1),title_string)
    
    drawnow
end

%-----------------------------------------------------------------
% ----- Figure: Fluorometer and Backscatter Sensors --------------
%-----------------------------------------------------------------

% - Find variables
var_list = { ...
    'Fluo', ...
    'Turbidity', ...
    'Chlorophyll', ...
    'BS', ...
    'Suna'};
all_variables = who('*','var'); 
index_to_variables = ismember(var_list, all_variables);
var_list = var_list(index_to_variables);


if  make_figures  && plot_sensors && ~ isempty(var_list)
    fig_num = fig_num +1;
    figure (fig_num); clf
    
    legend_string = [];
    Y_fast        = [];
    Y_slow        = [];
    Indep_data_fast   = [];
    Indep_data_slow   = [];
    
    % - populate vectors and legend entries
    for index = 1 : length(var_list)
        eval (['junk = ' var_list{index} ';']);
        if length(junk) == length(t_fast)
            Y_fast = [Y_fast junk];
            Indep_data_fast = Indep_var_fast;
            legend_string{end+1} = var_list{index};
        end
        if length(junk) == length(t_slow)
            Y_slow = [Y_slow junk];
            Indep_data_slow = Indep_var_slow;
            legend_string{end+1} = var_list{index};
        end
    end
    
    % - plot
    h = plot(Y_fast, Indep_data_fast, Y_slow, Indep_data_slow);grid on
    
    % - limits and labels
    set(gca,'ylim',y_limits)
    set(gca, 'ydir', 'rev')
    ylabel(Indep_labelstr)
    xlabel('[ ppb ]  [ FTU ]')
    
    % - legend and title
    legend(legend_string)
    title(title_string)
    
    drawnow
end

clear junk

%-----------------------------------------------------------------
% ----- Process data: Shear, gradT and gradC ---------------------
%-----------------------------------------------------------------
% - Apply gentle filtering to signals 
%      (Band pass shear, low pass grad T and grad C)
% - Despike the shear and microconductivity signals
% - Calculate and assemble data into lists
% ----------------------------------------------------------------

% - Filter parameters
[bh,ah] = butter(1, HP_cut/(fs_fast/2), 'high');    % high pass shear data
[bl,al] = butter(4, LP_cut/(fs_fast/2));            % low pass all signals (for profile plotting only)

% - wish lists
all_variables = who('*', 'var');

shear_wish_list         = {'sh1'      , 'sh2'      , 'sh3'      , 'sh4'};
shear_legend_wish_list  = {'\nablau_1', '\nablau_2', '\nablau_3', '\nablau_4'};

gradT_wish_list         = {'gradT'    , 'gradT1'   , 'gradT2'};
gradT_legend_wish_list  = {'\nablaT'  , '\nablaT_1', '\nablaT_2'};

gradC_wish_list         = {'gradC'    , 'gradC1'   , 'gradC2'};
gradC_legend_wish_list  = {'\nablaC'  , '\nablaC_1', '\nablaC_2'};

%%%%%%%%%%%%%%
%    Shear   %
%%%%%%%%%%%%%%

% - get shear data
index_to_shear    = ismember(shear_wish_list, all_variables);
shear_list        = shear_wish_list(index_to_shear);
shear_legend_list = shear_legend_wish_list(index_to_shear);
shear_count = length(shear_list);
if shear_count > 0 
    SH = zeros(size(eval(shear_list{1}),1), shear_count); % initialize
else
    SH = [];
end
for index = 1:length(shear_list)
    SH(:, index) = eval(shear_list{index});
end

% - Check length of shear and accelerometer data
L_shear = size(SH,     1);
L_accel = size(AA,     1);
L_time  = size(t_fast, 1);
if (L_shear ~= L_accel) || (L_shear ~= L_time)
    warning ('shear and acceleration signals are not of equal length')
    goodman = false;
end

% - Despike shear data 
spikes_sh     = {}; % they must exist
pass_count_sh = zeros(shear_count,1);
fraction_sh   = zeros(shear_count,1);
if despike_sh(1) ~= inf
    for index = 1:shear_count
        [SH(:,index), spikes_sh{index}, pass_count_sh(index), fraction_sh(index) ] =  ...
            despike(SH(:,index), despike_sh(1), despike_sh(2), fs_fast, round(despike_sh(3)*fs_fast));
    end
end

% - High-pass shear data
%----
% NOTE: 
%   SH_HP = filtfilt(bh, ah, SH);
% does not work well for high variance nearthe bottom of the profile.
%-----
SH_HP = filter(bh, ah, SH);
SH_HP = flipud(SH_HP);
SH_HP = filter(bh, ah, SH_HP);
SH_HP = flipud(SH_HP);

% - Low-pass the shear data (for display purposes)
SH_BP = filtfilt(bl, al, SH_HP); 

%%%%%%%%%%%%%%
%    gradT   %
%%%%%%%%%%%%%%

% - set scales
T_scale        = 1;
T_scale_string = '';
if tidal_channel
    T_scale        = 10;
    T_scale_string = '10\times';
end

% - get appropriate gradient data
index_to_gradT    = ismember(gradT_wish_list, all_variables);
gradT_list        = gradT_wish_list(index_to_gradT);
gradT_legend_list = gradT_legend_wish_list(index_to_gradT);
T_count           = length(gradT_list);
if T_count > 0
    GRADT = zeros(size(eval(gradT_list{1}), 1), T_count);
end
for index = 1:length(gradT_list)
	GRADT(:, index) = T_scale*eval(gradT_list{index});
end

% - Low-pass filter grad T
GRADT_LP = [];
if T_count>0, GRADT_LP = filtfilt(bl, al, GRADT);end

%%%%%%%%%%%%%%
%    gradC   %
%%%%%%%%%%%%%%

% - set scales
C_scale        = 1;
C_scale_string = '';
if tidal_channel
    C_scale        = 1;
    C_scale_string = '1\times';
end

% - get appropriate data
index_to_gradC    = ismember(gradC_wish_list, all_variables);
gradC_list        = gradC_wish_list(index_to_gradC);
gradC_legend_list = gradC_legend_wish_list(index_to_gradC);
C_count           = length(gradC_list);
if C_count > 0
    GRADC = zeros(size(eval(gradC_list{1}),1), C_count);
end

for index = 1:length(gradC_list)
    GRADC(:,index)            = C_scale*eval(gradC_list{index});
end

% - Despike
spikes_C     = {};                    
pass_count_C = zeros(C_count,1);
fraction_C   = zeros(C_count,1);
if despike_C(1) ~= inf
    for index = 1:C_count
        [GRADC(:,index), spikes_C{index}, pass_count_C(index), fraction_C(index)] =  ...
            despike(GRADC(:,index), despike_C(1), despike_C(2), fs_fast, round(despike_C(3)*fs_fast));
    end
end

% - Low pass filter
GRADC_LP = [];
if C_count>0, GRADC_LP = filtfilt(bl, al, GRADC); end

%-----------------------------------------------------------------
% ----- Figure: Shear, gradT and gradC ---------------------------
%-----------------------------------------------------------------

if make_figures && plot_sensors && (shear_count + T_count + C_count)>0
    
    fig_num = fig_num + 1;
    figure(fig_num); clf
    fig_aspectratio(gcf,2.5);
%     set(gcf,'Position',[100 100 1000 500])
    clear sp
    
    % - variables to plot
    plot_gradients_1 = {'SH_HP','SH_BP'};
    plot_gradients_2 = {'GRADT_LP','GRADC_LP'};
    
    % - subplots: 1) shear, 2) scalar gradients
    for ii = 1:2
        X_gradient_data  = [];
        X_gradient_names = {};
        
        if ii == 1
            plot_gradients = plot_gradients_1;
        elseif ii == 2
            plot_gradients = plot_gradients_2;
        end
        
        for kk = 1:length(plot_gradients)
            % - Data
            X_gradient_data = [X_gradient_data eval(plot_gradients{kk})];
            
            % - Names
            if strcmp(plot_gradients{kk},'SH_HP') || strcmp(plot_gradients{kk},'SH_BP')
                legend_list = strcat(shear_legend_list,{' '},plot_gradients{kk}(end-1:end));
            elseif strcmp(plot_gradients{kk},'GRADT_LP')
                legend_list = strcat(T_scale_string,{' '}, gradT_legend_list,{' '},plot_gradients{kk}(end-1:end));
            elseif strcmp(plot_gradients{kk},'GRADC_LP')
                legend_list = strcat(C_scale_string,{' '}, gradC_legend_list,{' '},plot_gradients{kk}(end-1:end));
            end
            X_gradient_names = [X_gradient_names legend_list];
            
        end
        
        % - Set offsets
        x_gradient_step   = 5;
        if tidal_channel, x_gradient_step = 20; end
        left_edge = -x_gradient_step;
        x_gradient_offset = left_edge + x_gradient_step;
        
        LL = size(X_gradient_data,2);
        for k = 1:LL
            X_gradient_data(:,k) = X_gradient_data(:,k) + x_gradient_offset;
            x_gradient_offset = x_gradient_offset + x_gradient_step;
        end
        
        % - Set up axes
        if ii == 1
            pos = [0.08, 0.11, 0.4, 0.75];
            sp(ii) = subplot('Position',pos); 
        else
            pos  = [0.5, 0.11, 0.2, 0.75];
            sp(ii) = subplot('Position',[0.5, 0.11, 0.2, 0.75]);
        end
        
        % - plot
        if ~isempty(X_gradient_data) % Might be empty for MicroPods
            plot(X_gradient_data, Indep_var_fast );
        end
        grid on
        
        % - limits and labels
        set(gca,'ylim',y_limits)
        set(gca,'ydir','rev')
        x_lim = get(gca, 'xlim');
        x_lim = [-x_gradient_step-1 LL*x_gradient_step+1];%x_lim(2)];
        if tidal_channel, x_lim = [-20 x_gradient_offset]; end
        set(gca,'xlim',x_lim)
        if ii == 1; xlabel('[s^{-1}]'); end
        if ii == 2; xlabel('[\circC m^{-1}]'); end
        
        % - legend
        legend(X_gradient_names, 'position',[pos(1)+pos(3)-0.07 .68 0.08 0.1]);
    end
    uistack(sp(1),'bottom')
    uistack(sp(2),'bottom')
    
    % - subplots (on the same y-axes): 3) thermistor, 4) microC 
    sp(3) = subplot('Position',[0.72, 0.11, 0.2, 0.75]);
    if ~isempty(C_data)
        sp(4)=axes('Position',get(sp(3),'Position'));
    end
        
    % - plot T
    if ~isempty(T_data)
        plot(sp(3),T_data,Indep_var_fast,'Parent',sp(3));
    end
    
    % - limits and labels
    set(sp(3),'YDir','rev','Ylim',y_limits)
    set(get(sp(3),'xlabel'),'string','T [\circ C]')

    % - legends
    lgnd0=legend(sp(3),T_names,'location','northeast');
    lgnd0_pos = get(lgnd0,'Position');
    
    % - plot C
    if ~isempty(C_data)
        col = ColOrd(T_count+1,:); %assumes there is only one conductivity sensor
        plot(sp(4),C_data,Indep_var_fast,'color',col);
        set(sp(4),'XAxisLocation','top',...
            'YAxisLocation','right','Color','none',...
            'XColor',col,'YColor','none',...
            'YDir','rev','Ylim',y_limits)
        % - limits
        x_lim = get(sp(4),'Xlim');
        if diff(x_lim)<0.1
            set(gca,'Xlim',[x_lim(1)-1 x_lim(2)+1])
        end
        set(sp(4),'box','off')
        grid off
        
        % - labels
        set(get(sp(4),'xlabel'),'string','[mS cm^{-1}]')
        
        % - legend
        lgnd=legend(C_names,'location','east');
%         lgnd_pos = get(lgnd,'Position');
%         set(lgnd,'Position',[lgnd0_pos(1) lgnd0_pos(2)+lgnd0_pos(4)-0.04 lgnd0_pos(3) lgnd_pos(4)])
         uistack(sp(3),'bottom')
%         uistack(sp(4),'bottom')


    end
    linkaxes(sp,'y')
    ylabel(sp(1),Indep_labelstr)
    set(sp(2),'Yticklabels',[])
    set(sp(3),'YAxisLocation','Right')
  
    
    if tidal_channel
        junk = ['HP = ' num2str(HP_cut) ' Hz'];
    else
        junk = ['HP = ' num2str(HP_cut) ' Hz, BP = ' num2str(HP_cut) ' \cdot\cdot\cdot ' num2str(LP_cut) ' Hz'];
    end
    title(sp(1),[title_string junk]) 
    drawnow
end

clear junk

%-----------------------------------------------------------------
% ----- Process data: Frequency and wavenumber spectra -----------
%-----------------------------------------------------------------

% - Get indices for selected range
range = [];
if ~isempty(P_start) && ~isempty(P_end) % will actually be t_start and t_end for HMP
    range = find((Indep_var_fast > P_start) & (Indep_var_fast <= P_end));
    if isempty(range)
        if ~isHMP
            warning(['The pressure range of ' num2str(P_start) ' to ' ...
                num2str(P_end) ' was not found in this profile'])
        else
            warning(['The time range of ' num2str(P_start) ' to ' ...
                num2str(P_end) ' was not found in this file'])
        end
    end
end

if ~isempty(range) && ~isempty(SH)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      Shear spectra      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % - extract range 
    y_start = round(Indep_var_fast(range(1))); 
    y_end   = round(Indep_var_fast(range(end)));
    
    % - set up info struct
    info.fft_length    = fft_num;
    info.diss_length   = length(range);
    info.overlap       = info.diss_length/2;
    info.fs_fast       = fs_fast;
    info.fs_slow       = fs_slow;
    info.speed         =       speed_fast(range);
    info.T             = temperature_fast(range);
    info.t             =           t_fast(range);
    info.P             =           P_fast(range);
    info.Data_fast     = [];
    info.Data_slow     = [];
    info.fit_order     = fit_order;
    info.f_AA          = f_AA;
    info.f_limit       = f_limit;
    info.fit_2_isr     = fit_2_isr;
    info.goodman       = goodman;
    
    % - Compute spectra and dissipation rate for selected range
    if exist('EMC_Cur','var') && (length(EMC_Cur) == size(AA,1)) % if there is an EMC
        diss = get_diss_odas(SH_HP(range,:), [AA(range,:) EMC_Cur(range)], info);
    else
        diss = get_diss_odas(SH_HP(range,:),  AA(range,:),                 info);
    end
    
    % - define simple variables (use Pk to denote wavenumber spectra)
    mean_speed = diss.speed;
    F = diss.F;
    K = diss.K;
    for index = 1:shear_count
        r = num2str(index); % put it into a string
        eval(['Pk_sh' r       ' = diss.sh(:,'       r ',' r ');'])
        eval(['Pk_sh' r '_clean = diss.sh_clean(:,' r ',' r ');'])
    end
        
    % - Convert to frequency space (use P to denote frequency spectra)
    for index = 1:shear_count
        r = num2str(index); % put it into a string
        junk = ['P_' shear_list{index}]; % The name of the spectrum
        eval([junk '       = diss.sh      (:,' r ',' r ') / mean_speed;' ])
        eval([junk '_clean = diss.sh_clean(:,' r ',' r ') / mean_speed;' ])
    end
   
    clear junk 
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %     scalar spectra   % 
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    % - initialize and setup parameters
    scalar_wish_list = {...
        'T_dT', 'T1_dT1', 'T2_dT2', 'C_dC', 'C1_dC1', 'C2_dC2'};
    scalar_vectors_wish_list = {...
        'gradT', 'gradT1', 'gradT2', 'gradC', 'gradC1', 'gradC2'};
    all_variables = who('*','var');
    index_to_scalar_list        = ismember(scalar_wish_list,         all_variables);
    index_to_scalar_vector_list = ismember(scalar_vectors_wish_list, all_variables);
    scalar_vector_names = scalar_vectors_wish_list(index_to_scalar_vector_list);
    scalar_names        = scalar_wish_list        (index_to_scalar_list);
    if index_to_scalar_list ~= index_to_scalar_vector_list
        warning('names of scalar-gradients and scalar-channels are not matched!')
    end
    L = 0;
    if ~isempty(scalar_vector_names)
        L = eval(['length(' scalar_vector_names{1} ');']);
    end
    
    % - identify method used to compute gradient
    if isfield(params, 'gradT_method') && isfield(params, 'gradC_method')
        gradient_method = params.gradT_method;
    else
        gradient_method = 'first_difference'; % for older mat-files
    end
    
    % - get diff_gain from setup string and put scalar data into matrix
    diff_gain      = zeros(length(scalar_names),1);
    scalar_vectors = zeros(L, length(scalar_vector_names));
    for index = 1:length(scalar_names) 
        scalar_vectors(:, index) = eval([scalar_vector_names{index} ';']);
        tmp = setupstr(cfgobj, scalar_names{index}, 'diff_gain');  % Find diff_gain
        diff_gain (index) = str2double(tmp{1});% Assume it is numeric
    end
    
    % - Set up info struct
    scalar_info.fft_length      = fft_num;
    scalar_info.spec_length     = length(range);
    scalar_info.overlap         = info.overlap;
    scalar_info.fs              = fs_fast;
    scalar_info.gradient_method = gradient_method;
    scalar_info.diff_gain       = diff_gain;
    scalar_info.f_AA            = f_AA;
    if exist('EMC_Cur','var') && (length(EMC_Cur) == size(scalar_vectors,1))
        % Check if we have a JAC electromagnetic current meter reference signal.
        scalar_info.goodman         = goodman;
    else
        scalar_info.goodman         = false;
    end
    
    % - Compute scalar spectra
    if ~isempty(scalar_vector_names)
        if scalar_info.goodman
            scalar_spectra = get_scalar_spectra_odas(...
                scalar_vectors(range,:), EMC_Cur(range), ...
                info.P, info.t, info.speed, scalar_info);
        else
            scalar_spectra = get_scalar_spectra_odas(...
                scalar_vectors(range,:), [], ...
                info.P, info.t, info.speed, scalar_info);
        end
    end
    
    % - Define simple variables (i.e. Pk for wavenumber spectra)
    if ~isempty(scalar_vector_names)
        for index = 1:length(scalar_vector_names)
            name = ['Pk_' scalar_vector_names{index}];
            spectra = scalar_spectra.scalar_spec(:,index);
            eval([name ' = spectra;']);
        end
    end
            
    % - Convert to frequency space
    for index = 1:length(scalar_vector_names)
        % extract the scalar spectra
        name = ['P_' scalar_vector_names{index}];
        spectra = scalar_spectra.scalar_spec(:,index) / scalar_spectra.speed;
        eval([name ' = spectra;']);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     correct gradT spectra    %
    %       (only wavenumber)      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % - frequency response after Vachon and Lueck (double pole)
    F_0 = 25*sqrt(mean_speed); 
    tau_therm = 2*pi*F_0 / sqrt(sqrt(2) - 1); 
    tau_therm = 1 / tau_therm;
    Hinv = (1 + (2*pi*tau_therm*F).^2).^2; % inverse of the frequency response
    
    % - correction (Hinv is nondimensional so can apply directly to Pk_gradT)
    if exist('Pk_gradT','var')  Pk_gradT  = Hinv.*Pk_gradT;  end
    if exist('Pk_gradT1','var') Pk_gradT1 = Hinv.*Pk_gradT1; end
    if exist('Pk_gradT2','var') Pk_gradT2 = Hinv.*Pk_gradT2; end

end

%-----------------------------------------------------------------
% ----- Figure: Frequency and wavenumber spectra -----------------
%-----------------------------------------------------------------

if ~isempty(range) && make_figures && plot_spectra
    
    % - Set up figure and axes
    fig_num = fig_num + 1;
    figure(fig_num); clf; clear sp
    fig_aspectratio(gcf,2.5);
    sp(1) = axes('Position',[0.06 0.12 0.45 .7]); % wavenumber spectra
    sp(2) = axes('Position',[0.58 0.12 0.4 .7]);  % frequency spectra
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  - Plot Wavenumber spectra - %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axes(sp(1))
    
    % - extract K_max and dissipation rates
    K_max = diss.K_max;
    nu = diss.nu;
    e = diss.e;
    e1_string = []; e2_string = []; e3_string = []; e4_string = []; % Initializing
    for index = 1:shear_count
        r = num2str(index); % Put it into a string
        eval(['e' r ' = e(' r ');'])
        eval(['K_max_index_' r ' = find (K == K_max(' r '));'])
        eval(['phi' r   ' = diss.Nasmyth_spec(:,' r ');'])
        e_value = eval(['e' r ';']);
        junk = ['\epsilon=' make_scientific(e_value,2) 'W kg^{-1}'];
        eval(['e' r '_string =  junk ;'])
    end
    
    clear junk 
    % - Settings for all channels
    plot_wish_list = {
        % Name          Legend Title         Colour           Size
        {'Pk_sh1_clean', '\nablau_1 clean',   ColOrd(1,:),    3   }
        {'Pk_sh1',       '\nablau_1',         ColOrd(1,:),    1.5 }
        {'phi1',        e1_string,           'k',             1.5 }
        {'Pk_sh2_clean', '\nablau_2 clean',   ColOrd(2,:),    3   }
        {'Pk_sh2',       '\nablau_2',         ColOrd(2,:),    1.5 }
        {'phi2',        e2_string,           'k',             1.5 }
        {'Pk_sh3_clean', '\nablau_3 clean',   gold,           3   }
        {'Pk_sh3',       '\nablau_3',         gold,           1.5 }
        {'phi3',        e3_string,           'k',             1.5 }
        {'Pk_sh4_clean', '\nablau_4 clean',   fuj,            3   }
        {'Pk_sh4',       '\nablau_4',         fuj,            1.5 }
        {'phi4',        e4_string,           'k',             1.5 }
        {'Pk_gradT',     '\nablaT',           ColOrd(3,:),    2   }
        {'Pk_gradT1',    '\nablaT_1',         ColOrd(3,:),    2   }
        {'Pk_gradT2',    '\nablaT_2',         ColOrd(5,:),    2   }
        {'Pk_gradC1',    '\nablaC_1',         ColOrd(4,:),    2   }
        {'Pk_gradC2',    '\nablaC_2',         ColOrd(6,:),    2   }};
    
    point_wish_list = {
        % X pos         Y pos                             Label                                                 Colour          Size
        {'K_max(1)',    'Pk_sh1_clean(K_max_index_1)',   '[''k_{max} u_1='' num2str(round(K_max(1))) ''cpm'']', ColOrd(1,:),    18  }
        {'K_max(2)',    'Pk_sh2_clean(K_max_index_2)',   '[''k_{max} u_2='' num2str(round(K_max(2))) ''cpm'']', ColOrd(2,:),    18  }
        {'K_max(3)',    'Pk_sh3_clean(K_max_index_3)',   '[''k_{max} u_3='' num2str(round(K_max(3))) ''cpm'']', gold,           18  }
        {'K_max(4)',    'Pk_sh4_clean(K_max_index_4)',   '[''k_{max} u_4='' num2str(round(K_max(4))) ''cpm'']', fuj,            18  }};
    
    % - Remove all non-applicable channels
    plot_list = {};
    for ch = plot_wish_list'
        if exist(ch{1}{1}, 'var'), plot_list{end+1,1} = ch{1}; end
    end
    
    % - Remove all non-applicable points
    point_list = {};
    for pt = point_wish_list'
        try
            point_list{end+1,1} = {eval(pt{1}{1}), eval(pt{1}{2}), eval(pt{1}{3}), pt{1}{4}, pt{1}{5}};
        catch, continue; end
    end
    
    % - Extract the required data
    Y_data = [];
    for ch = plot_list', Y_data(:,end+1) = eval( ch{1}{1} ); end
    
    % - Plot the spectral lines and set parameters
    log_plot = loglog(K, Y_data);
    for ch = 1:size(plot_list,1)
        set(log_plot(ch),'linewidth',plot_list{ch}{4},'color',plot_list{ch}{3});
    end
    
    % - Plot the requested points with correct colours.
    hold on;
    for pt = point_list'
        loglog( pt{1}{1}, pt{1}{2},         ...
            'Color',            'w',         ...
            'Marker',           '^',        ...
            'MarkerSize',       pt{1}{5},   ...
            'MarkerFaceColor',  pt{1}{4},   ...
            'MarkerEdgeColor',  'w');
    end
    hold off
    grid on;
    
    % - limits and labels
    set(gca, 'ylim',[1e-8 1e0], 'xlim', [K(2) K(end)])
    if tidal_channel, set(gca,'ylim',[1e-7 1e1]), end
    xlabel('\itk \rm [cpm]')
    ylabel('\Phi (\itk\rm)  [s^{-2} cpm^{-1}]')
    
    % - legend
    legend_list = {};
    for plt = plot_list',  legend_list{end+1} = plt{1}{2}; end
    for plt = point_list', legend_list{end+1} = plt{1}{3}; end
    legend(legend_list,'location','NorthEastOutside');
    
    % - title
    title_string_2_old = title_string{2};
    title_string{2} = [title_string{2} ',   Method = [' num2str(diss.method'),' ]' ];
    if ~isHMP
        range_string = sprintf('%d < P < %d m',min(y_start,y_end), max(y_start,y_end));
    else
        range_string = sprintf('%d < t < %d s,',min(y_start,y_end), max(y_start,y_end));
    end
    new_title_string = [title_string , strcat(range_string,...
            sprintf(' speed = %0.2f m s^{-1}, f_{HP} = %0.2f Hz',mean_speed, HP_cut))];
    title(new_title_string)
    title_string{2} = title_string_2_old;
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  - Plot Frequency spectra -  % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axes(sp(2))
    
    % - Scale down acceleration spectra to similar scale as shear
    if piezo
        diss.AA = 1e-6*diss.AA;
    end
    
    % - Remove EM current meter from acceleration matrix
    accel_count = size(diss.AA,2);
    if exist('EMC_Cur','var'), accel_count = accel_count - 1; end % In case we have an EM current meter, accel_count is too big by 1
    
    % - Define simple acceleration variables
    if accel_count > 0, P_Ax = diss.AA(:,1,1);end
    if accel_count > 1, P_Ay = diss.AA(:,2,2);end
    if accel_count > 2, P_Az = diss.AA(:,3,3);end
    
    % - Get Nasmyth curves for several dissipation rates
    epsilon = 1e-10*[1 10 100 1e3 1e4 1e5];
    if tidal_channel, epsilon = 10*epsilon; end
    nu = visc35(mean(temperature_fast(range)));
    [Pn, kn]=nasmyth(epsilon, nu);
    fn = kn*mean_speed; Pn = Pn/mean_speed;
    
    % - Line settings for all channels.
    plot_wish_list = {
        % Name          Legend Title         Colour         Size
        {'P_Ax',        'A_x',               cyan,          3 }
        {'P_Ay',        'A_y',               green,         3 }
        {'P_Az',        'A_z',               brown,         3 }
        {'P_sh1',       '\nablau_1',         ColOrd(1,:),   1.5 }
        {'P_sh2',       '\nablau_2',         ColOrd(2,:),   1.5 }
        {'P_sh3',       '\nablau_3',         gold,          1.5 }
        {'P_sh4',       '\nablau_4',         fuj,           1.5 }
        {'P_gradT',     '\nablaT',           ColOrd(3,:),   2   }
        {'P_gradT1',    '\nablaT_1',         ColOrd(3,:),   2   }
        {'P_gradT2',    '\nablaT_2',         ColOrd(5,:),   2   }
        {'P_gradC1',    '\nablaC_1',         ColOrd(4,:),   2   }
        {'P_gradC2',    '\nablaC_2',         ColOrd(6,:),   2   }};
    
    % - Remove all non-applicable channels
    plot_list = {};
    for ch = plot_wish_list'
        if exist(ch{1}{1}, 'var'), plot_list{end+1,1} = ch{1}; end
    end
    
    % - Extract the required data
    Y_data = [];
    for ch = plot_list', Y_data(:,end+1) = eval( ch{1}{1} ); end
    
    % - Plot the spectral lines and set parameters
    log_plot = loglog(F, Y_data);
    for ch = 1:size(plot_list,1)
        set(log_plot(ch),'linewidth',plot_list{ch}{4},'color',plot_list{ch}{3});
    end
    
    % - Plot Naysmyth curves
    hold all
    loglog(fn,Pn,'k')
    
    % - limits and labels
    f_upper = max(F);
    y_limit = [1e-9 10]; % limits for y-axis
    x_limit = [1./fft_length f_upper];% limits for x-axis
    set(gca, 'ylim', y_limit, 'xlim',x_limit)
    ylabel('[Variance Hz^{-1}]')
    xlabel('\it f \rm [Hz]')
     
    % - labels for Nasmyth curves    
    for index = 1: length (epsilon)
        f_location = find(fn(:,index) > 2*x_limit(1));
        if isempty(f_location)
            f_location = 1;
        else
            f_location = f_location(1);
        end
        h=text(fn(f_location,index), Pn(f_location,index)*1.3, ...
            num2str(log10(epsilon(index)),2), ...
            'fontsize', 20, 'HorizontalAlignment', 'right');
    end    
    
    
     % - legend
    legend_list = {};
    for plt = plot_list',  legend_list{end+1} = plt{1}{2}; end
    legend_list{end+1} = 'Nasmyth';
    legend(legend_list,'location','NorthEastOutside');
    
    % - Title
    new_title_string = title_string;
    new_title_string = [title_string , strcat(range_string,...
            sprintf('  f_{HP} = %0.2f Hz', HP_cut))];
    title(new_title_string)
    drawnow
end

%-------------------------------------------------------------------------
% ----- Process Data: Dissipation profile and 'diss' structure -----------
%--------------------------------------------------------------------------
%  Move vectors into matrices that will be passed to the get_diss_odas 
%  function, where they will be averaged over the range of each dissipation 
%  estimate. All slow channels go into matrix 'Data_slow', and the fast 
%  channels go into 'Data_fast'. The names of the vectors go ito the cell 
%  arrays 'fast_list' and 'slow_list'.
% -------------------------------------------------------------------------

% - Clear matrices (so that they don't get counted in the search of variables)
Data_fast = []; 
Data_slow = [];

% - All variables
S = whos; 

% - Get sizes and indices
junk = zeros(size(S,1),2);
for k = 1:size(S,1)
    junk(k,:) = S(k).size; % Get the size of every variable
end
index_fast = find((junk(:,1) == size(t_fast,1)) & junk(:,2) == 1); % index_fast points to the fast channels.
index_slow = find((junk(:,1) == size(t_slow,1)) & junk(:,2) == 1); % index_slow points to the slow channels.

% - Populate Data_fast matrix
fast_list = [];
Data_fast = zeros(size(t_fast,1), size(index_fast,1)); % prefill
for k = 1:size(index_fast,1),
    fast_list{k} = S(index_fast(k)).name;
    Data_fast(:,k) = eval( S(index_fast(k)).name );
end

% - Populate Data_slow matrix
slow_list = [];
Data_slow = zeros(size(t_slow,1), size(index_slow,1));
for k = 1:size(index_slow,1)
    slow_list{k} = S(index_slow(k)).name;
    Data_slow(:,k) = eval( S(index_slow(k)).name );
end

% - Set up info struct
info.fft_length    = fft_num;
info.diss_length   = round(diss_length*fs_fast);
info.overlap       = round(overlap    *fs_fast);
info.fs_fast       = fs_fast;
info.fs_slow       = fs_slow;
info.speed         = speed_fast;
info.T             = temperature_fast;
info.t             = t_fast;
info.P             = P_fast;
info.fit_2_isr     = fit_2_isr;
info.f_AA          = f_AA;
info.f_limit       = f_limit;
info.Data_fast     = Data_fast;
info.Data_slow     = Data_slow;
info.fast_list     = fast_list;
info.slow_list     = slow_list;
info.goodman       = goodman;

% - Add EMC_Cur to acceleration matrix (to be used for noise removal)
if exist('EMC_Cur','var') && (length(EMC_Cur) == size(AA,1))
    AAA = [ AA EMC_Cur];
else
    AAA = AA;
end

% - **** Compute dissipation rates ******
if ~isempty(SH_HP)
    diss = get_diss_odas(SH_HP, AAA, info);
end

%-----------------------------------------------------------------
% ----- Figure: Dissipation profile ------------------------------
%-----------------------------------------------------------------
if shear_count >0  && make_figures && plot_dissipation
    fig_num = fig_num + 1;
    figure(fig_num); clf
    
    e1 = diss.e(1,:)';
    if shear_count > 1
        e2 = diss.e(2,:)';
        ratio_limit = 2.5;
        e_ratio = e1 ./ e2; % dissipation ratio
        bad_ratio = find(e_ratio > ratio_limit);
        if ~isempty(bad_ratio),e1(bad_ratio) = nan(size(bad_ratio));end% e1 is too large
        bad_ratio = find(e_ratio < 1/ratio_limit);
        if ~isempty(bad_ratio),e2(bad_ratio) = nan(size(bad_ratio));end% e2 is too large
    end
    if ~isHMP
        y = diss.P;
    else
        y = diss.t;
    end
    
    h=semilogx(diss.e', y, '-o','linewidth', 2); hold on


    hold off
    grid on
    legend_string = [];
    for ii = 1:length(h)
        col = get(h(ii),'color');
        set(h(ii), 'markersize', 8, 'markerfacecolor', col)
        legend_string{ii} = ['\epsilon_',num2str(ii)];
    end
    set(gca,'ylim',y_limits)
    set(gca,'ydir','rev') 
%     set(gca,'view',[0,-90]) % For R2018b, this prevents labels from overlapping axes

    
    x_limits = get(gca,'xlim');
    if log10(x_limits(2) / x_limits(1)) < 1
        x_limits(1) = 10^(floor(log10(x_limits(1))));
        x_limits(2) = 10^( ceil(log10(x_limits(2))));
        set(gca,'xlim',x_limits)
    end
    
    
    ylabel (Indep_labelstr)
    xlabel ('\epsilon [W kg^{-1}]')

    
    legend(legend_string,'location','NorthEastOutside');
    
    title(title_string)
    drawnow
end

%-----------------------------------------------------------------
% ----- Process Data: Compute scalar spectra (if any) ------------
%-----------------------------------------------------------------

% - identify scalars
scalar_wish_list    = {'gradT', 'gradT1', 'gradT2', 'gradC', 'gradC1', 'gradC2'};
name_wish_list      = {'T_dT' , 'T1_dT1', 'T2_dT2', 'C_dC' , 'C1_dC1', 'C2_dC2'};
all_variables       = who('*', 'var');
index_to_list       = ismember(scalar_wish_list, all_variables);
scalar_vector_names = scalar_wish_list(index_to_list);
name_list           =   name_wish_list(index_to_list);
we_have_scalars     = false;

% - Populate matrices of scalar data
L = 0;
if ~isempty(scalar_vector_names)
    we_have_scalars = true;
    eval(['L  = length(' scalar_vector_names{1} ');']);
end
scalar_vectors  = zeros(L, length(scalar_vector_names)); % initialize
for index = 1:length(scalar_vector_names)
    scalar_vectors(:, index) = eval(scalar_vector_names{index}) ; 
end

% - Get differentiator gains
obj = setupstr(setupfilestr); 
diff_gain       = zeros(   length(scalar_vector_names),1);

for index = 1:length(scalar_vector_names)
    tmp = setupstr(obj, name_list{index}, 'diff_gain'); % Find diff_gain
    diff_gain(index) = str2double(tmp{1});              % Assume it is numeric
end

% - assign parameters to scalar_info (same as those for dissipation rate)
scalar_info.fft_length      = fft_num;
scalar_info.spec_length     = info.diss_length;
scalar_info.overlap         = info.overlap;
scalar_info.fs              = fs_fast;
scalar_info.diff_gain       = diff_gain;
scalar_info.f_AA            = f_AA;
if isfield(params, 'gradT_method') && isfield(params, 'gradC_method')
    gradient_method = params.gradT_method;
else
    gradient_method = 'first_difference'; % for older mat-files
end
scalar_info.gradient_method = gradient_method;
if exist('EMC_Cur','var') && (length(EMC_Cur) == size(scalar_vectors,1))
    % Check if we have a JAC electromagnetic current meter reference signal.
    scalar_info.goodman  = goodman;
else
    scalar_info.goodman  = false;
    EMC_Cur               = []; % empty so that the call has only one format.
end

% - Compute spectra -
scalar_spectra = []; % in case there are no scalars for spectra
if we_have_scalars
    scalar_spectra = get_scalar_spectra_odas(...
        scalar_vectors, EMC_Cur, info.P, info.t, info.speed, scalar_info);
end

%-----------------------------------------------------------------
% ----- Populate the diss structure ------------------------------
%-----------------------------------------------------------------

% - Add fields to returned 'diss' structure
diss.piezo         = piezo;

diss.fast_list     = fast_list; 
diss.slow_list     = slow_list;

diss.spikes_A      =     spikes_A; % The spikes found in piezo-accelerometer signals
diss.pass_count_A  = pass_count_A; % The number of passes of the despike function
diss.fraction_A    =   fraction_A; % The fraction of accelerometer data removed by the despike function

diss.spikes_sh     =     spikes_sh; % spikes found in the shear signals
diss.pass_count_sh = pass_count_sh;
diss.fraction_sh   =   fraction_sh;

diss.spikes_C      =     spikes_C; % The spikes found in micro_C signals
diss.pass_count_C  = pass_count_C;
diss.fraction_C    =   fraction_C;

diss.profiles_in_this_file = profiles_in_this_file;
diss.speed_source          = speed_source;
diss.params                = params;
diss.vehicle_info          = vehicle_info;

diss.scalar_spectra     = scalar_spectra; 
diss.scalar_vector_list = scalar_vector_names;
diss.scalar_info        = scalar_info;

diss.ql_info    = ql_info;
diss.ql_info_in = ql_info_in;

% - Set return value
result = diss;

%-----------------------------------------------------------------
% ----- Figure: wavenumber spectrogram of shear ------------------
%-----------------------------------------------------------------

% Assumes that sh1 and sh2 exist and does not handle the case of more
% shear probes 

if  make_figures && plot_spectrograms
    if ~isHMP
        y = diss.P;
    else
        y = diss.t;
    end
    
    fig_num = fig_num + 1;
    figure(fig_num); clf
    
    sh1_spec = squeeze(diss.sh_clean(:,1,1,:));
    if shear_count>1, sh2_spec = squeeze(diss.sh_clean(:,2,2,:)); end
    
    subplot(1,shear_count,1)
    pcolor(diss.K', y, log10(sh1_spec'));grid on
    %set(gca,'XScale','log')
    set(gca,'xlim',[1 150])
    %set(gca,'zlim',[-6 -1])
    caxis([-8 -1])
    
    set(gca,'ydir','rev')
    set(gca,'xdir','rev')
    set(gca,'tickdir','out')
    colorbar('location','eastoutside')
    shading interp
    ylabel(Indep_labelstr)
    xlabel('\it k \rm [cpm]')
    
    hold on
    h=plot(diss.K_max(1,:), y);
    set(h(1),'linewidth', 2, 'color','w')
    if any(diss.method(1,:) == 1)
        index = find(diss.method(1,:) == 1);
        h=plot(diss.K_max(1,index), y(index), '*');
        set(h(1),'linewidth', 2, 'color','w')
    end
    hold off
    set(gca,'layer','top')
    title(title_string)
    legend(['\Phi(\itk\rm ) ' shear_list{1}], 'location', 'southeast')
    set(gca,'ylim',y_limits)
    
    if shear_count > 1
        subplot(1,shear_count,2)
        pcolor(diss.K', y, log10(sh2_spec'));grid on
        %set(gca,'XScale','log')
        set(gca,'xlim',[1 150])
        %set(gca,'zlim',[-6 -1])
        caxis([-8 -1])
        
        set(gca,'ydir','rev')
        
        colorbar('location','eastoutside')
        shading interp
        xlabel('\it k \rm [cpm]')
        
        hold on
        h=plot(diss.K_max(2,:), y);
        set(h(1),'linewidth', 2, 'color','w')
        if any(diss.method(2,:) == 1)
            index = find(diss.method(2,:) == 1);
            h=plot(diss.K_max(2,index), y(index), '*');
            set(h(1),'linewidth', 2, 'color','w')
        end
        hold off
        set(gca,'tickdir','out')
        set(gca,'layer','top')
        set(gca,'ylim',y_limits)
        legend(['\Phi(\itk\rm ) ' shear_list{2}], 'location', 'southeast')
        %     title(['\Phi(\itk\rm) ' shear_list{2}])
    end
    
    drawnow
    
end

%-----------------------------------------------------------------
% ----- Figure: frequency spectrogram of shear -------------------
%-----------------------------------------------------------------

if  make_figures && plot_spectrograms

    fig_num = fig_num + 1;
    figure(fig_num); clf
    
    sh1_spec = sh1_spec ./ repmat(diss.speed', size(sh1_spec,1), 1); % scale to preserve variance
    if shear_count>1
        sh2_spec = sh2_spec ./ repmat(diss.speed', size(sh1_spec,1), 1);
    end
    
    x_lim = [1 150]; if tidal_channel, x_lim = [1 250]; end
    
    subplot(1,shear_count,1)
    pcolor(diss.F', y, log10(sh1_spec'));grid on
    %set(gca,'XScale','log')
    set(gca,'xlim',x_lim)
    %set(gca,'zlim',[-6 -1])
    caxis([-8 -1])
    
    set(gca,'ydir','rev')
    set(gca,'xdir','rev')
    set(gca,'tickdir','out')
    colorbar('location','eastoutside')
    shading interp
    ylabel(Indep_labelstr)
    xlabel('\it f \rm [Hz]')
    
    hold on
    h=plot(diss.K_max(1,:)' .* diss.speed, y);
    set(h(1),'linewidth', 2, 'color','w')
    if any(diss.method(1,:) == 1)
        index = find(diss.method(1,:) == 1);
        h=plot(diss.K_max(1,index)' .* diss.speed(index), y(index), '*');
        set(h(1),'linewidth', 2, 'color','w')
    end
    hold off
    set(gca,'layer','top')
    set(gca,'ylim',y_limits)
    title(title_string)
    legend(['\Phi(\itk\rm ) ' shear_list{1}], 'location', 'southeast')
    
    if shear_count >1
        subplot(1,shear_count,2)
        pcolor(diss.F', y, log10(sh2_spec'));grid on
        %set(gca,'XScale','log')
        set(gca,'xlim',x_lim)
        %set(gca,'zlim',[-6 -1])
        caxis([-8 -1])
        
        set(gca,'ydir','rev')
        set(gca,'tickdir','out')
        
        colorbar('location','eastoutside')
        shading interp
        xlabel('\it f \rm [Hz]')
        
        hold on
        h=plot(diss.K_max(2,:)' .* diss.speed, y);
        set(h(1),'linewidth', 2, 'color','w')
        if any(diss.method(2,:) == 1)
            index = find(diss.method(2,:) == 1);
            h=plot(diss.K_max(2,index)' .* diss.speed(index), y(index), '*');
            set(h(1),'linewidth', 2, 'color','w')
        end
        hold off
        set(gca,'layer','top')
        set(gca,'ylim',y_limits)
        %     title(['\Phi(\itk\rm) ' shear_list{2}])
        legend(['\Phi(\itk\rm ) ' shear_list{2}], 'location', 'southeast')
    end
    drawnow
    
end


%-----------------------------------------------------------------
%  ---- Figure: frequency spectrogram of acceleration   ----------
%-----------------------------------------------------------------
if  make_figures && plot_spectrograms

    AA_names = [];
    if ~isempty(AA_names_piezo)
        AA_names = AA_names_piezo;
    elseif ~isempty(AA_names_linear)
        AA_names = AA_names_linear;
    end
    
    fig_num = fig_num + 1;
    figure(fig_num); clf
    
    accel_count = size(diss.AA,2);
    if ~isempty(EMC_Cur), accel_count = accel_count - 1; end % In case we have an EM current meter, accel_count is too big by 1
    c_lim = zeros(  accel_count,2); % used to set color scale
    h     = zeros(1,accel_count);
    
    if accel_count <1, return, end
    if accel_count >0, Ax_spec = squeeze(diss.AA(:,1,1,:));end
    if accel_count >1, Ay_spec = squeeze(diss.AA(:,2,2,:));end
    if accel_count >2, Az_spec = squeeze(diss.AA(:,3,3,:));end
    
    h(1)=subplot(1,accel_count,1);
    pcolor(diss.F', y, log10(Ax_spec'));grid on
    set(gca,'xlim',x_lim)
    c_lim(1,:) = get(gca,'clim');
    
    set(gca,'ydir','rev')
    set(gca,'tickdir','out')
    colorbar('location','eastoutside')
    shading interp
    ylabel(Indep_labelstr)
    xlabel('\it f \rm [Hz]')
    
    set(gca,'layer','top')
    set(gca,'ylim',y_limits)
    legend(['\Phi(\itf\rm ) ' AA_names{1}], 'location', 'southeast')
    title(title_string)
    
    if accel_count > 1
        h(2)=subplot(1,accel_count,2);
        pcolor(diss.F', y, log10(Ay_spec'));grid on
        %set(gca,'XScale','log')
        set(gca,'xlim',x_lim)
        c_lim(2,:) = get(gca,'clim');
        
        %    caxis([0 6])
        
        set(gca,'ydir','rev')
        set(gca,'tickdir','out')
        
        colorbar('location','eastoutside')
        shading interp
        xlabel('\it f \rm [Hz]')
        
        set(gca,'layer','top')
        set(gca,'ylim',y_limits)
        legend(['\Phi(\itf\rm ) ' AA_names{2}], 'location', 'southeast')
        
        %     title(['\Phi(\itf\rm) ' accel_list{2}])
    end
    if accel_count > 2
        h(3)=subplot(1,accel_count,3);
        pcolor(diss.F', y, log10(Az_spec'));grid on
        %set(gca,'XScale','log')
        set(gca,'xlim',x_lim)
        c_lim(3,:) = get(gca,'clim');
        
        %    caxis([0 6])
        
        set(gca,'ydir','rev')
        
        shading interp
        xlabel('\it f \rm [Hz]')
        colorbar('location','eastoutside')
        
        set(gca,'layer','top')
        set(gca,'ylim',y_limits)
        set(gca,'tickdir','out')
        legend(['\Phi(\itf\rm ) ' AA_names{3}], 'location', 'southeast')
        
        %     title(['\Phi(\itf\rm) ' accel_list{3}])
    end
    c_lim = [min(c_lim(:,1)) , max(c_lim(:,2))];
    
    for index = 1:accel_count
        set(h(index),'clim', c_lim)
    end
    
    drawnow
    
end


end
% This is the end of the quick look function
%--------------------------------------------------------------------------

%-----------------------------------------------------------------
% ------ Function to Print figures -------------------------------
%-----------------------------------------------------------------

function printfile(fig_num, name, profile, eps, render)
if nargin < 5
    render = '-painters';
end
if eps
    if ~exist('export_fig','file')
        disp('You are missing the function "export_fig".');
        disp('It can be found in the Matlab file exchange.');
        disp('Please download this function and add to your path');
        disp('to enable exporting figures');
        disp('  ');
        return
    end
    handle = figure(fig_num);
    filename = sprintf('%s_P_%02d_Fig_%02d.pdf',name,profile,fig_num);
    export_fig(handle, filename, render, '-transparent');
    %saveas(fig_num, print_file)
    %fig2pdf(gcf,    print_file, [], '')
end
end


%------------------------------------------------------------------------
% ------  Function Assemble Input Argument ------------------------------
%------------------------------------------------------------------------
% Assemble input arguments into a 1X1 stracture.  Contents of the input
% struture are joined with parameter name / value pairs.  The resulting
% structure can then be saved for future reference.
%-------------------------------------------------------------------------
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