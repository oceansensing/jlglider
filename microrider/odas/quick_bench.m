%% quick_bench
% Quick evaluation of a data file collected while the instrument is on a bench.
%%
% <latex>\index{Functions!quick\_bench}</latex>
%
%%% Syntax
%   quick_bench( 'dataFileName', 'serialNumber', genPDF, genFIG )
%
% * [dataFileName] Name of the file to be processed. If not declared the
%        program will request a file. The default value is the most
%        recently changed file in the current directory.
% * [serialNumber] Serial number of the instrument as a string, or any
%        other useful information. This string is placed into a line of the
%        title of each figure.  If undeclared or left empty, the value of
%        the `SN' parameter within the `instrument_info' section of the
%        configuration file will be used, if available.
% * [genPDF] Logical statement indicating that PDF files of resulting
%        figures should be generated. Requires a ``fig2pdf" program be
%        installed.  Defaults to `true'.
% * [genFIG] Logical statement indicating that FIG files should be saved.
%        Defaults of `false'.
% * []
% * [empty] No return parameters but this function produces two figures.
%
%%% Description
%
% This function generates figures from data collected with a RSI
% instrument, that you wish to test. The instrument should be on a bench, or
% just standing in a laboratory. Dummy probes should be installed when
% collecting data. Data should be collected for a few minutes. This
% function processes the data to produce time-series and spectra of some of
% the channels in the instrument. The resulting graphs allow the user to
% determine if an instrument is working correctly.  
%
% The graphs are primarily used to detect excessive noise within an instrument.
% This helps identify problems such as corroded connections and other faults
% which would otherwise go unseen. For example, you can compare the spectra
% produced by $\texttt{quick\_bench}$ against the noise spectra in the
% calibration report for your instrument.
%
% For real profiles taken in the ocean, use quick_look.m to verify that
% your instrument is working correctly.
%
%%% Examples:
%
%    >> quick_bench( 'data_001.p', '43' )
%
% Plot the data in file $\texttt{data\_001.p}$ collected with the instrument
% that has serial number 43.  The serial number is not required, but
% the string will be added to the title of the figures.
%
% @image @images/quick_bench_1 @The time-series output from the
% $\texttt{quick\_bench}$ function. @The time-series output from the
% $\texttt{quick\_bench}$ function. Quick_bench tries to plot most of the
% variables in your raw data file. These include accelerometers, pressure
% signals, shear probes, thermistors, magnetometers, inclinometers,
% voltage-output oxygen sensors, micro-conductivity sensors, and JAC -T,
% -C, -Turbidity and -Chlorophyll sensors. For some signals, the function
% subtracts the mean and this is indicated in the legend on the right-hand
% side. Only the inclinometer signals are converted into physical units.
% All others remain in raw units of counts. 
%
% @image @images/quick_bench_2 @The spectra output from the
% $\texttt{quick\_bench}$ function.
% @Spectra of some of the signals shown in the time-series figure. The
% instrument should be well cushioned to minimize its vibrations. Even so,
% it is nearly impossible to suppress the output from the extremely
% sensitive accelerometers. AC power line frequency (50/60 Hz)
% contamination is also difficult to suppress and may show up as narrow
% spectral peaks. Dummy probes have been installed in place of the shear
% probes, thermistor and micro-conductivity sensors. Their spectra can be
% directly compared against those in your instrument calibration report to
% check if the noise level is close to that observed at RSI before your
% instrument was shipped.


% Revision History:
% 2007-01-01 (RGL) initial version
% 2007-11-05 (RGL) Added conditional plotting of magnetometer and oxygen sensor
% 2011-04-19 (AWS) support odas v6 and up, added tags for Doxygen, dynamically
%                  build list of vectors and legends.
% 2011-05-04 (RGL) removed most scaling to leave data in raw units.
% 2011-07-26 (RGL) added optional serial number for figure titles
% 2011-09-01 (AWS) added documentation tags for matlab publishing
% 2012-04-11 (WID) changed inifile_with_instring calls to setupstr
% 2012-04-23 (WID) changed plotting functions to provide improved output for
%                  latex input
% 2012-05-02 (RGL) Corrected the counting of windows for the first
%                  figure.
% 2012-11-07 (WID) documentation update
% 2013-01-15 (WID) changed dpres to P_dP for v3.1 files
% 2013-06-10 (WID) updated to use modified read_odas
% 2013-09-10 (WID) modified to make use of texstr
% 2014-03-07 (WID) Fixed bug for when 2 uC probes are present
% 2014-04-03 (WID) Allow working with only 1 shear probe by faking sh2
% 2014-12-15 (RGL) updated to make function test the existance of every
%                  variable to be plotted and to generate a time vector
%                  independent of the existance of any particular channel.
% 2015-01-23 (RGL) Added a figure for JAC sensors, if they exist.
% 2015-02-18 (RGL) Mystery change... variable "n_fft_slow" created.
% 2015-03-02 (WID) Removed file check.  Use check within read_odas.
% 2015-10-30 (RGL) Corrected spectrum of P_dP using fs_slow. Documentation
%                  changes. Added saving the figures as both Matlab figures
%                  and as pdf-files. Set the gca fontsize explicitly to 16.
% 2016-05-24 (WID) Allow slow channel thermistor data to still work.
% 2016-06-07 (WID, RGL) Corrected faulty logic of previous change. Added
%                   time labels to all x-axes.
% 2016-06-09 (WID) Another fix for previous change - fixed the spectra.
% 2016-06-09 (WID) Big changes. Reading most channels now occurs by type.
%                  SN is obtained from the configuration string is
%                  available.  PDF and FIG file generation can be turned
%                  off/on.
% 2016-06-22 (WID) Enable grids.
% 2016-10-14 (RGL) Fixed for Matlab2016a. Removed unncecessary subplot commands.
% 2016-11-21 (RGL) Modified accelerometer section so that this function can
%                  handle the case of both piezo-accelerometers and linear
%                  accelerometers, such as with tidal energy Nemo and for
%                  custom MRs.
% 2017-04-20 (RGL) Corrected ylabel on spectra.
% 2017-06-09 (RGL) Added support for Jac EMC. 
% 2018-02-21 (WL)  Plot JAC-C into two subfigures, JAC_C_V and JAC_C_I
% 2018-02-21 (JMM) Changed sizes of outputted figures
% 2019-05-06 (JMM) Changed figure sizes for EM Current Meter
% 2019-08-21 (JMM) Adjusted all figure sizes so that text size is more
%                  reasonable when exporting to PDF
% 2019-08-27 (JMM) Added plots for rsijac_t, rsijac_c and DO (Rinko)
% 2020-01-30 (JMM) Convert Rinko data to physical units (since the
%                  conversion isn't expected to change because the
%                  calibration coefficients are embedded in the sensor). 
% 2020-05-14 (JMM) Simplistic modification to subtract off mean from Ax
%                   channel for MicroPod systems (because channel is
%                   centered on -6500 instead of 0, as for other instr)
% 2020-06-16 (JMM) Modified ylimits of spectra plot, so that if spectral
%                   values exceed default limits, the axis will be scaled
%                   to fit all data. Added a warning message about axes
%                   limits changing. 
% 2020-06-17 (JMM) Added plots for sbt and sbc channels


function quick_bench(fname, SN, genPDF, genFIG)

fft_length_in_seconds = 2;

% Set the default serial number
if nargin < 2, SN = '___'; end

% Generate PDFs by default.
if nargin < 3, genPDF = true; end

% Do not generate FIGs by default.
if nargin < 3, genFIG = false; end

% Let read_odas sort out the file name.  This is easier and facilitates the
% use of wildcards in the file name.
if nargin < 1
    [variable_list, d] = read_odas(); % convert to a mat-file
else
    [variable_list, d] = read_odas(fname); % convert to a mat-file
end

if isempty(variable_list)
    error(['No data found in data file: ' fname]);
end

% When wildcards are used fname is meaningless.  Must use name returned
% from read_odas.
[filepath, filename, fileext] = fileparts( d.fullPath );


% If the serial number is defined within the instrument_info section, use
% it only if a SN is not explicitly provided.
sn = setupstr(d.cfgobj, 'instrument_info', 'SN');
if ~isempty(sn) && (isempty(SN) || strcmp(SN, '___'))
    SN = sn{1};
end


% __________________________________________________________________
%
% This is where we get some information about the data sampling

if isfield(d,'P_dP') && isfield(d,'P')
    P_hres = deconvolve('P_dP', d.P, d.P_dP, d.fs_slow, d.cfgobj);
end
if isfield(d,'Incl_X'), d.Incl_X = convert_odas(d.Incl_X, 'Incl_X', d.cfgobj);end
if isfield(d,'Incl_Y'), d.Incl_Y = convert_odas(d.Incl_Y, 'Incl_Y', d.cfgobj);end
if isfield(d,'Incl_T'), d.Incl_T = convert_odas(d.Incl_T, 'Incl_T', d.cfgobj);end


%_____________________________________________________________________
% Make figures

n_fft = round(fft_length_in_seconds*d.fs_fast); % length of fft for spectra in points
n_fft_slow = round(n_fft*d.fs_slow/d.fs_fast); % length of fft for spectra of slow channels


figure_file_name = ['QB_' SN '_' filename '_'] ;
figure_num = 0;

date_time_str = datestr(d.filetime, 'yyyy-mm-dd HH:MM UTC');
title_string = texstr( { sprintf('%s; %s', filename, date_time_str)
                         sprintf('SN_%s, Time Series', SN) } );

figure_num = figure_num + 1;
fig_handle = figure(figure_num);
clf(fig_handle)
plt_axes = [];

%-----------------------------------------------------
% -- Determine number of subplots --------------------
%-----------------------------------------------------
windows = 0;
accel = {};
for name = setupstr(d.cfgobj, '', 'type', 'accel')
    accel(end+1) = setupstr(d.cfgobj, name, 'name');
end
if ~isempty(accel), windows = windows + 1; end

piezo = {};
for name = setupstr(d.cfgobj, '', 'type', 'piezo')
    piezo(end+1) = setupstr(d.cfgobj, name, 'name');
end
if ~isempty(piezo), windows = windows + 1; end

shear = {};
for name = setupstr(d.cfgobj, '', 'type', 'shear')
    shear(end+1) = setupstr(d.cfgobj, name, 'name');
end
if ~isempty(shear), windows = windows + 1; end

if isfield(d,'P') || isfield(d,'P_dP') || isfield(d,'P_hres')
    windows = windows + 1;
end
if isfield(d,'Mx') || isfield(d,'My') || isfield(d,'Mz')
    windows = windows + 1;
end
if isfield(d,'O2')
    windows = windows + 1;
end
if isfield(d,'Incl_X') || isfield(d,'Incl_Y') || isfield(d,'Incl_T')
    windows = windows + 1;
end


ucond = {};
for name = setupstr(d.cfgobj, '', 'type', 'ucond')
    parts = strsplit(name{1}, '_d');
    if length(parts) == 2
        ucond(end+1) = setupstr(d.cfgobj, name, 'name');
    end
end
if ~isempty(ucond), windows = windows + 1; end


therm = {};
for name = setupstr(d.cfgobj, '', 'type', 'therm|t_ms')
    parts = strsplit(name{1}, '_d');
    if length(parts) == 2
        therm(end+1) = setupstr(d.cfgobj, name, 'name');
    end
end
if ~isempty(therm), windows = windows + 1; end

%-----------------------------------------------------
% -- Add Figures -------------------------------------
%-----------------------------------------------------
window_index = 0;

% Initialize as empty so plotting routine does not fail.
F                    = [];
F_slow               = [];
F_therm              = [];
F_ucond              = [];

spectra_fast         = [];
spectra_slow         = [];
spectra_therm        = [];
spectra_ucond        = [];

spectra_legend       = {};
spectra_slow_legend  = {};
spectra_therm_legend = {};
spectra_ucond_legend = {};

%
% if ~isempty(accel)
%     window_index = window_index + 1;
%     h(window_index)=subplot(windows,1,window_index);
%     plot_vector = [];
%     legend_string = {};
%     for ax = accel
%         plot_vector = [plot_vector d.(ax{1})];
%         legend_string{end+1} = ax{1};
%         [junk, F] = csd_odas(d.(ax{1}), d.(ax{1}), n_fft, d.fs_fast, [], n_fft/2,'linear');
%         spectra_fast = [spectra_fast junk];
%         spectra_legend{end+1} = ax{1};        
%     end
%     plot(d.t_fast, plot_vector); 
%     legend(legend_string,'location','eastoutside')
%     ylabel ('[counts]')
%     plt_axes = [plt_axes gca];
%     set(gca,'xticklabel',[])
% end
%
% check if accelerometers are fast or slow channels
if ~isempty(accel)
    window_index = window_index + 1;
    h(window_index)=subplot(windows,1,window_index);
    plot_vector = [];
    legend_string = {};
    it_is_fast = [];
    ax = accel;
    
    if length(d.(ax{1})) == length(d.t_fast)
        t = d.t_fast;
        it_is_fast = true;
    elseif length(d.(ax{1})) == length(d.t_slow)
        t = d.t_slow;
        it_is_fast = false;
    else
        t = [];
    end
    
    for ax = accel
        plot_vector = [plot_vector d.(ax{1})];
        legend_string{end+1} = ax{1};
        if ~isempty(it_is_fast) && it_is_fast
            [junk, F] = csd_odas(...
                d.(ax{1}), d.(ax{1}), n_fft, d.fs_fast, [], n_fft/2,'linear');        
            spectra_fast = [spectra_fast junk];
            spectra_legend{end+1} = ax{1};        
        elseif ~isempty(it_is_fast) && ~it_is_fast
            [junk, F_slow] = csd_odas(...
                d.(ax{1}), d.(ax{1}), n_fft_slow, d.fs_slow, [], n_fft_slow/2,'linear');
             spectra_slow = [spectra_slow junk]; 
             spectra_slow_legend{end+1} = ax{1};        
        end
    end
    plot(t, plot_vector); 
    legend(legend_string,'location','eastoutside')
    ylabel ('[counts]')
    plt_axes = [plt_axes gca];
    set(gca,'xticklabel',[])
end
%
% check if piezo-accelerometers are fast or slow channels
if ~isempty(piezo)
    window_index = window_index + 1;
    h(window_index)=subplot(windows,1,window_index);
    plot_vector = [];
    legend_string = {};
    it_is_fast = [];
    ax = piezo;
    
    if length(d.(ax{1})) == length(d.t_fast)
        t = d.t_fast;
        it_is_fast = true;
    elseif length(d.(ax{1})) == length(d.t_slow)
        t = d.t_slow;
        it_is_fast = false;
    else
        t = [];
    end
        
    for ax = piezo
        mean_A = mean(d.(ax{1}));
        if abs(mean_A)<500 % channel centered on zero
            plot_vector = [plot_vector d.(ax{1})];
            legend_string{end+1} = ax{1};
        else % channel not centered on zero (i.e. seaglider)
            plot_vector = [plot_vector d.(ax{1})-round(mean_A)];
            if mean_A>0
                legend_string{end+1} = [ax{1},' - ',num2str(round(mean_A))];
            else 
                legend_string{end+1} = [ax{1},' + ',num2str(round(-mean_A))];
            end
        end
        
        if ~isempty(it_is_fast) && it_is_fast
            [junk, F] = csd_odas(...
                d.(ax{1}), d.(ax{1}), n_fft, d.fs_fast, [], n_fft/2,'linear');        
            spectra_fast = [spectra_fast junk];
            spectra_legend{end+1} = ax{1};        
        elseif ~isempty(it_is_fast) && ~it_is_fast
            [junk, F_slow] = csd_odas(...
                d.(ax{1}), d.(ax{1}), n_fft_slow, d.fs_slow, [], n_fft_slow/2,'linear');
             spectra_slow = [spectra_slow junk]; 
             spectra_slow_legend{end+1} = ax{1};        
        end
    end
    plot(t, plot_vector); 
    legend(legend_string,'location','eastoutside')
    ylabel ('[counts]')
    plt_axes = [plt_axes gca];
    set(gca,'xticklabel',[])
end
%
%____________
if ~isempty(ucond)
    skip_size = length(d.t_fast) / length(d.(ucond{1}));
    
    if mod(skip_size,1) ~= 0
        warning(['Skipping micro conductivity channels. ' ...
                 'Invalid positionning of channels in address matrix.\n' ...
                 'uConductivity samples are not evenly spaced.']);
    else
        ucondfast = d.fs_fast / skip_size;
        ucondfft  = n_fft / skip_size;
                
        window_index = window_index + 1;
        h(window_index)=subplot(windows,1,window_index);
        plot_vector = [];
        legend_string = {};
        for uc = ucond
            mean_C = mean(d.(uc{1}));
            plot_vector = [plot_vector d.(uc{1}) - mean_C];
            legend_string{end+1} = texstr(sprintf('%s%+.0f', uc{1}, -mean_C));
            [junk, F_ucond]  = csd_odas(d.(uc{1}), d.(uc{1}), ucondfft, ucondfast, [], ucondfft/2,'linear');
            spectra_ucond = [spectra_ucond junk];
            name = strsplit(uc{1}, '_d');
            spectra_ucond_legend{end+1} = name{1};
        end
        
        t_fast_ucond = d.t_fast(1:skip_size:length(d.t_fast));
        plot(t_fast_ucond, plot_vector);
        legend(legend_string,'location','eastoutside')
        ylabel ('[counts]')
        plt_axes = [plt_axes gca];
        set(gca,'xticklabel',[])
    end
end


%
if isfield(d,'Mx') || isfield(d,'My') || isfield(d,'Mz')
    window_index = window_index + 1;
    h(window_index)=subplot(windows,1,window_index);
    plot_vector = [];
    legend_string = {};
    if isfield(d,'Mx')
        plot_vector = [plot_vector d.Mx];
        legend_string{end+1} = 'M_x';
    end
    if isfield(d,'My')
        plot_vector = [plot_vector d.My];
        legend_string{end+1} = 'M_y';
    end
    if isfield(d,'Mz')
        plot_vector = [plot_vector d.Mz];
        legend_string{end+1} = 'M_z';
    end
    plot(d.t_slow, plot_vector); 
    legend(legend_string, 'location','eastoutside')
    ylabel ('[counts]')
    plt_axes = [plt_axes gca];
    set(gca,'xticklabel',[])
end


%
if isfield(d,'O2')
    window_index = window_index + 1;
    h(window_index)=subplot(windows,1,window_index);
    plot(d.t_slow, d.O2); 
    legend('O_2 [Hz]', 'location','eastoutside')
    plt_axes = [plt_axes gca];
    set(gca,'xticklabel',[])
end


%
if isfield(d,'Incl_X') || isfield(d,'Incl_Y') || isfield(d,'Incl_T')
    window_index = window_index + 1;
    h(window_index)=subplot(windows,1,window_index);
    plot_vector = [];
    legend_string = {};
    if isfield(d,'Incl_X')
        plot_vector = [plot_vector d.Incl_X];
        legend_string{end+1} = 'Incl\_X';
    end
    if isfield(d,'Incl_Y')
        plot_vector = [plot_vector d.Incl_Y];
        legend_string{end+1} = 'Incl\_Y';
    end
    if isfield(d,'Incl_T')
        plot_vector = [plot_vector d.Incl_T];
        legend_string{end+1} = 'Incl\_T';
    end
    plot(d.t_slow, plot_vector); 
    legend(legend_string, 'location','eastoutside')
    ylabel ('[degrees]')
    plt_axes = [plt_axes gca];
    set(gca,'xticklabel',[])
end


%
if ~isempty(therm)
    % Sometimes the thermister channel is not fast.  Trim the time vector
    % accordingly.
    skip_size = length(d.t_fast) / length(d.(therm{1}));
    
    if mod(skip_size,1) ~= 0
        warning(['Skipping thermisters. ' ...
                 'Invalid positionning of channels in address matrix.\n' ...
                 'Thermister samples are not evenly spaced.']);
    else
        thermfast = d.fs_fast / skip_size;
        thermfft  = n_fft / skip_size;
                
        window_index = window_index + 1;
        h(window_index)=subplot(windows,1,window_index);
        plot_vector = [];
        legend_string = {};
        for t = therm
            mean_T = mean(d.(t{1}));
            plot_vector = [plot_vector d.(t{1}) - mean_T];
            legend_string{end+1} = texstr(sprintf('%s%+.0f', t{1}, -mean_T));
            [junk, F_therm]  = csd_odas(d.(t{1}), d.(t{1}), thermfft, thermfast, [], thermfft/2,'linear');
            spectra_therm = [spectra_therm junk];
            name = strsplit(t{1}, '_d');
            spectra_therm_legend{end+1} = name{1};
        end
        
        t_fast_therm = d.t_fast(1:skip_size:length(d.t_fast));
        plot(t_fast_therm, plot_vector);
        legend(legend_string,'location','eastoutside')
        ylabel ('[counts]')
        plt_axes = [plt_axes gca];
        set(gca,'xticklabel',[])
    end
end


%
if ~isempty(shear)
    window_index = window_index + 1;
    h(window_index)=subplot(windows,1,window_index);
    plot_vector = [];
    legend_string = {};
    for sh = shear
        plot_vector = [plot_vector d.(sh{1})];
        legend_string{end+1} = sh{1};
        [junk, F] = csd_odas(d.(sh{1}), d.(sh{1}), n_fft, d.fs_fast, [], n_fft/2,'linear');
        spectra_fast = [spectra_fast junk];
        spectra_legend{end+1} = sh{1};        
    end
    plot(d.t_fast, plot_vector); 
    legend(legend_string,'location','eastoutside')
    ylabel ('[counts]')
    plt_axes = [plt_axes gca];
    set(gca,'xticklabel',[])
end


%
if isfield(d,'P') || isfield(d,'P_dP')
    window_index = window_index + 1;
    h(window_index)=subplot(windows,1,window_index);
    plot_vector = [];
    legend_string = {};
    if isfield(d,'P')
        plot_vector = d.P;
        legend_string{end+1} = 'P';
    end
    if isfield(d,'P_dP')
        plot_vector = [plot_vector d.P_dP];
        legend_string{end+1} = 'P\_dP';
        [junk, F_slow]  = ...
            csd_odas(d.P_dP, d.P_dP, n_fft_slow, d.fs_slow, [], n_fft_slow/2,'linear');
        spectra_slow = [spectra_slow junk];
        spectra_slow_legend{end+1} = 'P\_dP';
    end

    plot(d.t_slow, plot_vector); 
    legend(legend_string, 'location','eastoutside')
    ylabel ('[counts]')
    plt_axes = [plt_axes gca];
    set(gca,'xticklabel',[])
end

title(h(1),title_string)

xlabel(h(end),'t [s]')
set(gca,'xticklabelmode','auto')
 
position = zeros(window_index,4);
for k = 1:window_index
    position(k,:) = get(h(k), 'position');
end

min_width = min(position(:,3));
position(:,3) = 0.95*min_width;% 100% causes problems on right edge

for k=1:window_index
    set(h(k), 'position',position(k,:))
end

% Turn on grids for each plot.
for i = 1:length(plt_axes)
    grid(plt_axes(i), 'on');
end

linkaxes(plt_axes,'x')

% min_extents = inf;
% for i = 1:length(plt_axes)
%     lim = get(plt_axes(i), 'xlim');
%     if lim(2) < min_extents, min_extents = lim(2); end
% end
% for i = 1:length(plt_axes)
%     set(plt_axes(i), 'xlim', [0 min_extents]);
%     set(plt_axes(i), 'ActivePositionProperty', 'outerposition');
% end

fig_name = [figure_file_name '_Fig_' num2str(figure_num)];
if genPDF && exist('fig2pdf', 'file')
    fig2pdf(fig_handle, fig_name, [10,6.25]);
end
if genFIG, saveas(fig_handle, fig_name, 'fig'); end

for i = 1:length(plt_axes)
    set(plt_axes(i), 'xticklabelmode', 'auto');
end

%
% Plot if JAC C, T Chlorophyll or Turbidity sensors exist
windows = 0;
if isfield(d,'Turbidity')
    windows = windows + 1;
end
if isfield(d,'Chlorophyll')
    windows = windows + 1;
end
if isfield(d,'JAC_T')
    windows = windows + 1;
end
if isfield(d,'JAC_C')
    windows = windows + 2;
end
window_index = 0;
if windows > 0
    figure_num = figure_num + 1;
    fig_handle = figure(figure_num);
    clf(fig_handle)
    plt_axes = [];
    
    if isfield(d,'Turbidity')
        window_index = window_index + 1;
        h(window_index)=subplot(windows,1,window_index);
        plot(d.t_fast, d.Turbidity); 
        legend('Turbidity', 'location','eastoutside')
        plt_axes = [plt_axes gca];
        set(gca,'xticklabel',[])
        ylabel('[counts]')
    end
    if isfield(d,'Chlorophyll')
        window_index = window_index + 1;
        h(window_index)=subplot(windows,1,window_index);
        plot(d.t_fast, d.Chlorophyll); 
        legend('Chlorophyll', 'location','eastoutside')
        plt_axes = [plt_axes gca];
        set(gca,'xticklabel',[])
        ylabel('[counts]')
    end
    if isfield(d,'JAC_T')
        window_index = window_index + 1;
        h(window_index)=subplot(windows,1,window_index);
        plot(d.t_slow, d.JAC_T); 
        legend('JAC\_T', 'location','eastoutside')
        plt_axes = [plt_axes gca];
        set(gca,'xticklabel',[])
        ylabel('[counts]')
    end
%     if isfield(d,'JAC_C')
%         window_index = window_index + 1;
%         h(window_index)=subplot(windows,1,window_index);
%         plot(d.t_slow, d.JAC_C); 
%         legend('JAC\_C', 'location','eastoutside')
%         plt_axes = [plt_axes gca];
%         set(gca,'xticklabel',[])
%         ylabel('[counts]')
%     end
    
    
     if isfield(d,'JAC_C')
            window_index = window_index + 1;
            h(window_index)=subplot(windows,1,window_index);
            plot(d.t_slow, bitshift(d.JAC_C, -16)); 

            legend('JAC\_C\_I', 'location','eastoutside')
            plt_axes = [plt_axes gca];
            set(gca,'xticklabel',[])
            ylabel('[counts]')
     end
     if isfield(d,'JAC_C')
            window_index = window_index + 1;
            h(window_index)=subplot(windows,1,window_index);
            plot(d.t_slow, rem(d.JAC_C,2^16)); 
            legend('JAC\_C\_V', 'location','eastoutside')
            plt_axes = [plt_axes gca];
            set(gca,'xticklabel',[])
            ylabel('[counts]')
     end
     
     
    subplot(windows,1,1)
    title(title_string)

    subplot(windows,1,windows)
    xlabel('t [s]')
    set(gca,'xticklabelmode','auto')

    position = zeros(window_index,4);
    for k = 1:window_index
        position(k,:) = get(h(k), 'position');
    end

    min_width = min(position(:,3));
    position(:,3) = 0.95*min_width;% 100% causes problems on right edge

    for k=1:window_index
        set(h(k), 'position',position(k,:))
    end
    
    % Turn on grids for each plot.
    for i = 1:length(plt_axes)
        grid(plt_axes(i), 'on');
    end


    fig_name = [figure_file_name '_Fig_' num2str(figure_num)];
    if genPDF && exist('fig2pdf', 'file')
        fig2pdf(fig_handle, fig_name, [10,6.25]);
    end
    if genFIG, saveas(fig_handle, fig_name, 'fig'); end
    
    for i = 1:length(plt_axes)
        set(plt_axes(i), 'xticklabelmode', 'auto');
    end
end

%----------------------------------------------------------------
% Plot if JAC C, T sampled using rsijac board
windows = 0;

if isfield(d,'RSIJAC_T')
    windows = windows + 1;
end
if isfield(d,'RSIJAC_C')
    windows = windows + 1;
end
window_index = 0;
if windows > 0
    figure_num = figure_num + 1;
    fig_handle = figure(figure_num);
    clf(fig_handle)
    plt_axes = [];
    
    if isfield(d,'RSIJAC_T')
        window_index = window_index + 1;
        h(window_index)=subplot(windows,1,window_index);
        plot(d.t_slow, d.RSIJAC_T);
        legend('JAC\_T', 'location','eastoutside')
        plt_axes = [plt_axes gca];
        set(gca,'xticklabel',[])
        ylabel('[counts]')
    end

    
    
     if isfield(d,'RSIJAC_C')
            window_index = window_index + 1;
            h(window_index)=subplot(windows,1,window_index);
            plot(d.t_slow, d.RSIJAC_C);
            legend('JAC\_C', 'location','eastoutside')
            plt_axes = [plt_axes gca];
            set(gca,'xticklabel',[])
            ylabel('[counts]')
     end
     
    subplot(windows,1,1)
    title(title_string)

    subplot(windows,1,windows)
    xlabel('t [s]')
    set(gca,'xticklabelmode','auto')

    position = zeros(window_index,4);
    for k = 1:window_index
        position(k,:) = get(h(k), 'position');
    end

    min_width = min(position(:,3));
    position(:,3) = 0.95*min_width;% 100% causes problems on right edge

    for k=1:window_index
        set(h(k), 'position',position(k,:))
    end
    
    % Turn on grids for each plot.
    for i = 1:length(plt_axes)
        grid(plt_axes(i), 'on');
    end


    fig_name = [figure_file_name '_Fig_' num2str(figure_num)];
    if genPDF && exist('fig2pdf', 'file')
        fig2pdf(fig_handle, fig_name, [10,6.25]);
    end
    if genFIG, saveas(fig_handle, fig_name, 'fig'); end
    
    for i = 1:length(plt_axes)
        set(plt_axes(i), 'xticklabelmode', 'auto');
    end
end

%----------------------------------------------------------------
% SBT / SBC fields
windows = 0;

if isfield(d,'SBT1') % Question: Is this always the channel name?
    windows = windows + 1;
end
if isfield(d,'SBC1')
    windows = windows + 1;
end
window_index = 0;
if windows > 0
    figure_num = figure_num + 1;
    fig_handle = figure(figure_num);
    clf(fig_handle)
    plt_axes = [];
    
    if isfield(d,'SBT1')
        window_index = window_index + 1;
        h(window_index)=subplot(windows,1,window_index);
        
        SBT1_physical = convert_odas(d.SBT1, 'SBT1', d.cfgobj);
        plot(d.t_slow, SBT1_physical);
        legend('SBT1', 'location','eastoutside')
        plt_axes = [plt_axes gca];
        set(gca,'xticklabel',[])
        ylabel('[\circ C]')
    end

    
    
     if isfield(d,'SBC1')
            window_index = window_index + 1;
            h(window_index)=subplot(windows,1,window_index);
            plot(d.t_slow, d.SBC1);
            legend('SBC1', 'location','eastoutside')
            plt_axes = [plt_axes gca];
            set(gca,'xticklabel',[])
            ylabel('[counts]')
     end
     
    subplot(windows,1,1)
    title(title_string)

    subplot(windows,1,windows)
    xlabel('t [s]')
    set(gca,'xticklabelmode','auto')

    position = zeros(window_index,4);
    for k = 1:window_index
        position(k,:) = get(h(k), 'position');
    end

    min_width = min(position(:,3));
    position(:,3) = 0.95*min_width;% 100% causes problems on right edge

    for k=1:window_index
        set(h(k), 'position',position(k,:))
    end
    
    % Turn on grids for each plot.
    for i = 1:length(plt_axes)
        grid(plt_axes(i), 'on');
    end


    fig_name = [figure_file_name '_Fig_' num2str(figure_num)];
    if genPDF && exist('fig2pdf', 'file')
        fig2pdf(fig_handle, fig_name, [10,6.25]);
    end
    if genFIG, saveas(fig_handle, fig_name, 'fig'); end
    
    for i = 1:length(plt_axes)
        set(plt_axes(i), 'xticklabelmode', 'auto');
    end
end

%---
% Handle a Rinko-FT 

if isfield(d,'DO')  
    % Convert to physical units
     [d.DO  units_DO] = convert_odas(d.DO, 'DO', d.cfgobj);
     [d.DO_T units_T] = convert_odas(d.DO_T, 'DO_T', d.cfgobj);

    % Make plot
    figure_num = figure_num + 1;
    fig_handle = figure(figure_num);
    clf(fig_handle)
    
    h(1) = subplot(2,1,1);
    plot(d.t_slow, d.DO);
    grid('on');
    legend('DO', 'location','eastoutside');
    title(title_string)
%     ylabel('[counts]')
    ylabel(units_DO)
    xlabel('\it t \rm [s]')
    set(gca,'xticklabelmode','auto')
    
    h(2) = subplot(2,1,2);
    plot(d.t_slow, d.DO_T);
    grid('on');
    legend('DO\_T', 'location','eastoutside');
%     ylabel('[counts]')
    ylabel(units_T)
    xlabel('\it t \rm [s]')
    set(gca,'xticklabelmode','auto')
    
    position = zeros(2,4);
    for k = 1:2
        position(k,:) = get(h(k), 'position');
    end
    
    min_width = min(position(:,3));
    position(:,3) = 0.95*min_width;% 100% causes problems on right edge
    
    for k=1:2
        set(h(k), 'position',position(k,:))
    end
    
    fig_name = [figure_file_name '_Fig_' num2str(figure_num)];
    if genPDF && exist('fig2pdf', 'file')
        fig2pdf(fig_handle, fig_name, [10,6.25]);
    end
    if genFIG, saveas(fig_handle, fig_name, 'fig'); end
end
%----------------------------------------------------------------
%spectra

title_string = texstr( { sprintf('%s; %s', filename, date_time_str)
                         sprintf('SN_%s, Spectra', SN) } );

figure_num = figure_num + 1;
fig_handle = figure(figure_num);
clf(fig_handle)

loglog(F, spectra_fast, ...
       F_therm, spectra_therm, ...
       F_ucond, spectra_ucond, ...
       F_slow, spectra_slow);
   
grid('on');

legend([spectra_legend spectra_therm_legend spectra_ucond_legend spectra_slow_legend], ...
       'location','eastoutside');

title(title_string)

xlim_lower = 0.9/fft_length_in_seconds;
xlim_upper = 1.1*d.fs_fast/2;
set(gca,'xlim',[xlim_lower xlim_upper])

ylim_lower = 1e-4;
ylim_upper = 1e2;
if max(max(spectra_fast))<ylim_upper % if spectra values are less than 100
    set(gca, 'ylim', [ylim_lower ylim_upper])
else % if spectral values are above 100, use default Matlab limit for ymax
    ylim_Matlab = get(gca,'ylim');
    set(gca, 'ylim', [ylim_lower ylim_Matlab(2)])
    if ylim_Matlab(2)/ylim_upper>10
        warning('ASTP spectra: Ylimits have been increased by more than a factor of 10 from default range. Spectral magnitudes may be higher than desired.')
    end
end
ylabel('[counts^2 Hz^{-1}]')
xlabel('\it f \rm [Hz]')

fig_name = [figure_file_name '_Fig_' num2str(figure_num)];
if genPDF && exist('fig2pdf', 'file')
    fig2pdf(fig_handle, fig_name, [10,6.25]);
end
if genFIG, saveas(fig_handle, fig_name, 'fig'); end

%----------------------------------------
% Handle a JAC Electro-magnetic current meter
% First the current to the sensor

if isfield(d,'EMC_Cur')  || isfield(d,'EM_Cur') % are measuring the current to the EMC transducer

    figure_num = figure_num + 1;
    fig_handle = figure(figure_num);
    clf(fig_handle)
    title_string = texstr( { sprintf('%s; %s', filename, date_time_str)
                         sprintf('SN_%s', SN) } );
    
    if isfield(d,'EM_Cur'), d.EMC_Cur = d.EM_Cur; end
%-----
    h(1) = subplot(3,1,1);
    
    plot(d.t_fast, d.EMC_Cur);
    
    grid('on');
    legend('EMC\_Cur', 'location','eastoutside');
    title(title_string)
    ylabel('[counts]')
    xlabel('\it t \rm [s]')
    set(gca,'xticklabelmode','auto')
    
%-----
    h(2) = subplot(3,1,2); % show only the first second of data
    
    plot(d.t_fast, d.EMC_Cur);
    
    grid('on');
    legend('EMC\_Cur', 'location','eastoutside');
    set(gca, 'xlim', [0 1])
    ylabel('[counts]')
    xlabel('\it t \rm [s]')
    set(gca,'xticklabelmode','auto')
    
%-----
    h(3) = subplot(3,1,3);
    
    fft_length_in_seconds = 4;
    n_fft = round(fft_length_in_seconds*d.fs_fast); % length of fft for spectra in points
    n_fft_slow = round(n_fft*d.fs_slow/d.fs_fast); % length of fft for spectra of slow channels

    d.EMC_Cur = d.EMC_Cur - mean(d.EMC_Cur); % subtract mean for better looking spectrum
    [P_EMC_Cur, F] = csd_odas(d.EMC_Cur, [], n_fft, d.fs_fast, [], n_fft/2,'linear');
    
    loglog(F, P_EMC_Cur)
    set(gca,'ylim', [1e-5 1e10])
    set(gca,'xlim', [0.9*F(2) 1.1*F(end)])
    
    grid('on');
    legend('EMC\_Cur', 'location','eastoutside');
    ylabel('[counts^2 Hz^{-1}]')
    xlabel('\it f \rm [Hz]')
    set(gca,'xticklabelmode','auto')
    
    
%     position = zeros(3,4);
%     for k = 1:3
%         position(k,:) = get(h(k), 'position');
%     end
%     
%     min_width = min(position(:,3));
%     position(:,3) = 0.95*min_width;% 100% causes problems on right edge
%     
%     for k=1:3
%         set(h(k), 'position',position(k,:))
%     end
    
    fig_name = [figure_file_name '_Fig_' num2str(figure_num)];
    if genPDF && exist('fig2pdf', 'file')
        fig2pdf(fig_handle, fig_name, [10,6.25]);
    end
    if genFIG, saveas(fig_handle, fig_name, 'fig'); end
end

% Handle a JAC Electro-magnetic current meter
% Second the signal from the sensor itself


if isfield(d,'U_EM')  % are measuring the current to the EMC transducer
    
    figure_num = figure_num + 1;
    fig_handle = figure(figure_num);
    clf(fig_handle)
    
    h(1) = subplot(2,1,1);
    
    plot(d.t_slow, d.U_EM);
    
    grid('on');
    legend('U\_EM', 'location','eastoutside');
    title(title_string)
    ylabel('[counts]')
    xlabel('\it t \rm [s]')
    set(gca,'xticklabelmode','auto')
    
    h(2) = subplot(2,1,2);
    
    [P_U_EM, F] = csd_odas(...
        d.U_EM, [], n_fft_slow, d.fs_slow, [], n_fft_slow/2,'linear');
    
    loglog(F, P_U_EM)
    
    grid('on');
    legend('U\_EM', 'location','eastoutside');
    ylabel('[counts^2 Hz^{-1}]')
    xlabel('\it f \rm [Hz]')
    set(gca,'xticklabelmode','auto')
    
    
    position = zeros(2,4);
    for k = 1:2
        position(k,:) = get(h(k), 'position');
    end
    
    min_width = min(position(:,3));
    position(:,3) = 0.95*min_width;% 100% causes problems on right edge
    
    for k=1:2
        set(h(k), 'position',position(k,:))
    end
    
    fig_name = [figure_file_name '_Fig_' num2str(figure_num)];
    if genPDF && exist('fig2pdf', 'file')
        fig2pdf(fig_handle, fig_name, [10,6.25]);
    end
    if genFIG, saveas(fig_handle, fig_name, 'fig'); end
end


