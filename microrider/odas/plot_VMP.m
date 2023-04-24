%% plot_VMP
% Simple real-time plots of raw data collected with a VMP.
%%
% <latex>\index{Functions!plot\_VMP}</latex>
%
%%% Syntax
%   plot_VMP( fileName )
%
% * [fileName] name of the data file to display
%
%%% WARNING
% This function is very old and does not make use of modern data files.  It
% is included for users who already use the function and is not recommended
% for new users.
%
%%% Description
% Provides a simple but effective preview of the data from a vertical 
% profile.  Plotted on a time vs count graph, the data forms descending 
% traces where each trace represents a channel.  Visualizing the data in 
% this manner lets one see what is happening to the instrument in either 
% real-time or during a previously recorded acquisition.
%
% When run, the user is asked for the figure duration.  The duration is the 
% length of the vertical axis in seconds.  The function will plot the 
% requested channels over this duration as a single figure.  When the end 
% of the figure is reached, the plot will pause before clearing the graph 
% and continuing from where it left off at the top of the graph.
%
% Data is plotted in units of counts - essentially raw values from the 
% analog to digital converter.  When plotted, the resulting traces tend to 
% overlap and are difficult to see.  To solve this problem one should 
% apply scalar and offset values to position traces to ensure they can be 
% viewed.
%
% The channels to display, along with their respective scalar and offset 
% modifiers, are controlled by variables defined near the top of this 
% function.  The section looks similar to what is shown below:
%
%    % Variables to be plotted, and their channel numbers
%    fast_vars     = {'Ax','Ay','Az','T1_dT1','T2_dT2','Sh1','Sh2','C1_dC1',
%                     'C2_dC2','Fluo','BS'};
%    fast_var_nums = [ 1  2  3  5  7  8  9  12  13  14  15 ];
%    slow_vars     = {'P','P_dP', 'Mx', 'My'};
%    slow_var_nums = [10   11      34    33];
%  
%    set(0,'Defaultaxesfontsize',10,'Defaulttextfontsize',10);
%    Ax_scale        = 1 ;    Ax_offset         =  -30000;         
%    Ay_scale        = 1 ;    Ay_offset         =  -25000;
%                Continued in function....
%
% @image @images/plot_vmp @Example output from $\texttt{plot\_VMP}$.
% @Example output from $\texttt{plot\_VMP}$. A detailed explanation of the
% plot is found within the function description. 
%
% The example plot shows the dark-blue pressure trace on the left side of 
% the figure.  The pre-emphasized pressure trace is blue and is next to the
% pressure trace. At t = 31 s, the pre-emphasized pressure signal separates
% from the pressure signal because the instrument is falling and has a 
% significant rate of change of pressure. This instrument stopped 
% descending at about t = 123 s.  This is evident by the change in the 
% pre-emphasized pressure, a step in the accelerometer signals (the bluish 
% and cyan traces around -30000), and a pulse in the signal from shear 
% probe 1. Shear probe #2 is not installed, nor is thermistor #2.
%
% The main purpose of this function is to give the user a real-time 
% graphical view of the data being collected by a profiler, with no 
% conversion to physical units or other types of processing that might 
% obscure potential problems. The data is shown with all of its warts, such
% as the frequent spikes in the shear signals due to collisions with 
% plankton (t = 76 s).  But the function can also be used to get a quick 
% view of data downloaded from an internally recording instrument. Data 
% courtesy of Manuel Figueroa.
%
% The function shows only the minimum and maximum of consecutive segments
% of the data, in order to plot the data rapidly. The length of the
% segments is based on the pixel resolution of your screen. What is shown
% is an accurate representation of the data. However, the zoom-in function
% will show details that are not real.

% *Version History:*
%
% * 2004-06-06 RGL
% * 2004-09-13 IG  modifications for speed and appearance
% * 2005-01-18 IG? minor modifications to the size & position of the figure
%                  & axes
% * 2005-02-15 IG? no longer asks the user whether or not to continue
%   plotting at the end of each window.  User can interrupt plotting in two
%   ways: Ctrl-C or simply close the window (the routine quits if it can't
%   find the window it created). If no extension is included on the file
%   name, now assume ".p" ending.
% * 2005-02-18 IG? Change routine so that channel names and numbers are
%   provided in the parameters, and locations of each channel within the
%   matrices read in from the file are determined automatically.
% * 2005-02-21 IG? Modified approach to plotting, plotting time for a 100s
%   window cut roughly in half. May be slightly more memory-intensitve,
%   however, and assumes that OpenGL graphics acceleration is available.
% * 2005-07-12 IG? assume 2 uC channels as the "default" instrument setup
% * 2005-07-15 IG? Added clause to deal with variables in the default setup
%   that are not actually present in the data file
% * 2006-05-14 RGL  Added little endian flag to fileopen. Corrected pause at
%   end of file while reading data.
% * 2009-03-06 RGL changed fopen to fopen_odas.
% * 2010-01-15 AWS support odas V6 and up
% * 2011-09-01 AWS added documentation tags for matlab publishing
% * 2012-04-11 WID changed inifile_with_instring calls to setupstr
% * 2012-07-17 WID performance optimizations and improved legend display
% * 2012-09-09 WID updated documentation for matlab publishing.
% * 2013-01-04 WID changed T1_dT1 to C1_dC1 - some instruments might require
%                  T2_dT2
% * 2013-03-01 WID include performance improvements and decimate fix.
% * 2013-06-29 WID real-time instrument fix - indexing byte/word problem
% * 2015-10-31 RGL. Changed documentation.
% * 2016-11-09 RGL. Removed setting of default font size.


function plot_VMP(fileName)

% parsing of command line parameters
if nargin < 1; fileName    = []; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Get Input from the User
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(fileName); fileName = getDataFileName; end;
if (fileName == 'q'); return; end

max_plot_length_in_records = input('How many seconds of data per plot? (default=100s) ');
if isempty(max_plot_length_in_records), max_plot_length_in_records=100; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialization of constants
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variables to be plotted: name, channel_number, scale, offset
% Note: offset values given in units of 1000 counts.
%       fast values:
tmp_f_vars   = {'Ax','Ay','T1_dT1','C1_dC1','Sh1','Sh2'}; 
tmp_f_nums   = [  1    2      5        65      8     9 ];
tmp_f_scale  = [  1    1      1        1       1     1 ];
tmp_f_offset = [-30  -25    -10      -10       0    10 ] * 1000;
%       slow values:
% tmp_s_vars   = { 'P','P_dP', 'PV', 'V_Bat', 'Incl_Y', 'Incl_X'};
% tmp_s_nums   = [ 10    11     12      32        40        41  ];
% tmp_s_scale  = [ 10    10      1       1         1        10  ];
% tmp_s_offset = [-50   -50    -20      20        30        30  ] * 1000;

tmp_s_vars   = { 'P','P_dP', 'JAC_C', 'JAC_T'};
tmp_s_nums   = [ 10    11     49         50];
tmp_s_scale  = [ 10    10      1          1];
tmp_s_offset = [-50   -50     10         20] * 1000;

for x = 1:length(tmp_f_vars),
    fast_vars(x).name   = tmp_f_vars{x};
    fast_vars(x).num    = tmp_f_nums(x);
    fast_vars(x).scale  = tmp_f_scale(x);
    fast_vars(x).offset = tmp_f_offset(x);
end
for x = 1:length(tmp_s_vars),
    slow_vars(x).name   = tmp_s_vars{x};
    slow_vars(x).num    = tmp_s_nums(x);
    slow_vars(x).scale  = tmp_s_scale(x);
    slow_vars(x).offset = tmp_s_offset(x);
end

%set(0,'Defaultaxesfontsize',10,'Defaulttextfontsize',10);

% Factors by which to decimate fast & slow data for plotting (decimation is
% used to speed up plotting of long records)
if max_plot_length_in_records>=500
    decimate_fast = 256;
    decimate_slow = 32;
elseif max_plot_length_in_records>=250
    decimate_fast = 128;
    decimate_slow = 16;
elseif max_plot_length_in_records>=100
	decimate_fast = 64;         
	decimate_slow = 8;
elseif max_plot_length_in_records>=50
    decimate_fast = 8;
    decimate_slow=1;
else
    decimate_fast=1;
    decimate_slow=1;
end
figurePos = [0.005 0.03 0.65 0.89];  % Figure position
figureUnits = 'normalized';
axesPos = [0.07 0.06 0.8 0.85];
xlims = [-60000 40000];

cmap = [0 0 0.5; 0 0 1; 0 1 1; 0 0.5 0; 0 1 0; 0.5 0 0; 1 0 0; 1 0 1; 1 1 0; 0 0 0; 0.5 0.5 0.5];
if length(fast_vars)<=8,
    cmap(9:11,:) = [0 0 0; 0.5 0.5 0.5; 1 1 0];
%           [0 0 0.5; 0 0 1; 0 1 1; 0 0.5 0; 0 1 0; 0.5 0 0; 1 0 0; 1 0 1; 0 0 0; 0.5 0.5 0.5; 1 1 0];
%else
%    cmap = [0 0 0.5; 0 0 1; 0 1 1; 0 0.5 0; 0 1 0; 0.5 0 0; 1 0 0; 1 0 1; 1 1 0; 0 0 0; 0.5 0.5 0.5];
end

% ODAS parameters
size_of_integer =  2;% size of integers from data files in bytes
header_size_i   = 18;% index of header_size in header record
block_size_i    = 19;% index of block_size in header record

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Open the file, get required information
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fileOpen(fileName); % open file to get some info out of its header

% extract some record parameters  
% assume that header has at least 64 2-byte words
HD = fread(fid, 64, 'ushort');	% read header
frewind(fid);			% back to the beginning of the file
header_size_in_bytes = HD(header_size_i); % This is the actual size of the header
header_size = header_size_in_bytes / size_of_integer;
data_record_size_in_bytes = HD(block_size_i)-header_size_in_bytes;
data_record_size = data_record_size_in_bytes / size_of_integer;
header_version = bitshift(HD(11), -8) + bitand(HD(11), 255) /1000;
if header_version >= 6
    first_record_size_in_bytes = header_size_in_bytes + HD(12);
else
    first_record_size_in_bytes = data_record_size_in_bytes + header_size_in_bytes;
end

% ODAS4IR has a bug. The clock frequency is not written into the header
% until the second header.

fseek (fid, first_record_size_in_bytes,'bof');
second_header = fread(fid, 64, 'ushort');
frewind(fid); % back to beginning of file
outbound_clock = second_header(21) + second_header(22)/1000; % Outbound address clock in Hz
record_duration = data_record_size / outbound_clock; % the length of a record in seconds

date_string =[num2str(HD(4)) '-' num2str(HD(5)) '-' num2str(HD(6)) '   ' num2str(HD(7)) ':' ...
        num2str(HD(8)) ':' num2str(HD(9)) '.' num2str(HD(10))];
% Need to get sampling rate

% load channel info by using the address matrix in record number zero of the data file
% compute sizes of setup matrix, total number of channels, etc.
[rows, cols, no_slow_cols, slow_ch, fast_ch] = load_ch_setup(fid);
fast_points_per_record = data_record_size / cols; % number of data points of a fast channel in a single record.
slow_points_per_record = fast_points_per_record / rows; % Number of points of a slow channel in a single record

max_length_of_fast_channels = max_plot_length_in_records * fast_points_per_record;% 
max_length_of_slow_channels = max_plot_length_in_records * slow_points_per_record;% 
for ii=1:length(fast_vars)      % Figure out positions of variables in fast & slow matrices
    fast_vars(ii).ch = find(fast_ch == fast_vars(ii).num);
%    eval([fast_vars{ii} '_ch = find(fast_ch==' num2str(fast_var_nums(ii)) ');']);
end
for ii=1:length(slow_vars)
    slow_vars(ii).ch = find(slow_ch == slow_vars(ii).num);
%    eval([slow_vars{ii} '_ch = find(slow_ch==' num2str(slow_var_nums(ii)) ');']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reading & Plotting
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hF=figure(1);
%set(hF,'units',figureUnits,'position',figurePos,'userdata','VMP_fig','renderer','openGL');
set(hF,'units',figureUnits,'position',figurePos,'userdata','VMP_fig');
clf
h =axes('position',axesPos);

% Setting up the legend.
max_length = 0;                % maximum length of a plot entry name
vars = [fast_vars slow_vars];  % all plot entries as one list

% Find the longest name so we can right justify the index
for var = vars,
    if length(var.name) > max_length, max_length = length(var.name); end
end

% Generate a name and add to the legend.
leg_text = {};
for var = vars,
    name = '';
    for x = 1:max_length-length(var.name), name = [' ' name]; end
    name = [name strrep(var.name, '_', '\_')];
    leg_text{end+1} = sprintf('%3d* %s %+6d', var.scale, name, var.offset);
end

% Next part...

t_f = (0:max_length_of_fast_channels -1)'/ fast_points_per_record;% time vector for plotting of fast channels
t_s = (0:max_length_of_slow_channels -1)'/ slow_points_per_record;% time vector for plotting of slow channels

fast_data = zeros(fast_points_per_record,length(fast_vars))*NaN;   % pre-assign matrices
slow_data = zeros(slow_points_per_record,length(slow_vars))*NaN;

status=fseek(fid,first_record_size_in_bytes,'bof'); % move to the beginning of the first real data record
if status==-1
    close(1)
    error(['Error trying to move to the beginning of the first real data record: ' ferror(fid)]);
end

buffer        = ones(data_record_size+header_size,1)*NaN;

fast_index = 1;
slow_index = 1;
record_counter = 0;

sset = 0;
fset = 0;
lasttime = [0,0,0,0,0,0];    % Empty time
refresh_count = 1;
overflow_count = 0;
draw_refresh_count = 0;

figure(1); cla
ylims = [0 max_plot_length_in_records];
set(h,'clipping', 'off','ylim',ylims,'ydir','rev','xlim',xlims);
set(h,'ColorOrder',cmap);
ylabel('\it t \rm [s]','fontsize',11); 
title([strrep(fileName,'_','\_') ';  ' date_string ' UT'],'fontsize',12,'fontweight','bold'); hold on; grid on
Y_all = zeros(floor(max_length_of_fast_channels/decimate_fast),length(fast_vars))*NaN;  % Initialize plotting matrices
t_all = zeros(floor(max_length_of_fast_channels/decimate_fast),1)*NaN;
Y2_all = zeros(floor(max_length_of_slow_channels/decimate_slow),length(slow_vars))*NaN;
t2_all = zeros(floor(max_length_of_slow_channels/decimate_slow),1)*NaN;

while 1
    % Set attept_reads to 4 records
    attempt_reads = 4*round(4*record_duration);
    
    while attempt_reads,
        % Read a record - header + data
        [buffer, count] = fread(fid, data_record_size+header_size, 'short');
        
        % Read successful, break out of loop.
        if count == data_record_size+header_size, break; end

        % If plotting a previously collected file, we have reached the end
        % of file and should plot and diplay the data.  Place in catch block
        % because variables might not exist yet.
        try
            for ii=1:size(Y_all,2), set(hh1(ii),'xdata',Y_all(:,ii),'ydata',t_all); end
            for ii=1:size(Y2_all,2), set(hh2(ii),'xdata',Y2_all(:,ii),'ydata',t2_all); end
            drawnow
        catch
        end
        
        % must be connected to a real-time instrument, so set the display
        % to refresh on every record.
        refresh_count = 0;
        
        attempt_reads = attempt_reads - 1;
        fseek(fid,-(count*2),'cof');% Move back to where we were in the file.
        pause(1/4);% Pause for 1/4 second.
    end
    
    % If set to 0, we have timed out.  Display an error message.  This is where
    % we exit the main loop.
    if ~attempt_reads
        disp('Pause for header + data exceeded the durtion of 4 records');
        break;
    end
    
    record_counter = record_counter + 1;
    
    [fast_data_all, slow_data_all] = demultiplex(buffer(header_size+1:end), ...
                                                 buffer(1:header_size)); % Put data into matrices
    for ii=1:length(fast_vars)
        junk = fast_data_all(:,fast_vars(ii).ch) * fast_vars(ii).scale + fast_vars(ii).offset;
        if ~isempty(junk)
            fast_data(:,ii) = junk;
        else
            fast_data(:,ii) = ones(size(fast_data_all,1),1)*NaN;
        end
    end
    for ii=1:length(slow_vars)
        junk = slow_data_all(:,slow_vars(ii).ch) * slow_vars(ii).scale + slow_vars(ii).offset;
        if ~isempty(junk)
            junk = junk(:,1); % in case of multiple samples in slow channels
            slow_data(:,ii) = junk;
        else
            slow_data(:,ii) = ones(size(slow_data_all,1),1)*NaN;
        end
    end
    clear fast_data_all slow_data_all
    ii=fast_index:fast_index+fast_points_per_record - 1;
    ii2=slow_index:slow_index+slow_points_per_record - 1;
    fast_index = fast_index + fast_points_per_record;
    slow_index = slow_index + slow_points_per_record;
    Y = decimate_me(fast_data,decimate_fast); % reduce the number of points
    t = decimate_me (t_f(ii),decimate_fast); % also for time vector
    Y2 = decimate_me(slow_data,decimate_slow);
    t2 = decimate_me(t_s(ii2),decimate_slow);
    
% Uncomment the following to disable the decimation of data.  This slows
% things down but provides real results.
%    Y  = fast_data;
%    t  = t_f(ii);
%    Y2 = slow_data;
%    t2 = t_s(ii2);

    fset = fset(end)+1:fset(end) + size(Y,1);
    sset = sset(end)+1:sset(end) + size(Y2,1);
    Y_all(fset,:)=Y;   t_all(fset)=t;     % Add data to the plotting matrices
    Y2_all(sset,:)=Y2; t2_all(sset)=t2;
        
    % On the first record, draw the plots and add the legend.
    if record_counter == 1
        hh1=plot(Y_all,t_all,'linewidth',1);
        hh2=plot(Y2_all,t2_all,'linewidth',1.5);
        l = legend(leg_text,'location','eastoutside');
        set(l, 'FontName', 'Courier');
    end
    
    % To minimize the work, only draw the plot occasionally.  Two methods are
    % used here, the draw_refresh_count and "clock" function.  The "clock"
    % function is slow so draw_refresh_count is used to minimize how many times
    % it is called.
    draw_refresh_count = draw_refresh_count - 1;
    if draw_refresh_count <= 0,
        draw_refresh_count = refresh_count;
        overflow_count = overflow_count + 1;
        tmp_time = clock;
        if etime(tmp_time, lasttime) > 0.5,      % Refresh 2 times a second.
            lasttime = tmp_time;
            % Check the overflow_count - use this to tune the refresh_count variable.
            % The optimal value will differ depending on comptuer speed.
            if overflow_count > 4,
                refresh_count = refresh_count * overflow_count/2;
                disp( [' OPTIMIZE: New refresh_count value: ' num2str(refresh_count)]);
            end
            overflow_count = 0;

            % Set the data to be plotted and then plot the data.
            for ii=1:size(Y_all,2), set(hh1(ii),'xdata',Y_all(:,ii),'ydata',t_all); end
            for ii=1:size(Y2_all,2), set(hh2(ii),'xdata',Y2_all(:,ii),'ydata',t2_all); end
            drawnow
        end
    end
    
    if record_counter == max_plot_length_in_records         % When done a plot, clean up the figure and prepare to plot some more
        % At the end of a screen, plot all remaining data and pause for a
        % second.
        for ii=1:size(Y_all,2), set(hh1(ii),'xdata',Y_all(:,ii),'ydata',t_all); end
        for ii=1:size(Y2_all,2), set(hh2(ii),'xdata',Y2_all(:,ii),'ydata',t2_all); end
        drawnow
        overflow_count = 0;
        pause(1);
        
        figure(1); cla
        Y_all(:,:)=NaN; t_all(:,:)=NaN; Y2_all(:,:)=NaN; t2_all(:,:)=NaN;
        ylims = ylims+max_plot_length_in_records;
        set(h,'clipping', 'off','ylim',ylims,'ydir','rev','xlim',xlims);
    % 	cmap = flipud(ColorCube(length(slow_vars)+length(fast_vars)+2)); cmap=cmap(2:end,:);
        set(h,'ColorOrder',cmap);
        ylabel('\it t \rm [s]','fontsize',11); 
        title([strrep(fileName,'_','\_') ';  ' date_string],'fontsize',12,'fontweight','bold'); hold on; grid on
        t_f = t_f + max_length_of_fast_channels/fast_points_per_record;
        t_s = t_s + max_length_of_slow_channels/slow_points_per_record;
        fast_index = 1;
        slow_index = 1;
        record_counter = 0;
    end    
end

%******************************************************************
function Y = decimate_me(X,R)
%
% function to speed up plotting by showing only max and min
% of successive data ensembles R long
% R is the reduction ratio, a positive integer and
% The number of rows in Y is R times smaller than the rows in X.
% X and Y have the same number of columns 
% 
% RGL Aug. 2004
% WID 2013-03-01  % Now works when R not a power of 2

R = floor(abs(R));% In case R is not a positive integer

% decimate vectors by factor R
len = size(X,1);

% Break matrix into row segments
seg = 1:R*2:len;

% If there is data at the end, use it too.  Now if R > size(X,1) it still works.
if seg(end) ~= len,
    seg(end+1) = len;
end

% Allocate matrix for result
Y = zeros(2*(length(seg)-1),size(X,2));

% Iterate over segments, adding the min and max values to the result.
for i = 1:length(seg)-1,
    Y(2*i-1,:) = min(X(seg(i):seg(i+1),:));
    Y(2*i-0,:) = max(X(seg(i):seg(i+1),:));
end

 
%******************************************************
function [fast_data, slow_data]=demultiplex(data,header)
% [fast_data, slow_data]=demultiplex(data,header)
%
% function to demultiplex data collected with the ODAS system.
% data is the vector of multiplexed data as recorded by ODAS.
% The header is used to figure out the number of slow and fast columns and the number of rows in 
% basic matrix used to multiplex the data.

if (~(nargin ==2 )); error ('input must have 2 arguments'),end;

fast_columns = header(29);% Matrix information is located in the header
slow_columns = header(30);
rows =         header(31);

columns = fast_columns + slow_columns ; % define total number of columns in data
data = reshape(data,columns,[]);
data = data';
fast_data = data(:,slow_columns+1:columns);% extract the fast channels

n_slow = rows*slow_columns; % Number of slow channels
if (slow_columns ~= 0)
    slow_data = data(:,1:slow_columns); % extract slow columns
    slow_data = slow_data';
    slow_data = slow_data(:);
    slow_data = reshape(slow_data,n_slow,[]);
    slow_data = slow_data'; % put slow vectors into columns
else
    slow_data = [];
end

%******************************************************************************
function fileName = getDataFileName
% fileName = getDataFileName returns the file name of the data file.
% Preforms error checking.
%
% Fab, March 1998.
% Modified by RGL to offer the latest '*.p' file for opening.
% 2004-06--6
% 2013-03-01 WID Now uses file_with_ext so case sensitive file systems will
% not cause any problems.

latest_file = get_latest_file;

while 1
   fileName = input(['Enter data file name (default: ' latest_file ' , ''q'' to quit): '], 's');
   if strcmp(fileName, 'q')
      fclose('all');
      break;
   elseif isempty(fileName)
      fileName = latest_file;
   end
   [P,N,E,fileName] = file_with_ext( fileName, {'' '.p' '.P'} );
   if isempty(N),
       fprintf( 1, 'Unable to find file: %s! Try again.\n', fileName );
   else
       break;
   end
end

%******************************************************************************
function fid = fileOpen(fileName)
% fid = fileOpen(fileName) returns the file ID for the file fileName.
% Preforms error checking.
%
% Fab, March 1998.

[fid, error_message] = fopen_odas(fileName, 'r');
if ~isempty(error_message), disp(error_message); end
if fid == -1, error('Error opening file %s !\n', fileName); end


%******************************************************************************
function [nRow, nCol, nSlowCol, slowCh, fastCh] = load_ch_setup(fid)
% [nRow, nCol, nSlowCol, slowCh, fastCh] = load_ch_setup(fid) 
% loads the channel matrix from the file fileName.
% Returns the number of rows in the matrix, the number of columns in the matrix,
% the number of slow columns, a vector containing the slow channels,
% and a vector containing the fast channels.
% 
% Replaces load_ch_setup.m written by L. Zhang.
% Fab, March 1998.
% AWS - 2010-01-14 changes for odas v6 and up

fseek(fid,0,'bof'); % Rewind to the beginning of the file
header = fread(fid, 64, 'ushort');%read the header
header_version = bitshift(header(11), -8) + bitand(header(11), 255) /1000;

if header_version >= 6,
    nFast = header(29);
    nSlowCol = header(30);
    nRow  = header(31);
    nCol = nSlowCol + nFast;
%    matrixSize = nRow * nCol; %use this to check what we get out of the setup file string
    setupfilestr = fread(fid, header(12), '*char*1');
    
    if isempty(setupfilestr)
        error('failed to extract setup file string from first data record');
    end
    setupfilestr = setupfilestr';

    cfg = setupstr(setupfilestr);
    rows = setupstr(cfg, 'matrix', 'row[0-9]+');
    matrix = [];
    for row = rows
      values = textscan(row{1}, '%d16');
      matrix = vertcat(matrix, values{1}');
    end
    
    slowCh = matrix(:, 1:nSlowCol)';
    slowCh = slowCh(:);
    if length(slowCh) ~= nSlowCol*nRow
        error('Error building channel matrix: number of slow channels does not agree.');
    end

    fastCh = matrix(1, nSlowCol+1:nCol)';
    if length(fastCh) ~= nFast
        error('Error building channel matrix: number of fast channels does not aggree.');
    end

else

    nFast = header(29);% extract info
    nSlowCol = header(30);
    nRow = header(31);

    % rebuild the channel matrix
    nCol = nSlowCol + nFast;
    matrixSize = nRow * (nSlowCol+nFast);
    matrix = fread(fid, matrixSize, 'ushort');
    matrix = reshape(matrix, nCol, matrixSize/nCol)';

    % build the slowCh vector
    slowCh = matrix(:, 1:nSlowCol)';
    slowCh = slowCh(:);
    if length(slowCh) ~= nSlowCol*nRow,
        error('Error building channel matrix: number of slow channels does not agree.');
    end

    %build the fastCh vector
    fastCh = matrix(1, nSlowCol+1:nCol)';
    if length(fastCh) ~= nFast,
       error('Error building channel matrix: number of fast channels does not aggree.');
    end
end
