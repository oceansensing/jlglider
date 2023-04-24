%% plot_freq_spec
% Plot frequency spectra from dissipation structure. This function is
% called by ``show_spec.m''
%%
% <latex>\index{Functions!show\_freq\_spec}</latex>
%
%%% Syntax
%   ax = plot_freq_spec( ax, slice, diss, titleString )
%
% * [ax] Axis onto which the spectra should be plotted.  If empty, create a
%       new figure for the plot.
% * [slice] Segment index from the diss structure that should be plotted.
% * [diss] Structure containing dissipation information.  All dissipation
%       information is included but only the segment identified by `slice'
%       is used for the plot.
% * [titleString] String to pre-pended to each figure title. Use it to help
%       identify the current profile.  Set to empty string if not needed.
% * []
% * [ax] Axis that contains the resulting plot.
%
%%% Description
% Plot the frequency spectra of shear and acceleration into the specified 
% axis.
%
% This function is called by ``show_spec''.  The axis input parameter is
% used to allow this function to update the "show_spec" figure.  When 
% called manually, the axis input parameter can be left blank.  This will
% result in a new figure for the plot.

% Version History
%
% 2015-01-30 RGL, Added the ability to handle only a single accelerometer
% 2015-04-21 WID, Updated documentation for publishing.
% 2015-07-29 WID, Modified input arguments to match plot_spec_gui
% 2015-08-27 RGL, Forced removal of zero frequency, set fontsize to 14 for
%                 readability.
% 2015-11-01 WID, Removed hardcoded font size - not cross platform.  Added
%                 note in documentation of plot_spec detailing how users
%                 can change the default font size to suit their specific
%                 computer.
% 2015-11-12 RGL, Removed artificial scaling of piezo acceleromters
% 2015-11-18 RGL, changed ylimits on spectra.
% 2015-11-20 RGL, changed scaling of accelerometers.
% 2015-11-23 RGL, Rename this function to plot_freq_spec to better reflect
%     its purpose.
% 2016-08-30 RGL, Took bold font out of title.
% 2016-08-31 RGL, Enabled detection of profile direction. For 'horizontal',
%     the pressure value is changed to a time (since start of file) value. 
% 2016-12-15 RGL, Changed legend to \nabla, the gradient symbol.
% 2017-04-28 RGL, Check on scalar spectra was inadequate. Need to also
%     check if it is empty. 

function ax = plot_freq_spec( ax, slice, diss, titleString )

titleS = titleString;
    
tiel  = [0   0.5 0.65]';
brown = [0.5 0.5 0]';
fuja  = [0.5 0   0.5]';
lime  = [0.5 0.8 0.5]';
my_colour = [tiel brown fuja lime];

y_lim = [1e-9 1e-1];
if isfield(diss.ql_info, 'op_area')
    if strcmpi(diss.ql_info.op_area,'open_ocean')
        y_lim = [1e-9 1e-1];
    elseif strcmpi(diss.ql_info.op_area,'tidal_ch')
        y_lim = [1e-6 1e1];
    end
end

if isempty(ax),
    f = figure();
    ax = axes('Parent',f);
end

scalar_names = [];
we_have_scalar_spectra = false;

if isfield(diss,'scalar_spectra')
    if ~isempty(diss.scalar_spectra)
        % we have scalar spectra as well
        we_have_scalar_spectra = true; % get the 3D matrix of scalar spectra
        % I assume that the length of the scalar spectra are identical to
        % those of shear, and that the wavenumber and frequencies are the
        % same.
        scalar_names = diss.scalar_vector_list; % The cell array of names for the scalar signals.
    end
end

%%%%
% The main loop

K           = diss.K(:,slice);
F           = diss.F(:,slice);
PD          = diss.P(slice); % pressure in diss structure [dbar].
t           = diss.t(slice); % time since start of file in diss structure [s].
speed       = diss.speed(slice);
K_max       = diss.K_max(:,slice);
T           = diss.T(slice);
profile_dir = diss.vehicle_info.profile_dir;

Nasmyth_spec = diss.Nasmyth_spec(:,:,slice);
e = diss.e(:,slice);

P_string = 'P';
unit_string = 'dBar';
if strcmpi(profile_dir, 'horizontal')
    P_string = 't'; 
    unit_string = 's'; 
    P = t; % use time
else
    P = PD; % use pressure
end

title_string_F = {};
if ~isempty(titleS), title_string_F{end+1} = titleS; end
title_string_F{end+1} = '\rmFrequecy Spectrum';
title_string_F{end+1} = sprintf('Index = %i, %s = %0.1f %s, speed = %0.2f m s^{-1}, T = %0.3f^{\\circ}C', ...
                            slice, P_string, P, unit_string, speed, T);
title_string_F{end+1} = '';
%%%%
%_____________________________
% First the wavenumber Spectra
% Spectra are already in wavenumber form
P_sh1       = diss.sh      (:,1,1,slice);
P_sh1_clean = diss.sh_clean(:,1,1,slice);
e1 = e(1);
phi1 = Nasmyth_spec(:,1);
e1_string = ['\epsilon=' make_scientific(e1,2) 'W kg^{-1}'];
K_max_index_1 = find (K == K_max(1)); % wavenumber limit of spectral integration

e2_string = []; % define but leave empty, if not used
e3_string = [];
e4_string = [];

if size(e,1) > 1
    P_sh2       = diss.sh      (:,2,2,slice);
    P_sh2_clean = diss.sh_clean(:,2,2,slice);
    phi2        = Nasmyth_spec(:,2);
    e2          = e(2);
    e2_string   = ['\epsilon=' make_scientific(e2,2) 'W kg^{-1}'];
    K_max_index_2 = find (K == K_max(2));
end
if size(e,1) > 2
    P_sh3       = diss.sh      (:,3,3,slice);
    P_sh3_clean = diss.sh_clean(:,3,3,slice);
    phi3        = Nasmyth_spec(:,3);
    e3          = e(3);
    e3_string   = ['\epsilon=' make_scientific(e3,2) 'W kg^{-1}'];
    K_max_index_3 = find (K == K_max(3));
end
if size(e,1) > 3
    P_sh4       = diss.sh      (:,4,4,slice);
    P_sh4_clean = diss.sh_clean(:,4,4,slice);
    phi4        = Nasmyth_spec(:,4);
    e4          = e(4);
    e4_string   = ['\epsilon=' make_scientific(e4,2) 'W kg^{-1}'];
    K_max_index_4 = find (K == K_max(4));
end

if we_have_scalar_spectra % extract the scalar spectra for this index
    P_scalar = diss.scalar_spectra.scalar_spec(:,:,slice); 
end

%%%%
% Second the frequency spectra
% get_diss returns wavenumber spectra. Convert to frequency spectra by
% dividing by the mean speed. The spectra are corrected for the wavenumber
% response of the shear probe.
P_sh1       = P_sh1       / speed;
P_sh1_clean = P_sh1_clean / speed;
if size(e,1) > 1
        P_sh2       = P_sh2       / speed;
        P_sh2_clean = P_sh2_clean / speed;
end
if size(e,1) > 2
        P_sh3       = P_sh3       / speed;
        P_sh3_clean = P_sh3_clean / speed;
end
if size(e,1) > 3
        P_sh4       = P_sh4       / speed;
        P_sh4_clean = P_sh4_clean / speed;
end

if isfield(diss, 'piezo') && diss.piezo
% scale down acceleration spectra to bring them on scale
    diss.AA = 1e-6*diss.AA;
end

if size(diss.AA,2)>=1
    P_Ax =diss.AA(:,1,1,slice); % assume not in physical units - uncalibrated.
end
if size(diss.AA,2)>=2
    P_Ay =diss.AA(:,2,2,slice);
end

if size(diss.AA,2)==3
    P_Az = diss.AA(:,3,3,slice);
end

plot_wish_list = {...
    'P_Ax', 'P_Ay', 'P_Az', ...
    'P_sh1',      'P_sh2',      'P_sh3',      'P_sh4'};
legend_wish_list = {...
    'A_x',  'A_y',  'A_z',...
    '\nablau_1',  '\nablau_2',  '\nablau_3',  '\nablau_4'};

plot_index = 0;
plot_list = [];
legend_list = [];
for k = 1:length(plot_wish_list)
    if exist(plot_wish_list{k},'var')
        plot_index = plot_index + 1;
        plot_list{plot_index}   = plot_wish_list{k};
        legend_list{plot_index} = legend_wish_list{k};
    end
end
legend_list{plot_index+1} = 'Nasmyth';

num_vectors = plot_index;
Y_data = zeros(size(F,1),num_vectors);

for k = 1: num_vectors
   eval(['Y_data(:,k) = '  plot_list{k} ';'])
end

F            = F(2:end,:); % Force removal of zero frequency for log-log plots
Y_data       = Y_data(2:end,:);
Nasmyth_spec = Nasmyth_spec(2:end,:);

h  = loglog(ax, F, Y_data, F, Nasmyth_spec/speed, 'k');
grid( ax, 'on' );
hh = legend(ax, legend_list,'location','NorthEastOutside');

for plot_index = 1: size(Y_data,2)
    set(h(plot_index), 'linewidth',1.5)
end

line_count = 0;
if exist('P_Ax','var')
    line_count = line_count + 1;
    set(h(line_count),'linewidth',2.5,'color','c'); % Ax
end
if exist('P_Ay', 'var')
    line_count = line_count + 1;
    set(h(line_count),'linewidth',2.5,'color',0.75*[1 0 1.2]); % Ay
end
shear_1_index = line_count + 1;
if size(diss.AA,2)==3,% we have Az
    shear_1_index = shear_1_index + 1;
    line_count = line_count + 1;
    set(h(line_count),'linewidth',3 ,'color',0.75*[1 1 0])
end
set(h(shear_1_index),    'linewidth',1.5,'color','b'); % sh1
set(h(shear_1_index + 1),'linewidth',1.5,'color','r'); % sh2

%y_limit = [1e-7 1e1]; % limits for y-axis, the way it was

x_limit = [min(F) max(F)];% limits for x-axis
set(ax, 'ylim', y_lim, 'xlim',x_limit)

ylabel(ax, '[Variance Hz^{-1}]')
xlabel(ax, '\it f \rm [Hz]')

%set(gca,'fontsize',14)

title(ax, title_string_F)


end

