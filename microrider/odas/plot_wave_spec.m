%% plot_wave_spec
% Plot wave spectra from dissipation structure
%%
% <latex>\index{Functions!show\_wave\_spec}</latex>
%
%%% Syntax
%   ax = plot_wave_spec( ax, slice, diss, titleString )
%
% * [ax] Axis onto which the spectra should be plotted.  If empty, create a
%       new figure for the plot.
% * [slice] Segment index from the diss structure that should be plotted.
% * [diss] Structure containing dissipation information.  All dissipation
%       information is included but only the segment identified by 'slice'
%       is used for the plot.
% * [titleString] String to pre-pended to each figure title. Use it to help
%       identify the current profile.  Set to empty string if not needed.
% * []
% * [ax] Axis that contains the resulting plot.
%
%%% Description
% Plot the wave spectra of shear and acceleration into the specified axis.
%
% This function is called by "show_spec".  The axis input parameter is
% used to allow this function to update the "show_spec" figure.  When 
% called manually, the axis input parameter can be left blank.  This will
% result in a new figure for the plot.

% Version History
%
% 2015-01-31 RGL, corrected bug with plotting of scalart spectra
% 2015-04-21 WID, Updated documentation for publishing.
% 2015-07-29 WID, Modified input arguments to match plot_spec_gui
% 2015-08-27 RGL, Forced removal of zero wavenumber, set fontsize to 14 for
%   readability. 
% 2015-11-01 WID, Removed hardcoded font size - not cross platform.  Added
%                 note in documentation of plot_spec detailing how users
%                 can change the default font size to suit their specific
%                 computer.
% 2015-11-19 RGL, Added MAD to title.
% 2015-11-23 RGL, Rename this function to plot_wave_spec to better reflect
%     its purpose.
% 2016-08-30 RGL, Replace MAD with FM (figure of merit) where 
%     FM = MAD*sqrt(DOF), which should be independent of the number of
%     fft-segments used to make a spectral estimate. Also, changed title
%     string for wavenumber spectrum.
% 2016-12-15 RGL, Changed legend to \nabla, the gradient symbol.
% 2017-04-28 RGL, Check on scalar spectra was inadequate. Need to also
%     check if it is empty. 

function  ax = plot_wave_spec( ax, slice, diss, titleString )

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
    elseif  strcmpi(diss.ql_info.op_area,'tidal_ch')
        y_lim = [1e-6 1e1];
    end
end

if isempty(ax),
    f = figure();
    ax = axes('Parent',f);
end

PD = diss.P; % pressure in diss structure
time = 0;

scalar_names = [];
we_have_scalar_spectra = false;

if isfield(diss,'scalar_spectra') % we have scalar spectra as well
    if ~isempty(diss.scalar_spectra)
        we_have_scalar_spectra = true; % get the 3D matrix of scalar spectra
        % I assume that the length of the scalar spectra are identical to
        % those of shear, and that the wavenumber and frequencies are the
        % same.
        scalar_names = diss.scalar_vector_list; % The cell array of names for the scalar signals.
    end
end

    
%%
% The main loop

e        = diss.e        (:,slice);
K        = diss.K        (:,slice);
F        = diss.F        (:,slice);
method   = diss.method   (:,slice);
mad      = diss.mad      (:,slice);
FM       = diss.FM       (:,slice);
dof_e    = diss.dof_e    (:,slice);
dof_spec = diss.dof_spec;
speed    = diss.speed    (slice);
P        = diss.P        (slice);
K_max    = diss.K_max    (:,slice);
T        = diss.T        (slice);
Nasmyth_spec = diss.Nasmyth_spec(:,:,slice);

P_string = 'P';
unit_string = 'dBar';
if time, P_string = 't'; unit_string = 's'; end

title_string_K = {};
if ~isempty(titleS), title_string_K{end+1} = titleS; end
title_string_K{end+1} = ['\rm Wavenumber Spectrum'];
title_string_K{end+1} = sprintf('Index = %i, %s = %0.1f %s, speed = %0.2f m s^{-1}, T = %0.3f^{\\circ}C', ...
                            slice, P_string, PD(slice), unit_string, speed, T);
title_string_K{end+1} = [' Method = ', num2str(method'),';',...
                        '   FM = [ ', num2str(FM',2),' ];',...
                        '   (approx. dof_\epsilon )^{1/2} = [ ', num2str(round(sqrt(dof_e')),4),' ]'];
% title_string_K{end+1} = [];


%%
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

plot_wish_list = {
    % Plot settings for all channels.
    % Name          Legend Title                        Colour  Size
    {'P_sh1_clean', '\nablau_1 clean', 'b',    3   }
    {'P_sh1',       '\nablau_1',       'b',    1.5 }
    {'phi1',        e1_string,                          'k',    1.5 }
    {'P_sh2_clean', '\nablau_2 clean', 'r',    3   }
    {'P_sh2',       '\nablau_2',       'r',    1.5 }
    {'phi2',        e2_string,                          'k',    1.5 }
    {'P_sh3_clean', '\nablau_3 clean', 'g',    3   }
    {'P_sh3',       '\nablau_3',       'g',    1.5 }
    {'phi3',        e3_string,                          'k',    1.5 }
    {'P_sh4_clean', '\nablau_4 clean', 'm',    3   }
    {'P_sh4',       '\nablau_4',       'm',    1.5 }
    {'phi4',        e4_string,                          'k',    1.5 }};

if we_have_scalar_spectra
    % We need to do this because vector names cannot have indicies.
    for i = 1:size(P_scalar,2)
        junk_names = scalar_names{i};
        
        if length(scalar_names{i}) > 1
            scalar_names{i} = ['\nabla' junk_names(end-1) '_' junk_names(end)];
        else
            scalar_names{i} = ['\nabla' junk_names(end)];
        end
        
        temp = P_scalar(:,i);
        %        eval(['d.diss.P_scalar_' num2str(i) ' = temp;']);
        eval(['P_scalar_' num2str(i) ' = temp;']);
        
        %        plot_wish_list{end+1} = ...
        %           {['d.diss.P_scalar_' num2str(i)], scalar_names{i}, my_colour(:,i), 2.0};
        plot_wish_list{end+1} = ...
            {['P_scalar_' num2str(i)], scalar_names{i}, my_colour(:,i), 2.0};
    end
end

point_wish_list = {
    % X pos         Y pos                           Label                                           Colour  Size
    {'K_max(1)',    'P_sh1_clean(K_max_index_1)',   '[''k_{max} u_1='' num2str(round(K_max(1))) ''cpm'']',  'b',    18  }
    {'K_max(2)',    'P_sh2_clean(K_max_index_2)',   '[''k_{max} u_2='' num2str(round(K_max(2))) ''cpm'']',  'r',    18  }
    {'K_max(3)',    'P_sh3_clean(K_max_index_3)',   '[''k_{max} u_3='' num2str(round(K_max(3))) ''cpm'']',  'g',    18  }
    {'K_max(4)',    'P_sh4_clean(K_max_index_4)',   '[''k_{max} u_4='' num2str(round(K_max(4))) ''cpm'']',  'm',    18  }};

% Remove all channels that are not available from the wish list.
plot_list = {};
for ch = plot_wish_list'
    if exist(ch{1}{1}, 'var') plot_list{end+1,1} = ch{1}; end
end

% Remove all points that are not available from the wish list.  Evaluate
% and replace the values.
point_list = {};
for pt = point_wish_list'
    try 
        point_list{end+1,1} = {eval(pt{1}{1}), eval(pt{1}{2}), eval(pt{1}{3}), pt{1}{4}, pt{1}{5}};
    catch continue; end
end

% Extract the required data then plot on a loglog plot.
Y_data = [];
for ch = plot_list', Y_data(:,end+1 : end + size(eval(ch{1}{1}),2)) = eval( ch{1}{1} ); end

K      = K(2:end,:); %Remove zero wavenumber for log-log plots
Y_data = Y_data(2:end,:);

log_plot = loglog(ax, K, Y_data);

% Set the requested line width and colour settings.
for ch = 1:size(plot_list,1)
    set(log_plot(ch),'linewidth',plot_list{ch}{4},'color',plot_list{ch}{3});
end

% Plot the requested points with correct colours.
hold(ax, 'on');
for pt = point_list'
    loglog( ax, pt{1}{1}, pt{1}{2},  ...
        'Marker',           '^',            ...
        'MarkerSize',       pt{1}{5},       ...
        'MarkerFaceColor',  pt{1}{4},       ...
        'MarkerEdgeColor',  'w');
end
hold( ax, 'off' );
grid( ax, 'on' );

set(ax, 'ylim',y_lim, 'xlim', [K(1) K(end)])

xlabel(ax, '\itk \rm [cpm]')
ylabel(ax, '\Phi (\itk\rm)  [s^{-2} cpm^{-1}]')

legend_list = {};
for plt = plot_list',  legend_list{end+1} = plt{1}{2}; end
for plt = point_list', legend_list{end+1} = plt{1}{3}; end
legend(ax, legend_list, 'location','NorthEastOutside');

title(ax, title_string_K)
%set(gca,'fontsize',14)

end

