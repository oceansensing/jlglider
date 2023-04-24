%% show_spec
% Plot frequency and wavenumber spectra of shear, scalar gradients and
% acceleration.
%%
% <latex>\index{Functions!show\_spec}</latex>
%
%%% Syntax
%   show_spec( diss, titleString )
%
% * [diss] Structure of dissipation and other information returned by,
%       for example, the functions quick_look, get_diss_odas, or
%       get_scalar_spectra_odas.
% * [titleString] Optional string to prepend to each figure title. Use
%       it to help identify the current profile.
%
%%% Description
% Function to plot the spectra of shear, scalar gradients and acceleration
% that are embedded in the structure $\texttt{diss}$.
%
% The structure $\texttt{diss}$ contains sets of spectra, one set for every
% dissipation estimate within the structure. This GUI plots all spectra
% (shear, Nasmyth, scalar-gradients, and acceleration)
% from a chosen set, and lets you tranverse among the sets of spectra using
% a slider bar, or the arrow keys. The spectra are displayed in a single
% figure, with frequency spectra on the right, and wavenumber specra on the
% left. There are options to export a single set into a pdf-file. Optimal
% settings for PDF export are pre-selected and can be changed.
%
%%% Quick Overview
% The $\texttt{show\_spec}$ function visualizes spectra that was generated
% using a function such as $\texttt{quick\_look}$. As such, before
% show_spec can be run the spectra must be generated.
%
%      >> diss = quick_look('DAT_001', 10, 20);
%
% With the $\texttt{diss}$ structure generated, $\texttt{show\_spec}$ can
% be launched. 
%
%      >> show_spec( diss )
%
% show_spec will now display the main window. The window displays the
% spectra from one segment used to estimate the rate of dissipation,
% $\epsilon$. The segment is identified by its index number in the title,
% along with the segment-average pressure, speed and temperature. Use the
% arrow keys or scrollbar, to navigate among all segments in the
% $\texttt{diss}$ structure. 
% 
% @image @images/plot_spec1 @View of the $\texttt{show\_spec}$ main window.
% @View of the $\texttt{show\_spec}$ main window. You navigate between
% segments by either using the arrow keys, changing the scroll bar at the
% top of the window, or by entering a specific index into the text field
% located in the top right corner of the window. Data courtesy of Ilker Fer,
% University of Bergen.
%
% Spectra of interest are exported using the ``Export'' button found at the
% top of the screen. This button opens a dialog box that simplifies the
% saving of a plot into a file.
%
% Plots are generated using the third-party package $\texttt{export\_fig}$.
% This package bypasses the MATLAB PDF rendering engine to use a locally
% installed version of Ghostscript to generate PDF files that are
% attractive and well suited for publication. 
%
% You must manually download and install $\texttt{export\_fig}$ before
% $\texttt{show\_spec}$ can export images. The package is free, open
% source, and available from the MATLAB file exchange website.
% $\texttt{export\_fig}$ requires a local version of Ghostscript. A dialog
% box with instructions is displayed if a local version of Ghostscript is
% not found. Users of LaTeX likely have Ghostscript already installed.
%
% @image @images/plot_spec2 @Exporting images with $\texttt{show\_spec}$.
% @Exporting images is simplified by using the export dialog - activated by
% the "Export" button found on the main window.
%
%%% Notes
%
% # Default font sizes are used when generating plots. On some computers
%   the default values might be too small. To change the default font size
%   for all MATLAB plots, put the following into your startup.m
%   configuration file.
%
%      >> set(0,'DefaultAxesFontSize',14)
%
% # This function requires Matlab version 2014b and higher.
%

% *Version History:*
%
% 2014-07-18 WD Original version
% 2015-04-20 WD Modified to make compatible with latest versions of Matlab.
%               Provided warning when PDF export library not found.
% 2015-05-28 WD Fixed bugs with latest versions of Matlab
% 2015-07-27 WD Modern Matlab rendering fix.
% 2015-07-29 WD Input arguments modified.
% 2015-07-30 WD Rename to plot_spec
% 2015-10-30 RGL Update documentation
% 2015-11-23 RGL Renamed function to show_spec. Renamed show_freq_spec and
%       show_wave_spect ot plot_xxx_spec.
% 2020-06-17 JMM Modified the size of the figure to be about 80% of the
%                screen size. Should be display independent and be centered
%                on the screen.

function show_spec(diss, titleString)

if nargin == 0
    error('Input dissipation structure required.')
elseif nargin == 1
    titleString = [];
end

% Check for old versions of MATLAB
for v = ver
    if strcmp(v(1).Name, 'MATLAB') && str2double(v(1).Version) < 8.4
        error('MATLAB R2014b or higher required for this function.');
    end
end

figSize = [0 0 0 0];

% Determine a good size for the figure
f0 = figure('Position',get(0,'screensize'),'Visible','off');
set(f0,'units','characters')
size_characters = get(gcf,'position');
width = size_characters(3);
height = size_characters(4);
close(f0)

size_characters = [0.1*width 0.1*height 0.8*width 0.8*height];
f = figure( 'Units','characters',...
            'Position',size_characters,...
            'Renderer','painters',...
            'PaperPositionMode', 'auto', ...
            'Name', 'Wavenumber and Frequency Spectra', ...
            'KeyPressFcn',@figKeyPressed, ...
            'Visible', 'off', ...
            'ResizeFcn',@resizeMaster);
            
sliderPanel = uipanel('Units','characters','Parent',f,...
                      'BorderType', 'none', ...
                      'BackgroundColor', 'w', ...
                      'ResizeFcn',@resizeSliderPanel);
leftPanel   = uipanel('Units','characters','Parent',f,...
                      'BorderType', 'none', ...
                      'BackgroundColor', 'w', ...
                      'ResizeFcn',@resizeLeftPanel);
rightPanel  = uipanel('Units','characters','Parent',f,...
                      'BorderType', 'none', ...
                      'BackgroundColor', 'w', ...
                      'ResizeFcn',@resizeRightPanel);
    
haxl = axes('Parent',leftPanel,'Units','characters');
haxr = axes('Parent',rightPanel,'Units','characters');

% Assign data to the graph - to be processed when plotting a slice.
data.diss = diss;
data.title_string = titleString;
data.nargin       = nargin;
data.slice        = 1;
    
guidata(f, data);

slices = size(diss.e, 2);
try step = 1 / (slices - 1);
catch, step = 1;
end
    
slider = uicontrol('Style', 'slider',...
    'Min',1,'Max',slices,'Value',data.slice,'SliderStep', [step,step], ...
    'Tag', 'SliceSlider', ...
    'Parent', sliderPanel, ...
    'Units', 'normalized', ...
    'Callback', {@untitled_SliderChange,haxl});
    
sliceText = uicontrol('Style', 'text', 'String', 'Index:  1', ...
    'BackgroundColor', 'w',...
    'Units','characters','Parent',sliderPanel);

exportButton = uicontrol('Style', 'pushbutton', 'String', 'Export', ...
    'Parent', sliderPanel, ...
    'Units', 'characters', ...
    'Callback', @untitled_ExportClicked);

% This function causes the figure to plot the first slice.
untitled_SliderChange(slider, '', haxl);

set(f, 'Visible', 'on');

% --- Executes just before untitled is made visible.
function untitled_SliderChange(hObject, event, ax)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to untitled (see VARARGIN)

    val = round(get(hObject,'Value'));

    % Update handles structure
    data = guidata(hObject);

    set(sliceText,'String',sprintf('Index: %2i',val));
    resizeSliderPanel();

    data.slice = val;
    guidata(hObject, data);

    plotData(hObject, haxl);

    if isfield(data, 'gui')
        if isfield( data.gui, 'refresh' )
            fun = data.gui.refresh;
            fun(hObject, event);
        end
    end
end


function plotData( hObject, ax )
    d = guidata(f);
    plot_wave_spec(haxl, d.slice, d.diss, d.title_string);
    plot_freq_spec(haxr, d.slice, d.diss, d.title_string);
end


function figKeyPressed( hObject, event )
%    c = event.Character;
    d = guidata(hObject);
    value = get( slider, 'Value' );
    oldValue = value;
    
    if strcmp('leftarrow', event.Key) || ...
       strcmp('uparrow', event.Key) || ...
       strcmp('backspace', event.Key),
        value = max( value-1, 1 );
    elseif strcmp('rightarrow', event.Key) || ...
       strcmp('downarrow', event.Key) || ...
       strcmp('space', event.Key),
        value = min( value+1, size(d.diss.e, 2) );
    end
    
    if oldValue == value, return; end

    set( slider, 'Value', value );
    untitled_SliderChange(slider, event, haxl);
end


function resizeMaster(src, evt)
    fpos = get(f,'Position');
    sh = 2;

    c = get(leftPanel, 'Children');
    cidx = strcmp(get(c, 'Tag'), 'legend');
    chPosl = get(c(cidx), 'Position');
    if isempty(chPosl), chPosl = [0,0,0,0]; end
    
    c = get(rightPanel, 'Children');
    cidx = strcmp(get(c, 'Tag'), 'legend');
    chPosr = get(c(cidx), 'Position');
    if isempty(chPosr), chPosr = [0,0,0,0]; end
    
    diff = (chPosl(3) - chPosr(3)) / 2;
    
    set(sliderPanel, 'Position', [0,fpos(4)-sh,fpos(3),sh]);
    set(leftPanel,  'Position', [0,0,fpos(3)/2+diff,fpos(4)-sh]);
    set(rightPanel, 'Position', [fpos(3)/2+diff,0,fpos(3)/2-diff,fpos(4)-sh]);
    
    d = guidata(f);
    if isfield(d, 'gui')
        if isfield( d.gui, 'refresh' )
            fun = d.gui.refresh;
            fun(src, evt);
        end
    end
end

function resizeSliderPanel(src, evt)
    pPos = get(sliderPanel, 'Position');
    sPos = get(sliceText, 'Extent');
    bPos = get(exportButton, 'Position');
    textWidthRatio = (sPos(3) + bPos(3) + 2) / pPos(3);
    set(slider, 'Position', [0,0,1-textWidthRatio,1]);
    set(sliceText, 'Position', [pPos(3) - sPos(3) - bPos(3) - 2, 0, sPos(3), sPos(4)]);
    set(exportButton, 'Position', [pPos(3) - bPos(3) - 1, 0, bPos(3), bPos(4)]);
end

function resizeLeftPanel(src, evt)
    cPos = get(leftPanel, 'Position');
    iPos = get(haxl, 'TightInset');
    
    c = get(leftPanel, 'Children');
    cidx = strcmp(get(c, 'Tag'), 'legend');
    chPos = get(c(cidx), 'Position');
    if isempty(chPos), chPos = [0,0,0,0]; end

    chPos = [chPos(1)*cPos(3), chPos(2)*cPos(4), chPos(3)*cPos(3), chPos(4)*cPos(4)];

    spacer = chPos(3) + 2.5;
    
    set(haxl, 'Position', [iPos(1), iPos(2), cPos(3)-iPos(1)-iPos(3)-spacer, cPos(4)-iPos(2)-iPos(4)-1]);
end

function resizeRightPanel(src, evt)
    cPos = get(rightPanel, 'Position');
    iPos = get(haxr, 'TightInset');

    c = get(rightPanel, 'Children');
    cidx = strcmp(get(c, 'Tag'), 'legend');
    chPos = get(c(cidx), 'Position');
    if isempty(chPos), chPos = [0,0,0,0]; end
    
    chPos = [chPos(1)*cPos(3), chPos(2)*cPos(4), chPos(3)*cPos(3), chPos(4)*cPos(4)];
    
    spacer = chPos(3) + 2;

    set(haxr, 'Position', [iPos(1), iPos(2), cPos(3)-iPos(1)-iPos(3)-spacer, cPos(4)-iPos(2)-iPos(4)-1]);
end


function untitled_ExportClicked(hObject,eventdata)
    nominalSizes = [0 0 467 389; 0 0 414 389; 0 0 881 389];

    plotLFig = true;
    plotRFig = true;

    d = guidata(f);
    if isfield(d, 'gui')
        if isfield(d.gui, 'currentFigID')
            try
                figure(d.gui.currentFigID);
                return;
            catch,
            end
        end
    end

    dBox = dialog('WindowStyle', 'normal', 'Name', 'Export Images', ...
            'Resize', 'on', 'Units', 'characters', ...
            'Visible', 'off', ...
            'ResizeFcn',@resizeHPanel);
        
    dBoxP = get(dBox, 'Position');
    set(dBox, 'Position', [dBoxP(1), dBoxP(2), dBoxP(3)+21, dBoxP(4)+8]);

    pan1 = uipanel('Units','characters','Parent',dBox,...
            'Title','Export Plot');
    pan2 = uipanel('Units','characters','Parent',dBox,...
            'Title','Image Type');
    pan5 = uipanel('Units','characters','Parent',dBox,...
            'Title','Image Size');
    pan3 = uipanel('Units','characters','Parent',dBox,...
            'Title','Output File');
    pan4 = uipanel('Units','characters','Parent',dBox,...
            'BorderType', 'none');
    
    h1 = uibuttongroup('Position', [0,0,1,1], 'Parent', pan1, ...
            'BorderType', 'none');
%    h2 = uibuttongroup('Position', [0,0,1,1], 'Parent', pan2, ...
%            'BorderType', 'none');
    h3 = uibuttongroup('Position', [0,0,1,1], 'Parent', pan3, ...
            'BorderType', 'none');
    h5 = uibuttongroup('Position', [0,0,1,1], 'Parent', pan5, ...
            'BorderType', 'none');
    
    u0 = uicontrol('Style','Radio','String',' Wavenumber Spectrum', ...
            'Units','characters', 'Tag', '1', ...
            'pos',[3 5 50 2],'parent',h1,'HandleVisibility','off');
    u1 = uicontrol('Style','Radio','String',' Frequency Spectrum', ...
            'Units','characters', 'Tag', '2', ...
            'pos',[3 3 50 2],'parent',h1,'HandleVisibility','off');
    u2 = uicontrol('Style','Radio','String',' Wavenumber + Frequency Spectrum', ...
            'Units','characters', 'Tag', '3',...
            'pos',[3 1 50 2],'parent',h1,'HandleVisibility','off');

    u3 = uicontrol('Style','checkbox','String',' EPS', ...
            'Units','characters', ...
            'pos',[3 5 35 2],'parent',pan2,'HandleVisibility','off');
    u4 = uicontrol('Style','checkbox','String',' PDF (Requires Ghostscript)', ...
            'Units','characters', ...
            'pos',[3 3 35 2],'parent',pan2,'HandleVisibility','off');
    u5 = uicontrol('Style','checkbox','String',' PNG (Requires Ghostscript)', ...
            'Units','characters',...
            'pos',[3 1 35 2],'parent',pan2,'HandleVisibility','off');

    u6 = uicontrol('Style','Radio','String',' Use default file name and path', ...
            'Units','characters', 'Tag', '1', ...
            'pos',[3 5 50 2],'parent',h3,'HandleVisibility','off');
    u7 = uicontrol('Style','Radio','String',' Specify full path:', ...
            'Units','characters', ...
            'pos',[3 3 50 2],'parent',h3,'HandleVisibility','off');


    u8 = uicontrol('Style','Radio','String',' Current Figure Size', ...
            'Units','characters', 'Tag', '1', ...
            'pos',[3 5 35 2],'parent',h5,'HandleVisibility','off');
    u9 = uicontrol('Style','Radio','String',' Nominal', ...
            'Units','characters', 'Tag', '2', ...
            'pos',[3 3 35 2],'parent',h5,'HandleVisibility','off');
    u10 = uicontrol('Style','Radio','String',' Custom:', ...
            'Units','characters', 'Tag', '3',...
            'pos',[3 1 14 2],'parent',h5,'HandleVisibility','off');
    xDim = uicontrol('Style','edit','String','6', ...
            'Units','characters','parent',h5, ...
            'Callback', @textF_selectSize, ...
            'pos',[18 1 8 2]);
    xTxt = uicontrol('Style','text','String','x', ...
            'Units','characters','parent',h5, ...
            'pos',[26.5 1.5 2 1]);
    yDim = uicontrol('Style','edit','String','4', ...
            'Units','characters','parent',h5, ...
            'Callback', @textF_selectSize, ...
            'pos',[29 1 8 2]);
    units = uicontrol('Style','popup','Units','characters','parent',h5, ...
            'pos',[38 1 13 2], ...
            'Callback', @radioB_selectSize, ...
            'String', 'inch|cm|mm|point');
        
    
    defaultPath = sprintf('%s%s%s', cd, filesep, 'defaultFile');
    pathTF = uicontrol('Style', 'edit', 'String', defaultPath, ...
            'Units','characters', 'HorizontalAlignment', 'left', ...
            'pos',[6 1 100 2],'parent',h3,'Enable','off');
        
    cancelButton = uicontrol('Style', 'pushbutton', 'String', 'Cancel', ...
        'Parent', pan4, ...
        'Units', 'characters', ...
        'Callback', @closeWindow);

    executeButton = uicontrol('Style', 'pushbutton', 'String', 'Execute', ...
        'Parent', pan4, ...
        'Units', 'characters', ...
        'Callback', @exportPlot);

    set(h1, 'SelectionChangeFcn', @radioB_selectPlot);
    set(h1, 'SelectedObject', u2);

%    set(h2, 'SelectedObject', u4);

    set(h3, 'SelectionChangeFcn', @radioB_fileName);
    set(h3, 'SelectedObject', u6);
    
    set(h5, 'SelectionChangeFcn', @radioB_selectPlot);
    set(h5, 'SelectedObject', u8);
    
    if isfield(d, 'gui') && isfield(d.gui, 'u0')
        set(u0, 'Value', d.gui.u0);
        set(u1, 'Value', d.gui.u1);
        set(u2, 'Value', d.gui.u2);
        set(u3, 'Value', d.gui.u3);
        set(u4, 'Value', d.gui.u4);
        set(u5, 'Value', d.gui.u5);
        set(u6, 'Value', d.gui.u6);
        set(u7, 'Value', d.gui.u7);
        set(u8, 'Value', d.gui.u8);
        set(u9, 'Value', d.gui.u9);
        set(u10, 'Value', d.gui.u10);
        set(xDim, 'String', d.gui.xDim);
        set(yDim, 'String', d.gui.yDim);
        set(units, 'Value', d.gui.units);
        set(pathTF, 'String', d.gui.pathTF);
        set(pathTF, 'Enable', d.gui.pathTF_en);
    end
    
    radioB_fileName(h3,[]);
    radioB_selectPlot(h5,[]);
    
    d.gui.refresh = @radioB_selectPlot;
    d.gui.currentFigID = dBox;
    guidata(f, d);
    
    set(dBox, 'Visible', 'on');

    
    function resizeHPanel(src, evt)
        bPos = get(dBox, 'Position');
        set(pan1, 'Position', [2,bPos(4)-11,bPos(3)-4,9]);
        set(pan2, 'Position', [2,bPos(4)-21,(bPos(3)-4)/2-1,9]);
        set(pan5, 'Position', [bPos(3)/2+1,bPos(4)-21,(bPos(3)-4)/2-1,9]);
        set(pan3, 'Position', [2,bPos(4)-31,bPos(3)-4,9]);
        set(pathTF, 'Position', [7 1 bPos(3)-13 2]);
        set(pan4, 'Position', [2,bPos(4)-35,bPos(3)-4,4]);
        set(cancelButton, 'Position', [bPos(3)-31 1 12 2]);
        set(executeButton, 'Position', [bPos(3)-17 1 12 2]);
    end
    
    function exportPlot(src, evt)
        str = get(units, 'String');
        unit_idx = get(units, 'Value');
        fig_units = strtrim(str(unit_idx,:));
        
        d = guidata(f);

        x = str2double(get(xDim, 'String'));
        y = str2double(get(yDim, 'String'));
        
        f1 = figure( 'Units',fig_units,'Position',[0 0 x y],...
            'Renderer','painters',...
            'ResizeFcn',@resizePrintMaster, ...
            'Visible', 'off', ...
            'PaperPositionMode', 'auto');
        
        set(f1, 'Units', 'characters');
        
        if plotLFig
            printPanelL = uipanel('Units','characters','Parent',f1,...
                                  'BorderType', 'none', ...
                                  'BackgroundColor', 'w', ...
                                  'ResizeFcn',@resizePrintPanelL);
            axL = axes('Parent',printPanelL,'Units','characters');
        end
        
        if plotRFig
            printPanelR = uipanel('Units','characters','Parent',f1,...
                          'BorderType', 'none', ...
                          'BackgroundColor', 'w', ...
                          'ResizeFcn',@resizePrintPanelR);
            axR = axes('Parent',printPanelR,'Units','characters');
        end    

        if plotLFig
            plot_wave_spec(axL, d.slice, d.diss, d.title_string);
            xlim = get(haxl, 'XLim');
            ylim = get(haxl, 'YLim');
            set(axL, 'XLim', xlim);
            set(axL, 'YLim', ylim);
        end
        if plotRFig
            plot_freq_spec(axR, d.slice, d.diss, d.title_string);
            xlim = get(haxr, 'XLim');
            ylim = get(haxr, 'YLim');
            set(axR, 'XLim', xlim);
            set(axR, 'YLim', ylim);
        end
        
        set(f1, 'Visible', 'on');
        resizePrintMaster();
        
        filePath = get(pathTF, 'String');
        
        eps_en = get(u3,'Value');
        pdf_en = get(u4,'Value');
        png_en = get(u5,'Value');
        
        % Issue warning if the export figure functions are not in the
        % current path.
        fPath = which('export_fig.m');
        if isempty(fPath)
            warning('Can not find the "export_fig" functions.  Export disabled.');
            disp('Note: Download and add the library to your path to remove');
            disp('      this warning.  Library is free and can be obtained from');
            disp('      the Mathworks FileExchange website.');
            return
        end
        
        if eps_en
            export_fig(filePath, '-eps', '-transparent', '-painters');
        end
        
        if pdf_en
            export_fig(filePath, '-pdf', '-transparent', '-painters');
        end

        if png_en
            export_fig(filePath, '-png', '-m2', '-painters');
        end
        
        closeWindow(src, evt);
        close(f1);

        function resizePrintMaster(src, evt)
            chPosL = [0 0 0 0];
            chPosR = [0 0 0 0];
            
            fpos = get(f1,'Position');
            
            if plotLFig
                c = get(printPanelL, 'Children');
                cidx = strcmp(get(c, 'Tag'), 'legend');
                for i = cidx'
                    if i
                        chPosL = get(c(cidx(i)), 'Position');
                    end
                end
            end
            
            if plotRFig
                c = get(printPanelR, 'Children');
                cidx = strcmp(get(c, 'Tag'), 'legend');
                for i = cidx'
                    if i
                        chPosR = get(c(cidx(i)), 'Position');
                    end
                end
            end
            
            diff = (chPosL(3) - chPosR(3)) / 2;

            lWidth = fpos(3)/2+diff;
            rWidth = fpos(3)/2-diff;
            rStart = lWidth;
            
            if ~plotLFig || ~plotRFig
                lWidth = fpos(3);
                rWidth = fpos(3);
                rStart = 0;
            end
            
            if plotLFig
                set(printPanelL, 'Position', [0,0,lWidth,fpos(4)]);
            end
            
            if plotRFig
                set(printPanelR, 'Position', [rStart,0,rWidth,fpos(4)]);
            end
        end
        
        function resizePrintPanelL(src, evt)
            if ~plotLFig, return; end

            chPos = [0 0 0 0];
            cPos = get(printPanelL, 'Position');
            iPos = get(axL, 'TightInset');
    
            c = get(printPanelL, 'Children');
            cidx = strcmp(get(c, 'Tag'), 'legend');
            for i = cidx'
                if i
                    chPos = get(c(cidx(i)), 'Position');
                end
            end
            chPos = [chPos(1)*cPos(3), chPos(2)*cPos(4), chPos(3)*cPos(3), chPos(4)*cPos(4)];                       
            spacer = chPos(3) + 2.5;
    
            set(axL, 'Position', [iPos(1), iPos(2), cPos(3)-iPos(1)-iPos(3)-spacer, cPos(4)-iPos(2)-iPos(4)-1]);
        end
        
        function resizePrintPanelR(src, evt)
            if ~plotRFig, return; end

            chPos = [0 0 0 0];
            cPos = get(printPanelR, 'Position');
            iPos = get(axR, 'TightInset');

            c = get(printPanelR, 'Children');
            cidx = strcmp(get(c, 'Tag'), 'legend');
            for i = cidx'
                if i
                    chPos = get(c(cidx(i)), 'Position');
                end
            end
            chPos = [chPos(1)*cPos(3), chPos(2)*cPos(4), chPos(3)*cPos(3), chPos(4)*cPos(4)];
            spacer = chPos(3) + 2.5;

            set(axR, 'Position', [iPos(1), iPos(2), cPos(3)-iPos(1)-iPos(3)-spacer, cPos(4)-iPos(2)-iPos(4)-1]);
        end
    end
    
    function closeWindow(src, evt)
        d = guidata(f);
        d.gui.u0 = get(u0, 'Value');
        d.gui.u1 = get(u1, 'Value');
        d.gui.u2 = get(u2, 'Value');
        d.gui.u3 = get(u3, 'Value');
        d.gui.u4 = get(u4, 'Value');
        d.gui.u5 = get(u5, 'Value');
        d.gui.u6 = get(u6, 'Value');
        d.gui.u7 = get(u7, 'Value');
        d.gui.u8 = get(u8, 'Value');
        d.gui.u9 = get(u9, 'Value');
        d.gui.u10 = get(u10, 'Value');
        d.gui.xDim = get(xDim, 'String');
        d.gui.yDim = get(yDim, 'String');
        d.gui.units = get(units, 'Value');
        d.gui.pathTF = get(pathTF, 'String');
        d.gui.pathTF_en = get(pathTF, 'Enable');        
        d.gui = rmfield(d.gui, 'refresh');
        d.gui = rmfield(d.gui, 'currentFigID');
        guidata(f, d);
        close(dBox);
    end
    
    function radioB_fileName(src, evt)
        tag = get(get(src,'SelectedObject'),'Tag');
        
        if strcmp(tag, '1')
            % Generate a suitable file name.  Just add the slice # until 
            % Rolf has a chance to review.
            d = guidata(f);
            if isfield(d, 'graphic_prefix')
                prefix = d.graphic_prefix;
            else
                prefix = 'data_';
            end
            fileName = [prefix num2str(d.slice)];
            
            defaultPath = sprintf('%s%s%s', cd, filesep, fileName);
            set(pathTF, 'String', defaultPath);
            set(pathTF, 'Enable', 'off');
        else
            set(pathTF, 'Enable', 'on');
        end
    end
    
    function radioB_selectPlot(src, evt)
        leftPanelUnits = get(leftPanel, 'Units');
        set(leftPanel, 'Units', 'points');
        rightPanelUnits = get(rightPanel, 'Units');
        set(rightPanel, 'Units', 'points');
        
        figSizeL = get(leftPanel, 'Position');
        figSizeR = get(rightPanel, 'Position');
        figSizeAll = [0,0,figSizeL(3)+figSizeR(3),figSizeR(4)];
        
        tag = get(get(h1,'SelectedObject'),'Tag');

        if strcmp(tag, '1')
            figSize = figSizeL;
            plotLFig = true;
            plotRFig = false;
        elseif strcmp(tag, '2')
            figSize = figSizeR;
            plotLFig = false;
            plotRFig = true;
        else
            figSize = figSizeAll;    
            plotLFig = true;
            plotRFig = true;
        end

        set(leftPanel, 'Units', leftPanelUnits);
        set(rightPanel, 'Units', rightPanelUnits);

        radioB_fileName(h3, evt);
        radioB_selectSize(src, evt);
    end

    function radioB_selectSize(src, evt)
        dispSize = figSize;
        unit_idx = get(units, 'Value');
        switch unit_idx
            case 1
                scale = 1/72;
            case 2
                scale = 2.54 * 1/72;
            case 3
                scale = 25.4 * 1/72;
            otherwise
                scale = 1;
        end
        
        if get(u9, 'Value')
            tag = str2double(get(get(h1,'SelectedObject'),'Tag'));
            dispSize = nominalSizes(tag,:);
        end

        if get(u10, 'Value')
            dispSize(3) = str2double(get(xDim, 'String')) / scale;
            dispSize(4) = str2double(get(yDim, 'String')) / scale;
        end

        set(xDim, 'String', sprintf('%.1f',dispSize(3)*scale));
        set(yDim, 'String', sprintf('%.1f',dispSize(4)*scale));
    end

    function textF_selectSize(src, evt)
        set(u10, 'Value', 1);
    end
end

end
