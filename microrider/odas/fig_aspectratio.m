%% fig_aspectratio.m
% Set the aspect ratio of a figure.
%%
% <latex>\index{Functions!fig\_aspectratio}</latex>
%
%%% Syntax
%   pos = fig_aspectratio ( h, AR )
%
% * [h] Figure handle
% * [AR] Desired aspect ratio
% * []
% * [pos] Figure position on screen
%
%%% Description
% A simple function to control the aspect ratios of figures

% Version History
%
% * 2017-01-16 (JMM) Written to better control figure sizes


function pos = fig_aspectratio(h,AR)
% set a figure aspect ratio
% h = figure handle
% AR = aspect ratio (width/height)

pos_orig = get(h,'Position');
pos = [pos_orig(1) pos_orig(2) AR*pos_orig(4) pos_orig(4)];
set(h,'Position',pos);

% Check to make sure figure is on the screen
pos = get(gcf,'Position');
screensize = get( groot, 'Screensize' );
if pos(3)>screensize(3);    warning('Figure too wide'); end
if pos(4)>screensize(4);    warning('Figure too wide'); end
if pos(1)+pos(3)>screensize(3)
    pos = [0 pos(2) pos(3) pos(4)];
end
if pos(2)+pos(4)>screensize(4)
    pos = [0 pos(2) pos(3) pos(4)];
end

