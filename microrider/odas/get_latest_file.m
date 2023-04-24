%% get_latest_file
% Get the name of the latest RSI raw binary data file in current folder.
%%
% <latex>\index{Functions!get\_latest\_file}</latex>
%
%%% Syntax
%   latest_file = get_latest_file( match_file )
%
% * [match_file] Expression to filter matching files.  Defaults to "*" for
%          accepting all files.
% * []
% * [latest_file] Name of the latest RSI raw binary data file in the current 
%          directory. Empty if no ODAS binary data files exist.
%
%%% Description
% Retrieve the name of the latest RSI raw binary data file in the local
% directory. It is useful for near real-time data processing and plotting
% because such operations are frequently conducted on the latest data file
% recorded in the local directory.  
%

% Version History
% * 2004-06-06 (RGL) initial
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-04-23 (WID) Updated to be case insensitive.  Simplified loops.
% * 2012-07-19 (WID) simplified loop
% * 2015-01-29 (WID) added match_file for use with calibration reports
% * 2017-05-14 (WID) add option to ignore extensions - for non-P files

function latest_file = get_latest_file( match_file, ignore_ext )

if nargin < 2, ignore_ext = false; end

% Set initial value to zero so an empty string is returned if no files exist.
newest.datenum = 0;
newest.name = [];

% Set default match_file string to anything if not provided.
if nargin < 1
    match_file = '*';
end

% Append the match string 
if ignore_ext
    for file = dir(match_file)
        if file.datenum > newest.datenum, newest = file; end
    end
else
    for file = [dir([match_file '.p'])' dir([match_file '.P'])']
        if file.datenum > newest.datenum, newest = file; end
    end
end

latest_file = newest.name;
