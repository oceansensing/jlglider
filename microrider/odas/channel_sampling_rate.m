%% channel_sampling_rate
% Calculate sampling rate of the requested channel
%%
% <latex>\index{Functions!channel\_sampling\_rate}</latex>
%
%%% Syntax
%   rate = channel_sampling_rate( channel, configStr, fs_fast )
%
% * [channel]   Name of the requested channel.
% * [configStr] Configuration string associated with the data.
% * [fs_fast]   Actual sampling rate of the instrument in rows / second.
%               This value is extracted from the data file header and 
%               typically saved within the variable 'fs_fast'.
% * []
% * [rate]      Sampling rate of specified channel in [samples/second].
%
%%% Description
% The sampling rate defined as the number of times a channel is sampled
% every second.  This function calculates the sampling rate based on the
% number of times the channel occurs within the channel matrix and the
% actual sampling rate of the instrument.
%
%%% Examples
%
%    >> rate = channel_sampling_rate( 'P', setupfilestr, fs_fast )
%
% Returns the sampling rate of channel 'P'.  The configuration string found
% in 'setupfilestr' is used to find the number of occurances of 'P' within 
% the channel matrix.  This value, along with 'fs_fast', is used to 
% calculate the resulting channel sampling rate.

% *Version History:*
%
% * 2013-05-16 (WID) initial version
% * 2013-07-29 (WID) error when requested channel not in setupfilestr
%%

function rate = channel_sampling_rate( channel, setupfilestr, fs_fast )

    % Pre-parse the configuration file for quick access.
    cfg = setupstr(setupfilestr);

    % Load rows from the configuration file
    rows = setupstr(cfg, 'matrix', 'row[0-9]+');
    ch_matrix = [];

    % Construct the channel matrix from the rows in the configuration file
    for row = rows
      values = textscan(row{1}, '%d');
      ch_matrix = vertcat(ch_matrix, values{1}');
    end
    
    % Find the channel id - needed when searching the matrix
    id = str2double(setupstr(cfg, channel, '(id)|(id_even)'));
    
    if isempty(id),
        error('Requested channel "%s" not in configuration string.', ...
            channel);
    end
    
    % Count the number of times the channel occurs in the matrix
    samples_per_matrix = length(find(ch_matrix==id(1)));
    
    % Sampling rate is calculated
    rate = fs_fast / length(rows) * samples_per_matrix;

end
