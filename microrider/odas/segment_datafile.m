%% segment_datafile
% Break data file into segments around bad buffers
%%
% <latex>\index{Functions!segment\_datafile}</latex>
% 
%%% Syntax
%   segment_datafile( file_name, min_record_length )
%
% * [file_name] Name of file with or without extension. Wildcard characters
%         accepted in the same format as the dir() function.
% * [min_record_length] Minimum size of a valid data segment. Default = 20.
%
%%% Description
% Segment a raw RSI data file with BAD BUFFERS into smaller data files 
% without BAD BUFFERS.
%
% Interpolation of bad buffer data works well for patching over short bad
% buffer segments within a data file. This is the default behaviour
% performed by the quick_look function. When bad buffers persist for
% multiple records, simply patching the bad buffers with interpolated
% values can induce notable errors in data channels with pre-emphasis.
%
% This function removes bad buffer data by segmenting the original data
% file into segments where each segment contains no bad buffers. This
% prevents interpolated values from being used when deconvolving
% pre-emphasised channel data.
%
% Each resulting data file has a minimum length of min_record_length 
% records (seconds). The new files are named by appending _XX to the
% original file name where _XX is a counter that increments for each
% observed segment. If the orginal file does not have a BAD BUFFER then the
% only file created is _00.  This allows subsequent functions and scripts
% to work regardless of the number of BAD BUFFERS in a file.

% 2015-12-31 RGL Original version
% 2016-01-04 WID Added input arguments and support for multiple input
%                files. Documentation update and function rename.
% 2016-01-15 WID Fixed ASCII bug with international versions of Matlab.

function segment_datafile( file_name, min_record_length )

if nargin < 1
    error('Invalid input. Please supply file name.');
end

if nargin < 2
    L = 20; % The minimum number of good records in a file
else
    L = min_record_length;
end

% The extension '.p' is optional.  Note that UNIX type systems will require
% an explicit test for '.p' in addition to '.P'.
if length(file_name) < 3 || ~strcmpi(file_name(end-1:end), '.p')
    file_name = [file_name '.P'];
end

file_list = struct2cell(dir(file_name));
file_list = file_list(1,:);

if size(file_list, 2) < 1
    error(['Warning;  File "' file_name '" does not exist']);
end

for file = file_list    
    [P,N,E] = fileparts(char(file));

     in_file_name = [N E];
    out_file_name = [N '_00' E];

    % check for bad buffers
    bad_buffers = check_bad_buffers(in_file_name);
    if isempty(bad_buffers) % no bad buffers in this file
        [SUCCESS,MESSAGE,MESSAGEID] = copyfile(in_file_name, out_file_name);
        if ~SUCCESS
            error('This directory may not be writable')
        end
    else
        segments = find_segments(bad_buffers, L);
        fid = fopen_odas(in_file_name, 'r');
        [header_length, record_length, cfg_length, total_records, ...
            endian_flag] = get_file_info(fid);
        
        % check how many record there are after the last bad record
        if total_records - segments(1,end) >= L
            segments(2,end) = total_records;
        else
            segments(:,end) = [];
        end
        
        % now read the input file
        fseek(fid, 0, 'bof'); % move to start of file
        record_hdr = fread(fid, header_length, 'short');
        record_cfg = fread(fid, cfg_length, 'char*1');
        records = fread(fid, record_length * total_records, 'short');
        fclose(fid);

        records = reshape(records, record_length, []); % every column is a record

        for index = 1:size(segments,2)
            out_file_name = sprintf('%s_%02d%s', N, index-1, E);
            fid = fopen(out_file_name, 'w', endian_flag, 'US-ASCII');
            fwrite(fid, record_hdr, 'short');
            fwrite(fid, record_cfg, 'char*1');
            fwrite(fid, records(:, segments(1,index):segments(2,index)), 'short');
            fclose(fid);
        end
    end
end


function segments = find_segments(bad_buffers, L)
% Function to identify the segments without bad buffers and a length
% greater than L records

segments = zeros(2, length(bad_buffers)+1);

bad_buffers = [0; bad_buffers(:)];
bad_buffer_difference = diff(bad_buffers);

index = 0;

while ~isempty(bad_buffer_difference)
    if (bad_buffer_difference(1) > L + 1)
        index = index + 1;
        segments(1, index) = bad_buffers(1) + 1;
        segments(2, index) = bad_buffers(2) - 1;
    end
    bad_buffers(1) = [];
    bad_buffer_difference = diff(bad_buffers);
end
% now there is only one remaining bad buffer
index = index + 1;
segments(1, index) = bad_buffers(1) + 1;
segments(2, index) = inf;
segments = segments(:,1:index); % trim off the zeros at the end


function [header_length, record_length, cfg_length, total_records, endian_flag] = ...
    get_file_info(fid)
fseek(fid, 0, 'bof'); % move to beginning of file
header        = fread(fid, 64, 'short');
header_length = header(18)/2; % in words
cfg_length    = header(12); % in bytes
record_length = header(19)/2; % in words
endian = header(64);
if endian == 1, endian_flag = 'l'; end
if endian == 2, endian_flag = 'b'; end
if endian == 0, endian_flag = 'l'; end

fseek(fid, 0, 'eof'); % move to end of files
total_length  = ftell(fid); % in bytes
total_length  = total_length - 2*header_length; % in bytes
total_records = total_length / 2; % in words
total_records = floor (total_records / record_length);
fseek(fid, 0, 'bof'); % return pointer to start of file







