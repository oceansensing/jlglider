%% read_odas
% Convert a RSI raw binary data file into a raw mat-file
%%
% <latex>\index{Functions!read\_odas}</latex>
%
%%% Syntax
%   [ch_list, outstruct] = read_odas( fname )
%
% * [fname] Name of the RSI raw binary data file. Optional, you are
%      prompted, if omitted. 
% * []
% * [ch_list] List of channels found within the resulting mat-file.
% * [outstruct] Structure holding contents of the mat-file. If this output
%               argument is used, the mat-file is not created.
%
%%% Description
% Convert a RSI raw binary data file (.p) into a Matlab mat-file, without
% conversion into physical units. The raw data file is read, demultiplexed,
% and assembled into Matlab data vectors. All channel data, along with
% other relavent information such as the configuration string, is converted
% into a Matlab-readable format.
%
% If created, the mat-file is given the same name as the input file with
% the $\texttt{.p}$ extension replaced by $\texttt{.mat}$.  If this file
% already exists, you are queried for an alternate name.
%
% When the $\texttt{outstruct}$ variable is used, a mat-file is not generated.
% Instead, all values that would have been placed into the mat-file are
% returned as fields in the structure $\texttt{outstruct}$. 
%
% The vectors are not converted into physical units. However, 1e-12 is
% added to every vector so that Matlab is forced to save them as true
% double-precision floating-point numbers.
%
%%% Examples
%
%    >> read_odas
%
% Extract all data from a $\texttt{.p}$ file and store inside a newly
% generated mat-file. Input file requested from the user with the most
% recent data file being the default value.
%
%    >> my_vars = read_odas('my_file.p')
%
% Extract data from $\texttt{my\_file.p}$ and store them in the file
% $\texttt{my\_file.mat}$. A list of the extracted channels is in
% $\texttt{my\_vars}$. The $\texttt{.p}$ extension is optional.
%
%    >> [myVar, d] = read_odas( '*03' )
%
% Read the first raw data file that has a name ending in $\texttt{03}$. Do
% not create a mat-file. All data are returned within the structure named
% $\texttt{d}$.

% *Version History:*
%
% * 2005-02-01 (IG) Based in part on misc. OTL/RGL routines, including
%   plot_tomi, demultiplex, read_tomi, convert_header, convert_mat_slow,
%   convert_mat_fast, and others
% * 2005-05-02 (IG) added second C_dC channel to default channel number/name
%   setup to deal with Peter Winsor's instrument
% * 2006-04-17 (RGL) added 1e-12 to data so that all variables are saved
%   as true 8-byte floats. Also forced a -v6 save to mat-file
% * 2007-11-05 (RGL) added SBE41F oxygen sensor to channels 48 and 49.
% * 2009-03-06 (RGL) changed fopen to fopen_odas.
% * 2009-12-22 (RGL) changed to demultiplexing record-by-record for better
%   memory efficiency.
% * 2010-01-15 (AWS) added support for odas v6 and up
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-04-11 (WID) changed inifile_with_instring calls to setupstr
% * 2012-07-19 (WID) speed improvements - reduced load from waitbar
% * 2012-10-23 (WID) documentation update and small logic change when opening
%                    the .mat file.
% * 2012-11-23 (WID) octave compatability change when setting output file name
% * 2013-06-10 (WID) modified to return data in place of tmp .mat filename
% * 2013-06-25 (WID) segmented the processing of modern and legacy files.
%                    Corrected bug with processing legacy files.
% * 2013-07-04 (WID) added support for interpreting unsigned values.
%                    added sorting and error reporting for multiple channel
%                    ids on a single line ie, "id = 16,17".
% * 2014-12-15 (WID) added t_fast and t_slow vectors
% * 2015-01-27 (WID) allow opening files in another directory while saving
%                    the resulting .mat files into the current directory.
% * 2015-02-26 (WID) Corrected bug when reading clipped / corrupted files.
% * 2015-03-12 (RGL) Forced Sea-Bird counter channels to be double
%                    precision by adding 1e-12 to every sample.
% * 2015-04-07 (WID) Added file name to the resulting MAT file / struct.
% * 2015-04-08 (WID) Cleaned up documentation, merged changes from Rolf.
% * 2015-07-30 (WID) Complete rewrite.  Legacy code removed.
% * 2015-07-31 (WID) Fixes for previous versions of Matlab.
% * 2015-08-01 (WID) Big speed increase - now same as before.  Bug fixes.
% * 2015-09-22 (WID) Print out the address matrix to the terminal.
% * 2015-11-01 (RGL) Update documentation
% * 2015-11-18 (RGL) Update documentation.
% * 2016-12-29 (WID) Updated logic for converting signed to/from unsigned.
% * 2017-05-09 (WID) Fixed reading of broken files - last segment now
%                    correctly trimmed.  Some minor logic simplification.
% * 2018-08-28 (JMM) Added a warning message if there is no channel section 
%                    in the setup file for an id number that is in the address 
%                    matrix.
% * 2019-10-03 (WID) Force jac_t types to be unsigned.
% * 2022-05-04 (WID) Fix start of file time to record 0, not record 1.

function [ch_list, outstruct] = read_odas( fname )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Check file name, check for a .mat file, and open the .p file %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ask the user for a file name if not provided as an input variable
if nargin < 1
    fname = input(['Enter data file name (default: ' get_latest_file '): '],'s');
    if isempty(fname), fname = get_latest_file; end
end

[P,filename,E,d.fullPath] = file_with_ext( fname, {'','.p','.P'} );
if isempty(filename), error('Unable to find file %s\n',fname); end

% Always save the '.mat' file into the current directory.  Do not use the
% fullPath variable.
while exist([filename '.mat'], 'file') && nargout ~= 2
	sub_name = input(['.mat file already exists. \nEnter a new name to ' ...
                      'save to a different file, "Enter" to overwrite, ' ...
                      'or "q" to quit: '], 's');

    if strcmpi(sub_name, 'q'), return; end
    if isempty(sub_name), break; else filename = sub_name; end
end

% Open the file
[fid, error_message] = fopen_odas(d.fullPath,'r');
if ~isempty(error_message), disp(error_message); end;
if fid<3
   error('Unable to open file %s\n', d.fullPath);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Constants & parameters                                   %%%%%%
%%%%%% These will not generally need to be changed by the user. %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
bytes_per_word      =  2; % Assume 2-byte integers
header_version_i    = 11; % this contains header version, always 1 prior to odas v6
setupfile_size_i    = 12; % odas v6 and higher only: size of ini string in first record
header_size_i       = 18; % Index of header_size in header record
block_size_i        = 19; % Index of block_size in header record
header_size_0       = 64; % Predicted size of the header in each block (integers)
fast_cols_i         = 29; % Header indices for the number of fast & slow columns
slow_cols_i         = 30;
rows_i              = 31; % Header index for the number of rows
clock_whole_i       = 21; %Index to integer part of clock frequency
clock_fraction_i    = 22; %Index to fractional part of clock frequency

% Parameters read in from the first header in the file
fseek(fid,0,'bof');
HD = fread(fid,header_size_0,'ushort');
header_size = HD(header_size_i)/bytes_per_word;
setupfile_size = HD(setupfile_size_i);
record_size = HD(block_size_i)/bytes_per_word;
data_size = record_size - header_size;
fast_cols = HD(fast_cols_i);
slow_cols = HD(slow_cols_i);
n_cols = fast_cols + slow_cols;
n_rows = HD(rows_i);
%rows_per_record = data_size/n_cols; % rows in a record
%slow_samples_per_record = rows_per_record / n_rows; % samples in a slow channel per record
%n_slow = n_rows*slow_cols; % number of slow channels
f_clock = HD(clock_whole_i) + HD(clock_fraction_i)/1000;
d.fs_fast = f_clock / n_cols; % the sampling rate of the fast channels
d.fs_slow = d.fs_fast / n_rows; % the sampling rate of the slow channels
fseek(fid,0,'eof');
filesize = ftell(fid);

%MSB has major version, LSB has minor version
d.header_version = bitshift(HD(header_version_i), -8) + ...
                     bitand(HD(header_version_i), 255) / 1000;


first_record_size = header_size*bytes_per_word + setupfile_size;
n_records = (filesize - first_record_size)/HD(block_size_i);

if mod(n_records, 1)
    n_records = floor(n_records);
    warning(['The file ' d.fullPath ' does not contain an integer number of records']);
end
if n_records <= 1
    error(['File ' d.fullPath ' contains no data'])
end

% Calculate a new file size after possibly trimming the data file.
filesize = first_record_size + n_records * HD(block_size_i);

matrix_count = floor( n_records * data_size / (n_rows * n_cols));

d.t_slow = (0:matrix_count-1)' / d.fs_slow;
d.t_fast = (0:matrix_count*n_rows-1)' / d.fs_fast;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Read the data file %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fseek(fid,header_size*bytes_per_word,'bof');

d.setupfilestr = fread(fid,setupfile_size,'*char*1')';
if isempty(d.setupfilestr)
    error('failed to extract setup file string from first record');
end

% Parse the configuration file.
cfg = setupstr(d.setupfilestr);
d.cfgobj = cfg;

% Load the matrix
rows = setupstr(cfg, 'matrix', 'row[0-9]*');

ch_matrix = zeros(n_rows, n_cols);
for i = 1:length(rows)
    values = textscan(rows{i}, '%d');
    ch_matrix(i,:) = values{1}';
end
ch_all = unique(ch_matrix);

% Print out the address matrix
fprintf('\nAddress Matrix:\n');
for addr = ch_matrix'
    for x = addr'
        fprintf('%4d', x);
    end
    fprintf('\n');
end
fprintf('\n');

% The special character is not listed as a channel.  Insert it manually if
% it exists within the channel matrix.
if ch_matrix(ch_matrix == 255)
    ch_nums = 255;
    ch_names = {'ch255'};
else
    ch_nums = [];
    ch_names = {};
end

missing_channels = [];
for section_name = setupstr(cfg, '')
    % No id, name, and type -- not a channel.
    if isempty(setupstr(cfg, section_name, 'id')),   continue; end
    if isempty(setupstr(cfg, section_name, 'type')), continue; end
    if isempty(setupstr(cfg, section_name, 'name')), continue; end

    tmpid = eval(sprintf('[%s]', char(setupstr(cfg, section_name, 'id'))));
    if ~isnumeric(tmpid), continue; end

    % Use parameter 'name' from the channel section
    namestr = char(setupstr(cfg, section_name, 'name'));

    if length(tmpid) == 1
        if ~isempty(ch_matrix(ch_matrix == tmpid))
            fprintf(1, '     channel: %2d = %s\n', tmpid, namestr);
            ch_nums(end+1) = tmpid;
            ch_names{end+1} = namestr;
        end
    end

    if length(tmpid) == 2
        if ~isempty(tmpid(1)) && ~isempty(ch_matrix(ch_matrix == tmpid(1)))
            fprintf(1, 'even channel: %2d = %s\n', tmpid(1), namestr);
            ch_nums(end+1) = tmpid(1);
            ch_names{end+1} = [namestr '_E'];
        end

        if ~isempty(tmpid(2)) && ~isempty(ch_matrix(ch_matrix == tmpid(2)))
            fprintf(1, ' odd channel: %2d = %s\n', tmpid(2), namestr);
            ch_nums(end+1) = tmpid(2);
            ch_names{end+1} = [namestr '_O'];
        end
    end
    
    if setdiff( tmpid, ch_all )
        missing_channels = [missing_channels tmpid];
    end
end

% Check if there is a channel section for each channel in address matrix.
% Append the special character to both sets so they cancel each other out.
if ~isempty(missing_channels)
    msg = [ newline ...
        '  The following channel(s) are sampled in the address matrix ' ...
        'but are not defined within a [channel] section. They will ' ...
        'not be processed.' newline];

    for ch_num = missing_channels
        msg = [ msg newline '    Channel ' num2str(ch_num) ];
    end

    msg = [ msg newline newline ...
        '  Modify and patch a configuration file to fix this error(s).' ...
        ' The configuraiton file should have a [channel] section for ' ...
        'each sampled address.' ];
    warning( msg )
end

%     fprintf(2,'*********************************** WARNING ***************************************\n')
%     fprintf(2,'*  Channel %d: Present in address matrix but NOT as an ''id'' in channel section    *\n',ch_num)
%     fprintf(2,'*              Data collected on this channel will not be processed              *\n')
%     fprintf(2,'*              Setup file should be patched with correct ''id'' number              *\n')
%     fprintf(2,'*  Hit any key to continue, CTRL+C to exit                                        *\n')
%     fprintf(2,'***********************************************************************************\n'),pause
% The file pointer should be at the beginning of the second record

% Read remailing file - use "filesize" to ensure we read an integer number
% of records.  Using "inf" will cause problems with corrupted data files.
fileData = fread(fid, (filesize - first_record_size)/bytes_per_word, 'short');
fclose(fid);

fileData = reshape(fileData, header_size + data_size, []);

% Convert the header elements into unsigned values - they were read as
% shorts.
header = fileData(1:64,:);
header(header < 0) = header(header < 0) + 2^16;
d.header = reshape(header, 64, []);

fileData = fileData(65:end,:);

% Extract the data for each specific file.
for i = 1:length(ch_names)
    indices = ch_matrix' == ch_nums(i);
    data = reshape(fileData, numel(ch_matrix), []);
    result = data(indices,:);
    ch.(ch_names{i}) = result(:);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Do some basic conversions: extract the time and date from the %%%%%%
%%%%%% header, and combine split channels                            %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract date/time from the header
d.filetime = datenum([HD(4:8)' HD(9)+HD(10)/1000]);
d.date = datestr(d.filetime,'yyyy-mm-dd');
d.time = datestr(d.filetime,'HH:MM:SS.FFF');


% Combine split channels (if requested in the setup portion of the routine)
% and slow channels that are sampled more than once per scan

for name = fieldnames(ch)'
    if length(name{1}) > 2 && ~strcmp(name{1}(end-1:end), '_E')
        % Join channels when we see an "_E" even channel.
        continue;
    end
    even_ch = name{1};
    odd_ch = [name{1}(1:end-2) '_O'];
    if ~isfield(ch, odd_ch), continue; end

    chE = ch.(even_ch);
    chE(chE < 0) = chE(chE < 0) + 2^16;
    chO = ch.(odd_ch);
    chO(chO < 0) = chO(chO < 0) + 2^16;

    joined = chO * 2^16 + chE;
    ch = rmfield( ch, [{even_ch}; {odd_ch}] );
    ch.(name{1}(1:end-2)) = joined;
end


% Correct for channels that should be signed.
for name = fieldnames(ch)'
    sign = char(setupstr(cfg, name, 'sign'));
    type = char(setupstr(cfg, name, 'type'));
    
    % JAC ADCs are always unsigned but people forget to mark channel
    % correctly.
    if strcmpi(type,'jac_t')
        sign = 'unsigned';
        
    % These types have already been converted and should be skipped.
    elseif strcmpi(type,'sbt') || ...
            strcmpi(type,'sbc') || ...
            strcmpi(type,'jac_c') || ...
            strcmpi(type,'o2_43f')
        continue
    end
    
    if strcmpi(sign, 'unsigned') 
        % Find and adjust unsigned "negative" values.
        ch.(name{1})(ch.(name{1})<0) = ch.(name{1})(ch.(name{1})<0) + 2^16;
    else
        % Find and adjust values too large to be signed.  This only applies
        % to 32 bit values - all 32 bit channels not identified above are,
        % by default, signed.
        ch.(name{1})(ch.(name{1})>=2^31) = ...
            ch.(name{1})(ch.(name{1})>=2^31) - 2^32;
    end
end


% Take the transpose of the header - for backwards compatibility.
d.header = d.header';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Save the variables to file. Note that this can be time- %%%%%%
%%%%%% consuming for a large file.                             %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ch_list = {};
for name = fieldnames(ch)'
    d.(name{1}) = ch.(name{1});
    ch_list{end+1} = name{1};
end


% -------------------------------------------------------------------------
% Add the version of the ODAS library used to create the resulting .MAT
% file.
% -------------------------------------------------------------------------
d.odas_version = odas_version_info();


if nargout < 2
    disp(['Saving raw data into file: ' filename]);
    save(filename, '-struct', 'd', '-v6' ); % save all variables.
else
    outstruct = d;
end


% Added to test the return of a struct containing the contents of the .MAT
% file.  This minimizes the number of file IO opperations by a factor of 5.
%
% To reconstruct the variables as though you loaded the mat file and
% assuming you called read_odas as follows, perform the next command.
%
% >> [fields, data] = read_odas( 'datafile.p' );
% >> for c=fieldnames(data)', eval([char(c) ' = data.' char(c) ';']); end

end
