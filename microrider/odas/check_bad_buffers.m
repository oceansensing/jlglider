%% check_bad_buffers
% Find bad buffers within a data file
%%
% <latex>\index{Functions!check\_bad\_buffers}</latex>
%
%%% Syntax
%   [badRecord] = check_bad_buffers( file )
%
% * [file] Name of data file with or without the '.p' extension.
% * []
% * [badRecords] Vector index to records with bad buffers.
%
%%% Description
% Scan a RSI raw binary data file for bad buffers. Communication problems can 
% sometimes cause the data acquisition software to drop data. These events are 
% detected and recorded within the data acquisition log file. The header in the
% corresponding data file is also flagged.
%
% This fuction searches the headers within a RSI binary data file for the bad 
% buffer flag. The indexes of these records are returned in BadRecord.
%
%%% Examples
%
%    >> badIndexes = check_bad_buffers( 'my_data_file.p' )
%
% Search for bad buffers within the file $\texttt{my\_data\_file.p}$. If
% $\texttt{badIndexes}$ is empty, no bad buffers were found.

% Version History
%
% * 2012-02-02 (RGL) initial
% * 2012-04-30 (WID) allow for extensions of either case (p|P)
% * 2012-04-30 (WID) comments added to allow for Matlab publishing
% * 2012-07-19 (WID) fix of cast and endian problems
% * 2012-08-31 (WID) use of file_with_ext function
% * 2012-12-11 (WID) Documentation update.
% * 2012-10-27 (RGL) Document corrections.
% * 2012-11-18 (RGL) Document corrections.

function  [bad_records] = check_bad_buffers(file)

if nargin < 1, file = get_latest_file; end

error_str = ['File "' file '" not found.  Only data files (.p) are valid.'];
[P,N,E,file] = file_with_ext( file, {'.p', '.P', ''}, error_str );

[fid, error_message] = fopen_odas(file,'r+');
if error_message, warning(error_message); end
if (fid < 3), error('Could not open input file = %s', file); end

% Set up some useful constants
bpw = 2; %Assume that the data is in 2-byte integers
junk = fread(fid, 64, 'ushort'); % read first header
header_size = junk(18)/bpw; % in units of 2 byte words
record_size = junk(19)/bpw;
data_size   = record_size-header_size;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ver = cast(junk(11), 'uint16');
header_version = cast(bitshift(ver,-8), 'double') + ...
                 cast(bitand(ver,255), 'double')/1000;

fseek(fid,0,'bof');

if header_version >= 6,
    first_record_size = junk(12) + header_size*2;
else
    first_record_size = record_size * bpw;
end

fseek(fid,0,'eof');
filesize_in_bytes = ftell(fid);
total_records = (filesize_in_bytes - first_record_size)/(record_size*bpw);

%position file pointer to second record that contains real data, first record contains
%setup file in string format
fseek(fid, first_record_size,'bof');

if total_records ~= floor(total_records),
    warning(['The file ' file ' does not contain an integer number of records.']);
end

header = zeros(floor(total_records)+1, header_size);

% Read the headers from the file
for index = 1:total_records
   % read the header
   header(index,:) = fread(fid, header_size,'ushort')';      
   %At this point a header has been read successfully, skip to next header
   fseek(fid, data_size*bpw, 'cof');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We now have all headers for this file. So look for failed special character check flag.
bad_records = find(header(:,16) ~= 0);% index to all failed special character checks

if isempty(bad_records),
   display(['No bad records were found in file ' file]);
end

fclose(fid);

