%% read_tomi
% Read a range of records from an ODAS data file.
%%
% <latex>\index{Functions!read\_tomi}</latex>
%
%%% Syntax
%
%   [blocks, header, data] = read_tomi( file, start_block, end_block, type )
%
% * [file] file descriptor for an open ODAS data file
% * [start_block] block index of the first segment block
% * [end_block] block index of the last segment block - inclusive.  Use
%               "NaN" to specify the end of file.
% * [type] data type contained within file.  Optional - defaults to 'int16'
%          for 16 bit data words.  Use 'uint8' to extract bytes.
% * []
% * [blocks] vector of the extracted record numbers
% * [header] record header
% * [data] requested record data
%
%%% Description
%
% Read a range of blocks from an ODAS data file. This function is primarily
% used by patch_odas but can also be used directly for the purpose of
% processing large data files into manageable portions. 
% 
% Also used to extract data from ODAS serial port data files.  When used
% with serial port data, specify type as 'uint8'.
%
%%% Examples
%
%    >> [fid, error] = fopen_odas( 'raw_data_file.p', 'r' );
%    >> if isempty( error ),
%    >>     [blks, header, data] = read_tomi( fid, 1, 1 );
%    >>     fclose( fid );
%    >> end
%
% Read the first block from the data file 'raw_data_file.p'.
%
%    >> [blks, header, data] = read_tomi( 'raw_data_file.p', 1, 1 );
%
% Read the first block from the data file 'raw_data_file.p' without having
% to open and close the data file.
%
%    >> [fid, error] = fopen_odas( 'serial_data_file.s', 'r' );
%    >> if isempty( error ),
%    >>     [blks, header, data] = read_tomi( fid, 1, NaN, 'uint8' );
%    >>     fclose( fid );
%    >>     data = char(data)';      % Convert into a character string.
%    >> end
%
% Extract all the data encoded within the 'serial_data_file.p' file.  The
% resulting 'data' variable will contain a copy of the data observed on the
% serial port without any ODAS added headers.

% Version History
%
% * 2010-01-15 (AWS) support for odas v6 and up
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-10-23 (WID) documentation added
% * 2013-03-21 (WID) added support for extracting data from serial files
% * 2017-03-23 (WID) allow file ID or file name input parameter

function [blocks, header, data] = read_tomi(file, start_block, end_block, type)

% set the default value for "type" - a 16bit signed word.
if nargin <= 3,
    type = 'int16';
end

if ~ischar(file)
    fid = file;
else
    if isempty(file)
        file = input(['Enter data file name (default: ' get_latest_file '): '],'s');
        if isempty(file), file = get_latest_file; end
    end
    
    [P,name,E,path] = file_with_ext( file, {'','.p','.P'} );
    if isempty(name), error('Unable to find file %s\n', file); end
    
    % Open the file
    [fid, error_message] = fopen_odas(path,'r');
    if ~isempty(error_message), warning(error_message);  end
    if fid<3, error('Unable to open file %s\n', path); end   
end

fseek(fid,0,'bof');

head = fread (fid,64,'ushort');
header_size = head(18)/2; % in units of 2 byte words
record_size = head(19)/2;
data_size = record_size - header_size;
header_version_i = 11;

header_version = bitshift(head(header_version_i), -8) + ...
                 bitand(head(header_version_i), 255) /1000;

if (header_version >=6)
    setupfile_size = head(12);
    if setupfile_size == 0,
       first_record_size = 0;
    else
       first_record_size = setupfile_size + header_size * 2;
    end
else
    first_record_size = record_size * 2;
end

% Allow "NaN" to be used to specify the entire file.  If specified,
% calculate the last block and use as end_block.
if isnan(end_block),
    fseek( fid, 0, 'eof' );
    file_length = ftell( fid );
    end_block = (file_length - first_record_size) / (2 * record_size);
end
    

data      = zeros(2*data_size*(end_block - start_block +1),1,'uint8');
d_buffer  = zeros(2*data_size,1,'uint8');
header    = ones (header_size*(end_block - start_block +1),1)*NaN;
h_buffer  = ones (header_size,1)*NaN;


fseek(fid, first_record_size + 2*(header_size+data_size)*(start_block-1),'bof');

% there is no data in block zero. 

% Read the file
for index = 1:(end_block - start_block + 1)

% Now read the header
   [h_buffer, count] = fread (fid, header_size, 'short');
   if count ~= header_size,
    disp(['Failed to read header at block = ' num2str( start_block + index -1)])
        blocks = index -1;
        header = header(1:(index-1)*header_size);%trim stuff we failed to read
        data   = data  (1:(index-1)*2*data_size);
        break
   end

% Now read the data
   [d_buffer, count] = fread (fid, 2*data_size, '*uint8');
   if count ~= 2*data_size
    disp(['Failed to read data at block = ' num2str(start_block +index -1)])
        blocks = index -1;
        header = header(1:(index-1)*header_size);%trim stuff we failed to read
        data   = data  (1:(index-1)*2*data_size);
        break
   end

   data (1+((index-1)*2*data_size):index*2*data_size) = d_buffer;      %fill in
   header (1+((index-1)*header_size):index*header_size) = h_buffer;%fill in

   blocks = index;
end

% Now must convert data into the requested type - account for endian format
[a,b,mendian] = computer;
dendian = 'B';
if head(64) == 1, dendian = 'L'; end

% Convert the data into the requested type.  This step is needed for proper
% use of the "swapbytes" function.
data = typecast(data, type);

% Swap the bytes if the machine that recorded the data and the current PC
% have different endian formats.
if mendian ~= dendian, data = swapbytes(data); end

% If not specifically requested, return the data as double values.
if nargin <= 3,
    data = cast(data, 'double');
end

if ischar(file)
    fclose(fid);
end

