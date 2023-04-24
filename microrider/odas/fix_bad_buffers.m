%% fix_bad_buffers
% Fix badly corrupted records in a RSI raw binary data file
%%
% <latex>\index{Functions!fix\_bad\_buffers}</latex>
%
%%% Syntax
%   [badRecords, chop, truncate] = fix_bad_buffers( fileName )
%
% * [fileName] Full name of a RSI raw binary data file including the '.p' extension.
% * []
% * [badRecords] Vector of the record numbers of bad records that were detected 
%         and fixed in the file 'fileName' after the fixing.  The record numbers
%         may shift if the program chops records off the front of the file.
% * [chop] Number of bad records chopped off the front of the file. Set to -1 if 
%         the file can not be opened.
% * [truncate] Number of bad records truncated from the file.  Set to -1 if the 
%         file can not be opened.
%
%%% Description
% Repair bad records in RSI raw binay data files. It sometimes happens that there is a 
% communication failure with an instrument. This is indicated by the 'Bad Buffer'
% message during data acquisition and its occurrence is flagged in the header of
% the affected record. The result is some skewed data in the binary file that can
% lead to huge errors if the skew is not fixed before processing the data.
%
% When errors are detected, this function replaces the entire record with 
% interpolated data.  All significant data in the record will be lost but the
% resulting record can still be graphed without disrupting neighboring
% records. However, because there will be a time gap between a record with
% a bad buffer and the next record, the deconvolution function may produce
% erroneous results. The error will be transitory and decay with an
% e-folding time equal to the gain of the differentiator used to 
% pre-emphasize a signal.
%
% If a file ends with bad records, those records are truncated. If a
% file starts with bad records, those records are also truncated. If a
% bad record has already been fixed successfuly using
% $\texttt{patch\_odas}$, the record is not changed. 
%
% IMPORTANT: Bad buffers should first be fixed using $\texttt{patch\_odas}$
% because it preserves more of the original data. Should this fail,
% $\texttt{fix\_bad\_buffers}$ can be used to patch the data file.
%
% @image @images/bad_buffers @Corrupted ODAS binary data file @A corrupted
% data file that has 'Bad Buffers' from about 45 to 48 seconds. This sample
% is from an Expendable Microstructure Profiler (XMP) which was notorious
% for communication failuresn in the early stages of its development.
% Records 45 through 48 are skewed and appear garbled as a result of
% communication failures. 
%
% @image @images/bad_buffers_fixed @Corrupted ODAS binary data file after
% being repaired. @Same data as in the previous figure after being
% processed by $\texttt{patch\_odas.m}$ and then by
% $\texttt{fix\_bad\_buffers.m}$. Three records have been replaced with
% interpolated data. Interpolated data should be used cautiously.
%
% This function will make a backup copy, but the user is advised to always
% make a backup of the data file in a secure directory before using this
% function. This function should be called *after* $\texttt{patch\_odas}$
% because less data is lost if $\texttt{patch\_odas}$ is run first.

% Version History
%
% * 2008-12-08 RGL
% * 2009-03-06 RGL Changed fopen to fopen_odas
% * 2010-01-14 AWS support odas v6 and higher
% * 2011-09-01 AWS added documentation tags for matlab publishing
% * 2012-04-11 WID replaced inifile_with_instring calls with setupstr
% * 2012-06-18 WID using 'file_with_ext' for determining file name
% * 2012-09-03 WID added check to prevent overwriting of backup file
% * 2012-11-05 WID documentation update
% * 2013-02-26 WID merged changes from Rolf
% * 2013-07-30 WID fix of s/file_name/filename/g bug
% * 2015-10-28 (RGL) Documentation changes.
% * 2015-11-18 (RGL) Documentation changes.

function [bad_records, chop, truncate] = fix_bad_buffers(file)

bad_records = []; chop = []; truncate = []; % pre-assign to nulls
[P N E filename]= file_with_ext( file, ...
                                 {'.p', '.P', ''}, ...
                                 ['Could not find file: "' file '"']);

if exist([filename '.bad_buffers'], 'file') ~= 0,
    warning_string = 'Backup file already present.  Do you wish to overwrite?  (y/N) :';
    reply = input( warning_string, 's' );
    if ~strncmp(reply, 'y', 1),
        return;
    end
end

success = copyfile(filename, [filename '.bad_buffers']); % make a back-up copy
if ~success, disp(['unable to back up ' filename]); return, end

% Open the 2-byte binary data file for reading
[fid, error_message] = fopen_odas(filename,'r+'); % open for rw but no truncate
if ~isempty(error_message), disp(error_message); end
if (fid < 3);
    chop = -1; truncate = -1; % these can be used to check for error conditions
%    error(['Could not open input file = ' filename']);
    return
end
fseek(fid,0,'bof'); % put pointer at the beginning of file

% Set up some useful constants
header_version_i = 11;   %offset to header version word

bytes_per_word = 2; %Assume that the data is in 2-byte integers
junk = fread (fid,64,'ushort'); % read first header
header_size  = junk(18)/bytes_per_word; % in units of 2 byte words
record_size  = junk(19)/bytes_per_word;
data_size    = record_size-header_size;
fast_columns = junk(29);% Matrix information is located in the header
slow_columns = junk(30);
columns      = fast_columns + slow_columns;
rows         = junk(31);
sp_char      = 32752; % numeric value of special character
scan_length  = columns*rows; %the number of elementsd in the basic scan matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

header_version = bitshift(junk(header_version_i), -8) + ...
                 bitand(junk(header_version_i), 255) /1000;

if (header_version >=6)

    if (junk(64) == 1)
        file_endian = 'l';
    else
        file_endian = 'b';
    end

    setupfile_size = junk(12);

    fseek(fid,0,'bof'); % put pointer at the beginning of file

    first_record = fread(fid,setupfile_size + header_size*2);
    first_record_length = length(first_record);

    setupfile_str = char(first_record(64:first_record_length))';
    
    cfg = setupstr(setupfile_str);

    rows = setupstr(cfg, 'matrix', 'row[0-9]+');
    address_matrix = [];
    
    for row = rows
      values = textscan(row{1}, '%f');
      address_matrix = vertcat(address_matrix, values{1}');
    end
% next line added 2013-02-26
    rows = size(address_matrix, 1);

else
    address_matrix = fread(fid, scan_length, 'short'); %get address matrix

    first_record_length = bytes_per_word * record_size;
end

index_to_sp_char = find(address_matrix == 255); % this is the location of
% the special character in the data file. It is usually at the beginning
% and equal to 1.

fseek(fid,0,'eof'); % put pointer at end of file
filesize_in_bytes = ftell(fid); % File size in bytes

total_records = (filesize_in_bytes - first_record_length)/(record_size*bytes_per_word);
%position file pointer to second record that contains real data, first record contains
%setup file in string format (v6 and up) or matrix (pre v6)
fseek(fid, first_record_length,'bof');

if (total_records - floor(total_records) ~= 0)
  warning('The file %s.p does not contain an integer number of records', filename);
end

header = zeros(total_records, header_size);

fseek(fid, first_record_length, 'bof'); %position file pointer to first real record
% as there is no data in block zero,  only the setup file in string form.

% Read the headers from the file
for index = 1:total_records
   junk = fread (fid, header_size,'ushort');% read the header
   header(index,:) = junk'; % move into matrix
   fseek(fid,data_size*bytes_per_word,'cof');% skip over data
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We now have all headers for this file. So look for failed special character check flag.

bad_records = find(header(:,16) ~= 0);% index to all failed special character checks
if isempty(bad_records)
   fprintf(1, 'No bad records were found in file %s', filename);
   return
end


% Now check that the records really are bad and have not been fixed
% before by this or some other program, like patch_tomi (or patch_odas)
remove_list = []; % The list of records flagged to be bad but are actually OK 
for index = 1:length(bad_records)
    [blocks, header, data]=read_tomi(fid,bad_records(index), bad_records(index));% get the bad record.
    data = reshape(data, scan_length, data_size/scan_length);
    n = find(data(index_to_sp_char,:) ~= sp_char); % empty if all columns have a special character
    if isempty(n)
        remove_list = [remove_list index]; % add to list of OK records.
    end
end
bad_records(remove_list) = []; % Remove the records that are OK because they may have been fixed already
if isempty(bad_records)
    display('No bad records found betwen second and second-last records');
    return;
end


% First trancate if bad blocks are at the end
index = find (bad_records == total_records); % Remember that the first record does not have data
if (~isempty(index))
%    'Bad record at end of file'
    truncate = index;
    % trancate the file at this record.
    got_more = 1; % check if the previous record is also bad
    while(got_more)
        if any(bad_records == (bad_records(truncate) - 1));
            truncate = truncate - 1;
        else
            got_more = 0;
        end
    end
    fseek(fid, 0 , 'bof'); %position file at begin

    junk = fread(fid, first_record_length + record_size*(bad_records(truncate)-1)*2); % read all except last bad record(s)

    fclose(fid);
    
    fid2 = fopen(filename ,'w+', file_endian); % Open and truncate
    if ~isempty(error_message)
        error(error_message);
    end
    
    if fid2 < 3
        error('unable to open file %s for writing', filename);
    end
    fwrite(fid2, junk); % write everything except last records
    fclose(fid2);

%    ['truncated file from record ' num2str(bad_records(truncate)) ' to end']
    truncate_record = bad_records(truncate); % save the record number
    bad_records(truncate:end) = []; % remove from list of bad_records
    truncate = truncate_record; % save for return value of this function
    
    if isempty(bad_records)
%        ['Got all bad records by truncating the tail of the file ' filename]
        return
    end
    [fid, error_message] = fopen_odas(filename ,'r+'); % open for rw and no truncate
    if ~isempty(error_message), error_message, end;
end



% Next chop off some of the beginning if bad_records occur in the first records
index = find (bad_records == 1);
if (~isempty(index))
%    'Bad records near beginning of file'
    chop = 1; 
    got_more = 1; % check if the next record is also bad
    while(got_more)
        if any(bad_records == (chop + 1))
            chop = chop + 1;
        else
            got_more = 0;
        end
    end
    % chop off all records up to and including this one.
    fseek(fid, 0 , 'bof'); %position file pointer to start of file
    junk = fread(fid, inf); % read all
    fclose(fid);

    [fid2, error_message] = fopen(filename,'w', file_endian); % Open and truncate
    if ~isempty(error_message), error_message, end;
    
    if fid2 < 3
        error('unable to open file %s for writing', filename);
    end

    fwrite(fid2, junk(1:first_record_length)); % write first record which contains ini string
    fwrite(fid2, junk((chop + 1)*record_size + 1 : end), 'short'); % write all except first few records
    fclose(fid2);

    [fid, error_message] = fopen_odas(filename,'r+');
    if ~isempty(error_message), error_message, end;
%    ['chopped off the first ' num2str(bad_records(chop)) ' from file']

    % Now re-establish the list of bad records because the indeces have been changed by chopping   
    fseek(fid,0,'eof'); % put pointer at end of file
    filesize_in_bytes = ftell(fid); % File size in bytes

    total_records = (filesize_in_bytes - length(first_record))/(record_size*bytes_per_word);
    total_records = floor(total_records);

% Predefine data buffer. Remember that the first record does not contain data
   header = zeros (total_records, header_size);
   fseek(fid,length(first_record), 'bof'); %position file pointer to first record

   % Read the headers from the file
    for index = 1:total_records
        junk = fread (fid, header_size,'ushort');% read the header
        header(index,:) = junk'; % move into matrix
        fseek(fid,data_size*bytes_per_word,'cof');% skip past data
    end
    bad_records = find(header(:,16) ~= 0);% index to all failed special character checks
    if (isempty(bad_records))
%        ['NO bad records left in file ' filename ' after chopping front-end of file.']
        fclose('all');
        return
    end
end

% The above proceedure quarantees that there is a good record at the start
% and end of the file which can, if necessary, be used for interpolation

% Now deal with the rest of the bad records
% Be sure to also handle consecutive bad records
bad_records_for_return = bad_records; % save teh list for the return of this function
while ~isempty(bad_records)
    before = bad_records(1) - 1; % the good record before the first bad one
    after = before + 2; % a possibly good record after the first bad one
    consecutives = 0;
    if (length(bad_records) > 1), consecutives = 1; end % check for consecutive bad_records
    while(consecutives)
        if any(bad_records == after)
            after = after + 1;
        else
            consecutives = 0;
        end
    end
    n = find(bad_records <= after-1);
    bad_records(n) = []; % remove them from the list
    
    [blocks, header,data]=read_tomi(fid,before, before);% get the record before it
    scan_before = data(end - scan_length + 1 : end); % this is the last sweep through the address matrix before the bad record
    [blocks, header,data]=read_tomi(fid,after, after);% get the record after the bad_record(s)
    scan_after = data(1:scan_length); % this is the first sweep through the address matrix after the bad record
    
    % reshape to find average values of fast channels
    scan_before = reshape(scan_before, columns, rows);
    scan_before(slow_columns+1:end,:) = repmat(mean(scan_before(slow_columns+1:end,:),2),1,rows);
    scan_before = scan_before(:); % Turn back into a column vector
    scan_after = reshape(scan_after, columns, rows);
    scan_after(slow_columns+1:end,:) = repmat(mean(scan_after(slow_columns+1:end,:),2),1,rows);
    scan_after = scan_after(:); % Turn back into a column vector
        
    new_data = zeros(scan_length, (after-before-1)*data_size/scan_length); % a data matrix for filling in
    delta_data = (scan_after - scan_before);% The change across the bad record(s)
    data_gradient = delta_data/(size(new_data,2)); % the change from column to column
    new_data(:,1) = scan_before + data_gradient;
    for pointer=2:size(new_data,2)
        new_data(:,pointer) = new_data(:,pointer-1) + data_gradient;
    end
    
    fseek(fid, first_record_length + bytes_per_word*record_size*(before), 'bof'); %position file pointer to start of bad record
    new_data = new_data(:); % turn it into a vector for easier indexing
    for index = 1 : (after-before-1) % this is the number of records that get filled
        fseek(fid,bytes_per_word*header_size,'cof'); % move past header
        buffer = new_data(1:data_size);
        new_data(1:data_size) = [];
        fwrite(fid, buffer , 'short'); % fill in new interpolated data
    end
end
bad_records = bad_records_for_return; % tell the user which ones were fixed

        
fclose all;

