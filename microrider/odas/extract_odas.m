%% extract_odas
% Extract a range of records from a RSI raw binary data file.
%%
% <latex>\index{Functions!extract\_odas}</latex>
%
%%% Syntax
%   [newFile, successFlag] = extract_odas( fileName, startRecord, endRecord )
%
% * [fileName]    Name of the binary data file, extension optional.
% * [startRecord] Record number from which to start copying records.
% * [endRecord]   Record number of the last required record.
% * []
% * [newFile] Name of the file created by this function. The first record
%        contains the configuration string. This allows the extracted file
%        to be processed like any other raw data file. All records from 
%        startRecord to endRecord are included in the new file.
% * [successFlag] Integer representing status.  0 if operation successful,
%        otherwise <= -1.
%
%%% Description
% Extract a range of records from a RSI raw binary data file and move them
% into a new file. The new file has the same name as the old file except
% that the range of the records is appended to the name. The original file
% is not altered.
%
% Binary data files are sometimes overwhelmingly long. This function
% provides a means to segment a file into shorter pieces to speed up the
% data processing.  The companion function 'show_P' can be used to extract
% the record-average pressure from a binary file which might be useful for
% determining the appropriate range(s) to be extracted.
%
% This function can be used iteratively to segment a raw binary data file
% into an arbitrary number of smaller files.  See also,
% $\texttt{show\_P.m}$. 
%
%%% Examples
%
% The following example shows how the data file ``DAT001.p" can be 
% segmented to produce a smaller file. For example, if the region of 
% interest is between records 5400 and 6900.
%
%    >> extract_odas( 'DAT001', 5400, 6900 )
%    ans = DAT001_5400_6900.p
%
% The resulting file can be processed normally and, because it is smaller, will
% process faster.
%
% Additional examples are in the top of the function itself.
 
% The following example shows the first and first + second records being
% extracted from the data file titled "DAT001.p".
%
%    >> extract_odas( 'DAT001', 1, 1 )
%    ans = DAT001_1_1.p
%    >> extract_odas( 'DAT001', 1, 2 )
%    ans = DAT001_1_2.p
%    >> files = dir( 'DAT001*' );
%    >> for file = files', disp([file.name '  ---  ' ...
%                                num2str(file.bytes)]); end
%    DAT001.p  ---  7082595
%    DAT001_1_1.p  ---  10595
%    DAT001_1_2.p  ---  18915
%
% Two data files are created from the original data file "DAT001.p".  The file
% "DAT001_1_1.p" contains only the first record.  The second file "DAT001_1_2.p"
% contains the first two records.  The file sizes for the three files are then 
% printed to the console.
%
% To better understand how a data file is constructed, one can analyse the 
% resulting file sizes of the extract_odas function. The size of a data record 
% is observed to be the difference between a file with one record and two 
% records.  From the above example;
%
%                      Data Record Size = 18915 - 10595 
%                                       = 8320 [bytes]
%
% Each data file has a prepended header and configuration string. The size of
% the prepended header and configuration string can be calculated as;
%
%    Header + Configuration String Size = 10595 - 8320 
%                                       = 2275 [bytes]
%
% And finally, we observe that the original data file consists of a whole number
% of data records plus the prepended header and configuration string.  Notice 
% the result is a whole number.
%
%              Record Count in DAT001.p = (7082595 - 2275) / 8320 
%                                       = 851 [records]

% Version History
%
% * 2009-03-31 RGL
% * 2009-11-03 AWS adapted to handle odas v6 and higher
% * 2011-09-01 AWS added documentation tags for matlab publishing
% * 2012-05-07 WID update to documentation
% * 2012-08-30 WID make use of file_with_ext function
% * 2012-11-05 WID documentation update
% * 2013-05-06 WID added bounds check for when start_record == 0
% * 2015-07-27 WID Documentation update.
% * 2015-10-28 (RGL) Documentation corrections.
% * 2015-11-18 (RGL) Documentation corrections.
% * 2017-01-27 (WID) Added hidden option (4th argument) for optional
%                    destination file.

function [new_file, success_flag] = extract_odas(file_name, start_record, end_record, destination)

success_flag = -1; % preset to failure

% Look for the .p extension and deal with problems
[P,N,E,file_name] = file_with_ext( file_name,           ...
                                   {'' '.p' '.P'},      ...
                                   ['Unable to find file: ' file_name] );

[fid_in, error_message] = fopen_odas(file_name,'r');
if ~isempty(error_message)
    error(error_message);
end
if (fid_in < 3); error(['Cannot open ' file_name]); end
fseek(fid_in, 0, 'bof'); % Move to start of file



header = fread (fid_in, 64, 'ushort');
header_version_i = 11;
endian_flag = header(64);
header_size = 64;
record_size = header(19)/2;
setupfile_size = header(12);


if ~start_record,
    warning('The first data record is at index 1.  Using 1 in place of %d\n', start_record);
    start_record = 1;
end


%MSB has major version, LSB has minor version
header_version = bitshift(header(header_version_i), -8) + bitand(header(header_version_i), 255) /1000;

if(header_version >= 6)
    setupfile_size = header(12);
    length_of_first_record_in_bytes = setupfile_size + header_size*2;
else
    length_of_first_record_in_bytes = 2*record_size;
end

fseek(fid_in, 0, 'bof'); % Move to start of file

%this record is special
[first_record, count] = fread(fid_in, length_of_first_record_in_bytes, 'uchar');

if (count ~= length_of_first_record_in_bytes)
    error('Unable to read first record in data file = %s\n',file_name);
end

% How long is this file??
fseek(fid_in, 0, 'eof'); % move to end of file
length_of_file = ftell(fid_in); % length in bytes

n_records = (length_of_file - length_of_first_record_in_bytes)/(record_size*2);
n_records = floor(n_records); % use floor in case there is a runt record.

% Adjust end_record in case the user mistakenly wants to go beyond the
% end-of-file.
if (end_record > n_records), end_record = n_records; end

fseek(fid_in, length_of_first_record_in_bytes + 2*record_size*(start_record-1), 'bof'); %position file pointer to start_record
junk = fread(fid_in, record_size*(1+end_record - start_record), 'short'); % read the records
fclose (fid_in);

if nargin >= 4
    new_file = destination;
else
    new_file = [N '_' num2str(start_record) '_' num2str(end_record) '.p'];
end
if(endian_flag == 0 || endian_flag == 1), endian_string = 'l'; end;
if(endian_flag == 2),                     endian_string = 'b'; end;
fid_out = fopen (new_file, 'w', endian_string); % give OP file the same endian as IP file
if (fid_out < 3); error(['Cannot create the new file, ' new_file]); end
fseek(fid_out, 0, 'bof');
fwrite(fid_out, first_record, 'uchar');
fwrite(fid_out, junk, 'short');
success_flag = 0;
fclose(fid_out);

