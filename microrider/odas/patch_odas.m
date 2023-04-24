%% patch_odas
% Find and fix bad buffers in a RSI raw binary data file that can be
% corrected by a simple shift and interpolation 
%%
% <latex>\index{Functions!patch\_odas}</latex>
%
%%% Syntax
%   [bad_record, fix_manually] = patch_odas( file )
%
% * [file] Name of RSI raw binary data file with bad records.
% * []
% * [bad_record] Vector of indices to bad records.
% * [fix_manually] Vector of indices to records that can not be fixed.
%
%%% Description
% It sometimes happens that a communication failure with an 
% instrument results in some data loss. This is indicated by 
% "Bad Buffer" messages during data acquisition and its
% occurrence is flagged in the header of the effected record. The result is
% some skewed data in the binary file that can lead to huge errors if the
% skew is not fixed before applying the standard data processing tools.
%
% This function is the first step when attempting to correct bad records
% in a RSI binary data file. It corrects small errors in a 
% record without damaging the surrounding data, by locating 
% the missing value and inserting a "best guess" approximation for that
% datum.
%
% This only works when damage to a data file is minor. If the damage in a
% record is extensive, the algorithm leaves such records unchanged. The
% indices to such records are returned in the $\texttt{fix\_manually}$
% vector. 
%
% Records that can not be fixed with this function can be fixed with the
% $\texttt{fix\_bad\_buffers}$ function. The $\texttt{fix\_bad\_buffers}$
% function replaces all data within a record using linear interpretation.
%
% Warning, this function changes the data file.  Make a back-up copy of
% your data file before using this function.
% 
% @image @images/bad_buffers2 @Example data file with four bad records @An 
% example of a raw binary data file that has 4 bad records [70 153 159 160],
% plotted using $\texttt{plot\_VMP}$. This data is from a prototype
% Expendable Microstructure Profiler. The legend is a wish list of
% channels. Only a subset exist in this instrument.
%
% @image @images/bad_buffers_fixed2 @Repaired data file @Plot of the 
% previous figure after using patch_odas. The data in records [153 159 160]
% had their skew corrected by inserting the missing datum. Corrected
% records can now be used for data processing. However, record 70 was not
% fixed because it had too many missing data points. The function
% fix_bad_buffers can patch this record with interpolated data but must be
% used cautiously. Such records should not be used for dissipation
% calculations. 
%
%%% Examples
%
%    >> copyfile( dataFile, 'original_data_file.p' );
%    >> [bad_records, fix_manually] = patch_odas( dataFile );
%    >> if ~isempty( fix_manually ),
%    >>     fix_bad_buffers( dataFile );
%    >> end
%
% Repair a data file with bad buffers. First repair the file with
% $\texttt{patch\_odas}$, and then with $\texttt{fix\_bad\_buffers}$, if
% necessary. 

% Version History
%
% * 2000-05-22 (RGL)
% * 2001-12-12 (TR) No longer assume only slips by one. Corrects multiple
%  slips in each block. Shouldn't need to be run more than once.
%  fix_manually contains blocks with abnormally large slips which need to 
%  be looked at manually. They are the only blocks not fixed by one pass of
%  this version of patch_odas.
% * 2003-09-15 (IG) Added a case to deal with situations where a very large
%  shift occurs at the beginning  of a file (i.e. i(1)\>scan_length below).
%  Such records must be fixed manually.
% * 2009-03-06 (RGL) changed fopen to fopen_odas.
% * 2010-01-15 (AWS) adapted to handle new odas v6 data format
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-06-18 (WID) using 'file_with_ext' for determining file name
% * 2012-09-07 (WID) update of documentation for matlab publishing
% * 2015-10-30 (RGL) Changed documentation

function  [bad_record, fix_manually] = patch_odas(file)

[P N E file] = file_with_ext( file, ...
                              {'.p', '.P', ''}, ...
                              ['Could not find file: "' file '"']);

[fid, error_message] = fopen_odas(file,'r+');% add .p extension for parallel interface and open file
if ~isempty(error_message), warning(error_message), end;
if (fid < 3);error(['Could not open input file = ' file]);end
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

header_version = bitshift(junk(header_version_i), -8) + bitand(junk(header_version_i), 255) /1000;

fseek(fid,0,'bof');

if (header_version >=6)
    setupfile_size = junk(12);
    first_record_size = setupfile_size + header_size*2;
else
    first_record_size = record_size*bytes_per_word;
end

fseek(fid,0,'eof');
filesize_in_bytes = ftell(fid);
total_records = (filesize_in_bytes - first_record_size)/(record_size*bytes_per_word);

%position file pointer to second record that contains real data, first record contains
%setup file in string format
fseek(fid, first_record_size,'bof');

if (total_records - floor(total_records) ~= 0)
    display(['WARNING: The file ' file ' does not contain an integer number of records']);
end

header = zeros(total_records, header_size);

% Read the headers from the file
for index = 1:total_records            % Start of Main loop

   junk = fread (fid, header_size,'ushort');% read the header
   header(index,:) = junk';
      
   %At this point a header has been read successfully, skip to next header

   fseek(fid,data_size*bytes_per_word,'cof');
   blocks = index;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We now have all headers for this file. So look for failed special character check flag.

fix_manually = []; % vector to track bad records with very large slips. Happens rarely
bad_record = find(header(:,16) ~= 0);% index to all failed special character checks

if (isempty(bad_record))
   display(['No bad records were found in file ' file]);
   return
end

for index = 1: length(bad_record)
    [blocks, junk,data]=read_tomi(fid,bad_record(index),bad_record(index));% get the bad record
    i=find(data==sp_char);
    if isempty(i)
        fix_manually=[fix_manually; bad_record(index)];
    else
        junk=find(diff(i)~=scan_length); %find instances where sp_char spacing incorrect (i.e. not at beginning of scan)
        if data(1)~=sp_char % slip occured in last scan of previous record.
            shift=scan_length+1-i(1); %need to shift by number of samples skipped
            if shift>=0              % If the shift is reasonable, then proceed to fix the data
                data(1+shift:data_size)=data(1:data_size-shift); %shift record
                data(1:shift)=data(scan_length+1:scan_length+shift); %patch in missed samples by copying next scan
                %               (usually just special character)
                fseek(fid,-bytes_per_word*data_size,'cof'); %move to start of data for this record
                fwrite(fid,data,'short'); % write it back into file
                % now must fix end of previous record
                [blocks, junk,data]=read_tomi(fid,bad_record(index)-1,bad_record(index)-1);
                data = reshape(data,scan_length,data_size/scan_length); % one scan per column.
                data(:,size(data,2)) = data(:,size(data,2)-1);%copy second last scan over last scan
                fseek(fid,-bytes_per_word*data_size,'cof'); %move to start of data for this record
                fwrite(fid,data,'short'); % write it back into file
                % reopen record to check for further slips
                [blocks, junk,data]=read_tomi(fid,bad_record(index),bad_record(index));
                i=find(data==sp_char);
                junk=find(diff(i)~=scan_length);
            else                     % If the shift is negative, then there is probably a problem, and will have to be fixed manually
                junk=[];
                fix_manually=[fix_manually; bad_record(index)];
            end
        end
    end
    
   
   while ~isempty(junk)
        junk=junk(1);
        if (i(junk+1)-i(junk))<scan_length&(i(junk+1)-i(junk))>(scan_length-5) %slips of more than 3 very unlikely
            shift=scan_length-i(junk+1)+i(junk); %shift by number slipped
            data(i(junk+1)+shift:data_size)=data(i(junk+1):data_size-shift); %shift record
            if junk~=1
                data(i(junk):i(junk)+scan_length-1)=data(i(junk)-scan_length:i(junk)-1);%copy previous scan into slipped scan
            else
                data(i(junk):i(junk)+scan_length-1)=data(i(junk)+scan_length:i(junk)+2*scan_length-1);%copy following scan into slipped scan
            end
            i=find(data==sp_char);
            junk=find(diff(i)~=scan_length);
        else % if slip appears to be large, check if perhaps a special character has shown up in the data
            data = reshape(data,scan_length,data_size/scan_length); % one scan per column
            n = find(data(1,:) ~=sp_char);% index to failed special characters
            data=data(:);
            if isempty(n)
                junk=[];
            elseif (n(1)-1)>junk %there wasn't a slip in that scan - simply a special character in the data
                i=i([1:junk junk+2:end]);
                junk=find(diff(i)~=scan_length);
            else %something odd is going on here and needs to be checked out manually
                junk=[];
                fix_manually=[fix_manually; bad_record(index)];
            end
        end
    end
 
    fseek(fid,-bytes_per_word*data_size,'cof'); %move to start of data for this record
    fwrite(fid,data,'short'); % write it back into file
    
end
fclose(fid);

