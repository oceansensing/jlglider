%% show_P
% Extract the record-average pressure from a RSI raw binary data file.
%%
% <latex>\index{Functions!show\_P}</latex>
%
%%% Syntax
%
%   [P, record] = show_P( fileName )
%
% * [fileName] - String containing the name of a RSI raw binary data file.
% * []
% * [P] - Vector of the record-average pressure in units of dbar.
% * [record] - Vector of the record numbers corresponding to P.
%
%%% Description
%
% Calculate the record-average pressure from a RSI raw binary data file.
% This function is used to quickly glimpse the pressure-time history of an
% instrument. The pressure channel address must be present in the address
% $\texttt{[matrix]}$ and there must be a $\texttt{[channel]}$ section
% containing the coefficients for converting raw pressure data into
% physical units.
%
% The returned vectors are suitable for plotting the pressure history in a
% data file.
%
%%% Examples
%
%     >> [P records] = show_P( 'my_data_file.p' );
%     >> plot(P); set(gca, 'YDir', 'reverse'); grid on;
%
% Generate a pressure (depth) vs. record (time) line plot to reveal the
% descents and ascents of your instrument.
%
% Use
%
%     >> plot(diff(P)); grid on;
%
% to see the rate-of-change of preesure (vertical velocity).

% *Version History:*
%
% * 2010-01-15 (AWS) support for odas v6 and up
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-04-11 (WID) changed inifile_with_instring calls to setupstr
% * 2012-10-24 (WID) documentation update
% * 2013-02-26 (WID) fix name of pressure channel with setupstr
% * 2013-06-29 (WID) removed references to XMP instruments
% * 2015-11-02 (RGL) Documentation updates

function [P, record] = show_P(fileName)

% parsing of command line parameters
if nargin < 1; fileName    = []; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Get Input from the User
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(fileName); fileName = getDataFileName; end; % get data file name from the user
if (fileName == 'q'); return;end

% Look for the .p and .mat files, deal with problems
[P,N,E]=fileparts(fileName);
if isempty(E)
   fileName = [fileName '.p'];
elseif (strcmp(E, '.mat'))
    fileName = [N '.p']; % use the binary file
end

% Now find the pressure data (channel 10) in the file
[ch_slow,ch_fast] = query_odas(fileName);

address_matrix = [ch_slow, repmat(ch_fast, size(ch_slow,1), 1)];
address_matrix = address_matrix'; address_matrix = address_matrix(:);
P_index = find(address_matrix == 10);
P_index = P_index(1);
P_step = length(address_matrix);

% We now have all the information needed to extract the pressure data

fid = fopen_odas(fileName,'r'); % safe to assume that file exists.
junk = fread (fid,64,'ushort');
header_size_i    = 18;
record_size_i    = 19;
header_version_i = 11;
setupfile_size_i = 12;
header_size = junk(header_size_i)/2; % in units of 2-byte words
record_size = junk(record_size_i)/2;

data_size = record_size-header_size;

header_version = bitshift(junk(header_version_i), -8) + bitand(junk(header_version_i), 255) /1000;

if (header_version >=6)
    setupfile_size = junk(setupfile_size_i);

    %sanity check to see if there is a valid setup file size in header
    if (setupfile_size <= 0)
        error('incorrect setup file size extracted from first data record');
    end

    first_record_size = header_size * 2 + setupfile_size;
    setupfilestr = char(fread(fid, setupfile_size))';
else
    first_record_size = record_size*2;
end

fseek(fid, 0, 'eof'); %move to end of file
length_in_bytes = ftell(fid);
total_records = floor((length_in_bytes - first_record_size) / (record_size*2)); % use floor in case their is a rundat the end.

P = zeros(total_records,1); record = P; % pre-assign vectors

fseek(fid, 0, 'bof'); % move to the beginning of the file
fseek(fid, first_record_size, 'bof'); % move to the second record.

for index = (1:total_records)
    fseek(fid, 2*header_size, 'cof'); % skip over the header
    data = fread(fid, data_size, 'short'); 
    P(index) = mean(data(P_index: P_step: end)); % record-average pressure
    record(index) = index;
end

% Now convert into physical units.
if header_version >= 6
    P = convert_odas(P,'P', 'string', setupfilestr, header_version);
else
    P = convert_odas(P,'pres','file','setup.txt');
end

function fileName = getDataFileName;
% fileName = getDataFileName returns the file name of the data file.
% Preforms error checking.
%
% Fab, March 1998.
% Modified by RGL to offer the latest '*.p' file for opening.
% 2004-06--6


fileName='';
while isempty(fileName); 
    test_string = get_latest_file;
   fileName = input(['Enter data file name (default: ' test_string ' , ''q'' to quit): '], 's');
   if strcmp(fileName, 'q')
      fclose('all');return;
   elseif isempty(fileName)
      fileName = test_string;
   end
   if ~exist(fileName, 'file');
      fprintf('Can''t open file %s! Try again.\n', fileName);
      fileName = '';
   end
end

