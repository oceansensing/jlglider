%% query_odas
% Report which channels are in a RSI raw binary data file.
%%
% <latex>\index{Functions!query\_odas}</latex>
%
%%% Syntax
%
%   [ch_slow, ch_fast] = query_odas( fileName )
%
% * [fileName] RSI raw binary data file name (.p) with optional extension
% * []
% * [ch_slow] vector of slow channels from the address matrix 
% * [ch_fast] vector of fast channels from the address matrix 
%
%%% Description
%
% Return a list of channel $\texttt{id}$ numbers of the data in a RSI raw
% binary data file, a $\texttt{p}$-file. This function reads the address
% $\texttt{[matrix]}$ in the configuration string, and returns a list of the
% channels in the matrix. If two output variables are specified, then
% there is a separate list for slow and fast channels. If you specify a
% single output variable, then the slow and fast channels are listed
% together.
% 
%%% Examples
%
%    >> query_odas
%
% Extract the channel $\texttt{id}$ numbers from a data file. The user will
% be prompted for the file name.
%
%    >> ch_list = query_odas( 'myfile' );
%
% Channel $\texttt{id}$ numbers for all of the channels within the file
% $\texttt{myfile}$ are stored into the $\texttt{ch\_list}$ variable.
%
%    >> [ch_slow,ch_fast] = query_odas( 'myfile.p' );
%
% Channel $\texttt{id}$ numbers for the slow and fast channels within the
% file $\texttt{myfile}$ are listed separately in $\texttt{ch\_slow}$ and
% $\texttt{ch\_fast}$, respectively. 

% *Version History:*
%
% * 2005-02-08 IG  initial
% * 2009-03-06 RGL changed fopen to fopen_odas.
% * 2009-11-30 AWS adapted to handle odas v6 and up
% * 2011-09-01 AWS added documentation tags for matlab publishing
% * 2012-04-11 WID changed inifile_with_instring calls to setupstr
% * 2012-09-10 WID updated documentation for publishing
% * 2015-10-31 RGL Changed documentation.

function [ch_slow,ch_fast] = query_odas(fname)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Check the file name and open the file %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ask the user for a file name if not provided as an input variable
if nargin<1, fname = input('Enter the ODAS file name: ','s'); end

% Open the file
[P,N,E]=fileparts(fname);
if isempty(E)
   fname = [fname '.p'];
end
[fid, error_message] = fopen_odas(fname,'r');
if ~isempty(error_message), error_message, end;
if fid<3
   error(sprintf('Unable to open file %s\n',fname))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Constants & parameters %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
bytes_per_word = 2;     % Assume 2-byte integers
header_size_i = 18;     % Index of header_size in header record
block_size_i = 19;      % Index of block_size in header record
header_size_0 = 64;     % Predicted size of the header in each block (integers)
fast_cols_i = 29;       % Header indices for the number of fast & slow columns
slow_cols_i = 30;
setupfile_size_i = 12;  % Index of setup file size (header v6 and up)
header_version_i = 11;  % Index of header version, prior to odas v6, this is always 1
rows_i = 31;            % Header index for the number of rows

% Parameters read in from the first header in the file
fseek(fid,0,'bof');
HD = fread(fid,header_size_0,'ushort');
header_size = HD(header_size_i)/bytes_per_word;
record_size = HD(block_size_i)/bytes_per_word;
data_size = record_size - header_size;
fast_cols = HD(fast_cols_i);
slow_cols = HD(slow_cols_i);
n_cols = fast_cols + slow_cols;
n_rows = HD(rows_i);
n_slow = n_rows*slow_cols;

%MSB has major version, LSB has minor version
header_version = bitshift(HD(header_version_i), -8) + bitand(HD(header_version_i), 255) /1000;

if(header_version >= 6)
    setupfile_size = HD(setupfile_size_i);
    if setupfile_size <= 0
        fclose(fid);
        error('incorrect setup file size extracted from first data record');
    end
    setupfilestr = char(fread(fid,setupfile_size, 'char'))';
    if(isempty(setupfilestr))
        fclose(fid);
        error('failed to extract setup file string from first data record');
    end

    cfg = setupstr(setupfilestr);
    rows = setupstr(cfg, 'matrix', 'row[0-9]+');
    ch_matrix = [];
    for row = rows
      values = textscan(row{1}, '%d16');
      ch_matrix = vertcat(ch_matrix, values{1}');
    end
else
    d_buffer = fread(fid,data_size, 'short');
    ch_matrix = reshape(d_buffer(1:n_cols*n_rows), n_cols, n_rows)';
end

if fast_cols~=0
    ch_fast = ch_matrix(:,slow_cols+1:n_cols);
    ch_fast = ch_fast(1,:);
else
    ch_fast=[];
end
if slow_cols~=0
    ch_slow = ch_matrix(:,1:slow_cols)';
    ch_slow = ch_slow(:)';
    ch_slow = reshape(ch_slow,slow_cols,n_rows)';
else
    ch_slow=[];
end
clear ch_matrix
fclose(fid);

if nargout<=1, ch_slow = sort([ch_slow(:)' ch_fast]); end




