%% patch_setupstr
% Patch a configuration string into an existing RSI raw binary data file
%%
% <latex>\index{Functions!patch\_setupstr}</latex>
%
%%% Syntax
%   patch_setupstr( rawDataFile, configFile, ['-force'], ['-revert'] )
%
% * [rawDataFile] Name of the data file into which 'configFile' should
%         be embedded.
% * [configFile] Name of the configuration file to embed.
% * ['-force']  Force the patch.  Use with extreme caution.
% * ['-revert'] Revert to the original configuration file. All previously 
%         applied patches will be lost.
%
%%% Description
%
% Patch an external configuration string into the first record of a RSI raw
% binary data file.
%
% This function is used to modify the configuration string in a data file.
% The configuration file may have erroneous conversion parameters that must
% be corrected, before the file is turned into a mat-file in physical
% units. It can also be used to upgrade an older-style data file that does not
% have a configuration string, so that it can be used with this Matlab
% Library.
%
% When updating a data file with a new configuration string, the
% original configuration string is commented out and appended to the new
% configuration string. This retains the contents of the original string,
% making it possible to revert to the original one.
%
% When patching a data file, the parameters that directly affect data 
% acquisition *MUST NOT* be changed.  These values include;
%
% # rate,
% # recsize,
% # no-fast,
% # no-slow,
% # matrix.
%
% The $\texttt{'-force'}$ bypasses all safety checks. Use with caution.
%
% The $\texttt{'-revert'}$ option reverts the configuration string to its
% original value. 
%
%%% Examples
%
%    >> patch_setupstr( 'data_005.p', 'setup.cfg' );
%
% Embed the configuration string in the file $\texttt{setup.cfg}$ into the data file 
% $\texttt{data\_005.p}$.
%
%    >> extract_setupstr( 'data_005.p', 'data_005.cfg' );
%    >> edit data_005.cfg
%    >> patch_setupstr( 'data_005.p', 'data_005.cfg' );
%
% Example of how $\texttt{patch\_setupstr}$ can be used with $\texttt{extract\_setupstr}$ to
% modify the calibration coefficients embedded within a data file. It is
% assumed that changes will be made by the $\texttt{edit}$ command.
%
%    >> patch_setupstr( 'data*', '-revert' )
%
% Removes all patches to the specified data files. Note the use of a wild
% card character to specify more than a single file. 
%
% Older-style data files, that do not have a configuration string, can be
% upgraded with a configuration-string, so that the data file can be
% processed with this version of the ODAS Matlab Library. The upgrade
% replaces the first record with a new header and the configuration string.
% No data is lost because there is no data in the first record of
% older-style data files. 
%

% *Version History:*
%
% * 2012-03-14 (WID) initial version
% * 2012-03-22 (WID) extensive rewrite, bugs fixed + the following improvements
%                     - error checking when reading, writing, and opening files
% * 2012-03-26 (WID) added support for patching v1 data files with v6 headers
% * 2012-04-11 (WID) replaced inifile_with_instring calls with setupstr
% * 2012-04-25 (WID) updated documentation and added changes for R2007b
% * 2012-09-09 (WID) modified to utilize file_with_ext.  Updated documentation.
% * 2012-11-05 (WID) updated documentation
% * 2015-02-16 (WID) include the force option
% * 2015-07-17 (WID) added support for revert option.  Now saves the
%                    original configuration file within the comments.
% * 2015-07-27 (WID) documentaion update.
% * 2015-10-30 (RGL) documentaion update.
% * 2016-02-03 (WID) fix for older versions of Matlab.
% * 2018-10-26 (WID) fix for v1.0 data files


function patch_setupstr(datafname, varargin)

% Constants:
CONFIG_IDENTIFIER = '; ### Original Configuration String Below ###';

% Process input options...
% Currently only look for the '-force' option.
force = false;
revert = false;
setupfname = [];

for arg = varargin
    if strcmpi(arg, '-force') || strcmpi(arg, '-f')
        force = true;
        warning('Forcing Insertion of new configuration file.');
    elseif strcmpi(arg, '-revert') || strcmpi(arg, '-r')
        revert = true;
        disp('Reverting to the original configuration file.');
    else
        setupfname = arg{1};
    end
end


% Attempt to open the data file.
[N,M,E,datafname] = file_with_ext( datafname, ...
                       {'' '.p' '.P'}, ...
                       ['Unable to find file: ' datafname] );

[data_fid, msg] = fopen_odas(datafname, 'r');
if (data_fid < 3 || ~isempty(msg))
    error('Unable to open file: %s\n%s', datafname, msg);
end


% If configuration file not specified, it is an error - unless we are
% reverting to the original version.
if isempty(setupfname) && ~revert
    error('Configuation file not specified.');
    
elseif ~revert
    % Attempt to open the setup file.
    [N,M,E,setupfname] = file_with_ext( setupfname, ...
        {'' '.cfg' '.CFG'}, ...
        'ODAS v6 configuration file not found.' );
    
    [setup_fid, msg] = fopen(setupfname, 'r');
    if (setup_fid < 3 || ~isempty(msg))
        error('Unable to open file: %s\n%s', datafname, msg);
    end

end


% Read the header, process errors.
[header, size] = fread(data_fid, 64, '*ushort');
if size ~= 64, error(['Unable to read header of file: ' datafname]); end

% set endian format.
if (header(64) == 1), endian = 'l'; else, endian = 'b'; end

% set version of the data file.  Should be 0x0001 for legacy data files
if header(11) == 1, version = 1; else, version = header(11); end

% Construct new header for v1 data files.
if version == 1, header = convert_header(header); end

% The new generated data file will be created as a temporary file - make 
% sure this temp file doesn't exist - Matlab doesn't guarantee the function call
% 'tempname' returns a unique name.  Loop until one is found.
new_file = tempname;
while exist(new_file, 'file'), new_file = tempname; end
    
[patched_fid, msg] = fopen(new_file, 'w', endian);
if (patched_fid < 3 || ~isempty(msg))
    error('Unable to create temporary file: %s\n%s', new_file, msg);
end
    
% the size in bytes of the setup file string contained in the first record
setupFileSize = header(12);

% read the old and new setup file into single strings
[oldstring, size] = fread(data_fid, double(setupFileSize), '*char*1');
if size ~= setupFileSize, error('Unable to extract setup file.'); end

if ~revert
    [newstr, size] = fread(setup_fid, inf, '*char*1');
    if size == 0, error('Unable to read provided setup file.'); end
    newstr = newstr';
    
else
    % If reverting, we grap a copy of the old string, strip out the comment
    % marks, and use it as the input string.
    idx = regexp(oldstring', CONFIG_IDENTIFIER, 'split');
    
    if length(idx) == 2
        idx = regexp(idx{2}, '((\n\r)|(\r\n)|\n|\r); ', 'split');
        newstr = strjoin(idx, '\n');
    else
        error('Unable to extract original configuration string.');
    end
end

oldstr = oldstring';

% skip different safetly checks for v1 data files.  The observed setup file for 
% such files will be invalid - it's the first data block.
if version == 1 && ~force
    
  % check that the following values are the same:  rate, no-fast, no-slow, 
  % and the matrix
  cfg = setupstr(newstr);
  n_rate = str2double(setupstr(cfg,'root','rate'));
  n_fast = str2double(setupstr(cfg,'root','no-fast'));
  n_slow = str2double(setupstr(cfg,'root','no-slow'));
  n_rows = str2double(setupstr(cfg,'matrix','num_rows'));
  
  % Skip a record. The correct rate is stored in the second record.
  [header, size] = fread(data_fid, 64, '*ushort');
  
  % derive old rate from header...
  freq = header(21) + header(22)/1000;
  o_rate = freq / (header(29) + header(30));
  if abs(n_rate - o_rate) > 0.5
    error('\ndifference in rate detected: old = %d  new = %d', uint16(o_rate), n_rate);
  end
  
  % compare no-fast
  if n_fast ~= header(29)
    error('\ndifference in no-fast detected: old = %d  new = %d', header(29), n_fast);
  end
  
  % compare no-slow
  if n_slow ~= header(30)
    error('\ndifference in no-slow detected: old = %d  new = %d', header(30), n_slow);
  end
  
  % read in the existing matrix.  It is stored as characters in the oldstr 
  % variable.  Remember, endian format matters...
  c = header(29) + header(30);          % c for columns
  r = header(31);                       % r for rows
  oldstr = oldstr(1:c*r*2);             % trim to matrix length (in bytes)
  oldstr = cast(oldstr, 'uint8');       % we want integers, not characters
  oldstr = reshape(oldstr, 2, c*r);     % prepare for conversion into words
  if endian == 'l', msb = 2; lsb = 1; else msb = 1; lsb = 2; end
  matrix = oldstr(msb,:).*256 + oldstr(lsb,:);    % convert bytes to words
  matrix = reshape(matrix, c, r);                 % turn into matrix
  
  % compare the number of rows
  if n_rows ~= r
    error('Invalid number of matrix rows, old: %d new: %d', r, n_rows);
  end
  
  % now read the matrix from the setup file and compare
  rows = setupstr(cfg, 'matrix', 'row[0-9]+');
  
  for r=1:n_rows
    % parse row into integers
    C = textscan(rows{r}, '%d16');

    % compare the number of elements in each row
    if c ~= length(C{1})
      error('Matrix row #%d has wrong number of columns.', r);
    end

    % compare matrix rows.  Note that the matrix is the transpose of what you
    % would expect.
    index = find(C{1} ~= matrix(:,r), 1);
    if ~isempty(index)
      error('Matrix(%d,%d) values differ - old/new = %d/%d', r,index, matrix(index,r), C{1}(index));
    end
  end

elseif ~force
  % This is a v6 (or larger) data file

  cfg_new = setupstr(newstr);
  cfg_old = setupstr(oldstr);

  n_rate = str2double(setupstr(cfg_new,'root','rate'));
  n_recs = str2double(setupstr(cfg_new,'root','recsize'));
  n_fast = str2double(setupstr(cfg_new,'root','no-fast'));
  n_slow = str2double(setupstr(cfg_new,'root','no-slow'));
  n_rows = str2double(setupstr(cfg_new,'matrix','num_rows'));

  o_rate = str2double(setupstr(cfg_old,'root','rate'));
  o_recs = str2double(setupstr(cfg_old,'root','recsize'));
  o_fast = str2double(setupstr(cfg_old,'root','no-fast'));
  o_slow = str2double(setupstr(cfg_old,'root','no-slow'));
  o_rows = str2double(setupstr(cfg_old,'matrix','num_rows'));

  if n_rate ~= o_rate
    error('Difference in "rate" detected: old=%d  new=%d\n',o_rate,n_rate);
  end  
  if n_recs ~= o_recs
    error('Difference in "recsize" detected: old=%d  new=%d\n',o_recs,n_recs);
  end  
  if n_fast ~= o_fast
    error('Difference in "no-fast" detected: old=%d  new=%d\n',o_fast,n_fast);
  end
  if n_slow ~= o_slow
    error('Difference in "no-slow" detected: old=%d  new=%d\n',o_slow,n_slow);
  end
  if n_rows ~= o_rows
    error('Difference in "num_rows" detected: old=%d  new=%d\n',o_rows,n_rows);
  end
    
  % now make sure that the matrix rows are identical 
  n_rows = setupstr(cfg_new,'matrix','row[0-9]+');
  o_rows = setupstr(cfg_old,'matrix','row[0-9]+');

  for i=1:length(n_rows)
    n_row = textscan(n_rows{i}, '%f');
    o_row = textscan(o_rows{i}, '%f');
    if find(n_row{1} ~= o_row{1})
      error('Matrix differs: #%d, old_row=%s  new_row=%s\n',i,char(o_row),char(n_row));
    end
  end
end     % end safety checks

% Find the original configuration string from the file's embedded
% configuration string.  If not specified, convert the entire embedded
% configuration string into a new "original" configuration string.

if version ~= 1
    idx = strfind(oldstring', CONFIG_IDENTIFIER);
    if isempty(idx)
        result = '';
        parts_old = regexp(oldstring', '(\n\r)|(\r\n)|\n|\r', 'split');
        for part = parts_old
            result = sprintf('%s\n; %s', result, part{1});
        end
        oldstring = sprintf('%s%s', CONFIG_IDENTIFIER, result);
    else
        oldstring = oldstring(idx(end):end)';
    end
    
    % Anything below the configuration identifier has to be removed from the
    % input configuration string.
    idx = regexp(newstr, CONFIG_IDENTIFIER, 'split');
    newstr = idx{1};
else
    oldstring = '\n';
end

% show patch date/time at the top of the file to remind user that we are dealing with a patched setup file.
if ~revert
    tmpstr = ['; ' datestr(now,'yyyy-mm-dd HH:MM:SS') ' patched configuration string' char(10) ];
    newstr = [tmpstr newstr];
end

parts_new = regexp(newstr,    '(\n\r)|(\r\n)|\n|\r', 'split');
parts_old = regexp(oldstring, '(\n\r)|(\r\n)|\n|\r', 'split');

newstr = strjoin(parts_new, '\n');
oldstr = strjoin(parts_old, '\n');

if ~revert
    % Only use the old string when not reverting.
    newstr = strjoin( {newstr oldstr}, '\n\n' );
end

header(11) = 6 * 256;
header(12) = length(newstr);

% write the header, return error if not fully written
count = fwrite(patched_fid, header, 'ushort');
if count ~= 64, error('Unable to write header to new data file.'); end
   
% write new setup string, return error if not fully written
count = fwrite(patched_fid, newstr);
if count ~= length(newstr), error('Error writing setup into data file.'); end

% version 1 files now require the next header be written.
if version == 1
    count = fwrite(patched_fid, header, 'ushort');
    if count ~= 64, error('Unable to write header to new data file.'); end
end

ONEMB = 1024*1024;

% It is assumed that the file pointer is positioned at the second record
% before entering this loop.
while ~feof(data_fid)
    [buf, numread] = fread(data_fid, ONEMB);
    fwrite(patched_fid,buf(1:numread));
end

fclose('all');

% Move the new data file into the location of the old data file
[success, msg] = movefile(new_file, datafname);
if ~success, error('Error writing to datafile: %s\n%s', datafname, msg); end

end



function new = convert_header(header)

    % by default, copy all the values.
    new = header;
    
    % data file format to v6.0
    new(11) = bitshift(6, 8);
    
    % setup file size not yet known.  Set to size of data block so we skip over
    % the first empty block when extracting the setup file.
    new(12) = header(19) - header(18);
    
    % product ID is unknown - leave as legacy so we always know this was a 
    % converted file.
    new(13) = 0;
    
    % unknown build number - just leave as 0.
    new(14) = 0;
    
end


% For those with older versions of Matlab, this function is provided.  It 
% overrides the function provided by Matlab on newer versions of Matlab.
function result = strjoin(str, delim)

    result = str{1};
    
    if nargin == 1, delim = ' '; end
    if length(str) <= 1, return; end
    
    for i = 2:length(str)
        result = sprintf([texstr(result,'') delim texstr(str{i},'')]);
    end
    
end
