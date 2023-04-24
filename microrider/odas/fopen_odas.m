%% fopen_odas
% Open a RSI raw binary data file of unknown endian type
%%
% <latex>\index{Functions!fopen\_odas}</latex>
%
%%% Syntax
%   [fid, errorMsg] = fopen_odas(fileName, permission)
%
% * [fileName]   Name of the requested RSI raw data file.
% * [permission] Permission string, usually 'r', or 'r+'. Consult Matlab 
%                documentation for the command 'fopen'.
% * []
% * [fid] File pointer to `fileName'.  Values less then 3 indicate 	errors.  For
%         error values, consult Matlab documentation for the command 'fopen'.
% * [errorMsg] String indicating the endian flag error.  Empty if file failed to 
%         open or no error was observed.
%
%%% Description
% Used to open a RSI raw binary data file of unknown endian type.  This function
% observes the endian flag in the header to determin the correct endian format 
% of the file.
% 
% When the endian flag matches the endian format deduced by reading the file, 
% the file is opened without returning any error messages.
%
% When the endian flag is empty (equal to 0), we assume the file is in 'little' 
% endian format.  We open the file and return the appropriate error message.
% 
% Should the endian flag indicate the opposite endian value as that deduced by
% reading the header, the file is opened using the deduced value.  An error 
% message is returned indicating the error.
%
% All other values of the endian flag are considered invalid and result in the
% file not being opened and the appropriate error message being returned.
%
% Acceptable endian flags in the header are as follows:
%
% # 0 = unknown
% # 1 = little endian, the usual value with Intel and compatible processors
% # 2 = big endian, the usual value with ODAS5-IR
%
% Possible _error_msg_ Values:
%
% # 0: endian flag not found, assume 'little'.
% # 1: endian flag does not match file - ignore flag and assume 'little'.
% # 2: endian flag does not match file - ignore flag and assume 'big'.
% # 3: endian flag invalid - byte offsets (127,128) = |byte127, byte128|.
%
% Error string 0 occurs when the endian flag is zero.  When this happens we
% assume the file to be in 'little' endian format, open the file, and report the
% error.
%
% Error strings 1 and 2 occur when the endian flag in the header does not 
% agree with the observed endian format of the file.  In these scenarios we 
% open the file with the observed endian format and report the error.
%
% Error string 3 occurs when the endian flag is invalid.  This is most likely 
% the result of a corrupted or invalid data file.  The two bytes representing
% the endian flag are displayed.  The file is closed.
%

% *Version History:*
%
% * 2009-01-21 (RGL) initial
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-03-21 (WID) rewritten - simpler, added documentation.
% * 2012-05-09 (WID) simplified and fixed.
% * 2015-10-28 (RGL) Documentation changes.
% * 2016-01-15 (WID) Forced to open in "US-ASCII" format.

% =========================================================================
function [fid, error_msg] = fopen_odas(file_name, permission)

flagBEndian = 2;
flagLEndian = 1;

error_msg = [];

error_msg0 = '0: endian flag not found, assume ''little''.';
error_msg1 = ['1: endian flag does not match file - ignore flag and ' ...
			 'assume ''little''.'];
error_msg2 = ['2: endian flag does not match file - ignore flag and ' ...
			 'assume ''big''.'];
% Only a partial error message - values have to be appended if used.
error_msg3 = '3: endian flag invalid - big endian word = ';

% open the requested file.  Read header in both big and little endian formats.
fid = fopen(file_name, permission, 'b');
if (fid < 3), return, end % file open failed
headerB = fread(fid, 64, '*uint16');
fclose(fid);

fid = fopen(file_name, permission, 'l');
if (fid < 3), return, end % file open failed
headerL = fread(fid, 64, '*uint16');
fclose(fid);

% no endian flag was in the header, assume little endian.
if ~headerB(64)
	endian = 'l';
	error_msg = error_msg0;

% header correctly contained the 'big' endian flag
elseif headerB(64) == flagBEndian
	endian = 'b';

% header correctly contained the 'little' endian flag
elseif headerL(64) == flagLEndian
	endian = 'l';

% header contained 'big' endian but was of type 'little' - assume 'little' 
elseif headerL(64) == flagBEndian
	endian = 'l';
    error_msg = error_msg1;

% header contained 'little' endian but was of type 'big' - assume 'big' 
elseif headerB(64) == flagLEndian
	endian = 'b';
    error_msg = error_msg2;

% endian is invalid - set error msg, close file, set fid error.
else
  error_msg = [error_msg3 headerB(64) '.'];
  fid = -1;
  return
end

% Now re-open the file with the correct endian flag
fid = fopen(file_name, permission, endian, 'US-ASCII');

