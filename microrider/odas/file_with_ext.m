%% file_with_ext
% Find files with default extensions
%%
% <latex>\index{Functions!file\_with\_ext}</latex>
%
%%% Syntax
%   [P,N,E,fullName] = file_with_ext( file, extensions, errorStr )
%
% * [file] Name of file with or without extension.
% * [extensions] Cell array of acceptable extensions. ex, {'.p' '.P' ''}
% * [errorStr] Error string used if error is to be triggered.  (optional)
% * []
% * [P] Full path to the requested file.
% * [N] Name of requested file.
% * [E] Extension of requested file.
% * [fullName] Full path and file name with path seperators. (optional)
%
%%% Description
% Find files with default extensions. Convenience function that allows
% users to input file names without an extension.
%
% This function searches for a file with the name $\texttt{file}$ and an
% extension from the array $\texttt{extensions}$.  When found, it returns
% the file components $\texttt{[P,N,E]}$ along with the fully qualified
% name, $\texttt{fullName}$. If not found, the function triggers the error
% $\texttt{errorStr}$. Should $\texttt{errorStr}$ be empty then no error
% is triggered and empty values are returned for
% $\texttt{[P,N,E,fullName]}$. 
%
% The return value $\texttt{fullName}$ is provided to simplify addressing
% the resulting file. Alternatively, this value can be constructed from the
% $\texttt{[P,N,E]}$ values.
%
% The $\texttt{extensions}$ array lists the acceptable file extensions for
% the requested file. Each extension should consist of a string value
% prepended with a period ($\texttt{.}$). To allow for the extension to be
% declared within the $\texttt{file}$ variable, prepend the array with an
% empty string. See the example for details.
%
% Be sure to include both possible case values within $\texttt{extensions}$
% - for example, $\texttt{\{'.p', '.P'\}}$. Some operating systems are case
% sensitive while others are not.
%
%%% Examples
%
%    >> file_err = 'No valid input file was found..';
%    >> [P,N,E] = file_with_ext( file, {'' '.p' '.P'}, file_err );
%
% Typical usage where $\texttt{file}$ is a string provided by the user.
%
%    >> [P,N,E] = file_with_ext( file, {'' '.p' '.P'} );
%    >> if isempty(N), %processing error...;  end
%
% Errors are optional.  If $\texttt{errorStr}$ is not provided, this
% function will not trigger an error if the file is not found.
%
%    >> [P,N,E,fName] = file_with_ext( fName, 
%                                     {'' '.p' '.P'},
%                                     ['Unable to find file: ' fName] );
%    >> dos(['mv ' fName ' /dev/null']);
%
% Record the full matching file identifier in $\texttt{fName}$. This
% example demonstrates how the function can be utilized in a single call.
% The subsequent $\texttt{dos}$ command will effectively delete the file on
% a UNIX system. Using the full name and path ensures the wrong file does
% not get deleted.

% Version History
%
% * 2012-04-30 (WID) initial
% * 2012-05-14 (WID) updated documentation
% * 2012-07-19 (WID) fixed iteration over cell array + support path inputs
% * 2012-08-14 (WID) Documentation update + some simplified code
% * 2012-09-09 (WID) Bug fix for testing if files are present.
% * 2012-11-05 (WID) Documentation update.
% * 2014-04-15 (WID) Return '.' for P if P is empty - current directory.
% * 2015-01-26 (WID) Modified to return the absolute path in fullName.
% * 2015-01-29 (WID) Expand astrix characters in returned values.
% * 2015-04-27 (WID) Added check for empty file name / extensions.
%                    Return proper error identifiers.
% * 2015-10-28 (RGL) Documentation changes.
% * 2015-11-18 (RGL) Documentation changes.

function [P,N,E,fullName] = file_with_ext(file, extensions, error_str)

if nargin < 2
  error('ODAS:argCount', 'Invalid number of arguments.');
end

if nargin == 2
  error_str = [];
end

if isempty(file) || isempty(extensions)
  error('ODAS:argEmpty', 'Empty input arguments.');
end

% Open the 2-byte binary data file for reading
[P,N,E] = fileparts(file);

% Test extensions - in order from first to last.
for ext = extensions
    
  % Generate file name using one of the possible extensions.
  fullName = fullfile(P, [N E char(ext)]);
  
  % Test if the file name points to a real file.
  [found results] = fileattrib(fullName);
  
  % If found, assign output variables then exit.
  if found
      [P,N,E] = fileparts(fullName);
      % If P not specified it is empty.  Set to current directory.
      if isempty(P), P = '.'; end
      % Extract absolute path generated from the "fileattrib" function.
      fullName = results.Name;
      % Expand out wild characters -- like "*"
      [jnk,N,E] = fileparts(fullName);
      return
  end
  
end

% No extensions were found.  Return an error if requested.
P = ''; N = ''; E = ''; fullName = '';
if ~isempty(error_str), error('ODAS:fileNotFound', error_str); end

end
