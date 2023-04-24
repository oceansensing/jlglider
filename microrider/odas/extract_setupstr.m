%% extract_setupstr
% Extract configuration string from a RSI raw binary data file.
%%
% <latex>\index{Functions!extract\_setupstr}</latex>
%
%%% Syntax
%   configStr = extract_setupstr( dataFileName, configFileName )
%
% * [dataFileName]   Name of a version raw binary data file (version >= 6).
%       (extension optional) 
% * [configFileName] Optional name of configuration file that should be
%        generated.  When excluded, no file is created but the
%        configuration string is still returned.
% * []
% * [configStr]      Resulting configuration string.
%
%%% Description
% Extract the configuration string from a RSI raw binary data file.  Use
% this function when one needs a copy of the configuration file used when
% the data file was created. Note that this function will only work on a v6
% or greater data file. 
%
% If the $\texttt{configFileName}$ is not provided and configStr not
% requested, this function will display the resulting configuration string.
% When $\texttt{configFileName}$ is provided, a new configuration file
% named $\texttt{configFileName}$ is created with the contents of the
% configuration string. When $\texttt{configStr}$ is requested, the
% variable $\texttt{configStr}$ is set to the configuration string. 
%
% Used in conjunction with $\texttt{patch\_setupstr}$, this functions allows one to
% modify the calibration coefficients for the conversion of the data into
% physical units, and for other forms of data processing. First, extract 
% the configuration string into a new file.  Edit the file as required.
% Finally, use $\texttt{patch\_setupstr}$ to update the data file with the new
% configuration file. 
%
%%% Examples
%
%    >> configFile = extract_setupstr( 'data_005.p' )
%
% Extract the configuration string from the data file $\texttt{data\_005.p}$ and store in
% the variable $\texttt{configFile}$.
%
%    >> extract_setupstr( 'data_005.p', 'setup.cfg' )
%
% Extract the configuration string from $\texttt{data\_005.p}$ and store it in the file
% $\texttt{setup.cfg}$.

% *Version History:*
%
% * 2012-03-13 WID initial version
% * 2012-09-03 WID made use of file_with_ext for finding the file
% * 2012-11-05 WID documentation update
% * 2012-11-27 WID fixed call to file_with_ext
% * 2012-11-27 WID replaced fclose(all) with proper file closing.
% * 2015-10-28 RGL Documentation changes.
% * 2015-11-18 RGL Documentation changes.

function out_file = extract_setupstr(data_file, setup_file)

if nargin == 0
    disp('extract_setupstr v1.0.0');
    disp('Copyright (c) 2012 RSI');
    disp('Type <help extract_setupstr> for a complete description.');
    return
elseif nargin == 1
	% Nothing to do here..
elseif nargin == 2
    if exist(setup_file,'file')
       error('File "%s" already exists, please select another file.', setup_file);
    end
     % Open/create the new data file...
    fdout = fopen(setup_file,'w+');
else
	disp('Invalid number of input parameters.');
    disp('Type <help extract_setupstr> to get more help on its usage.');
	error('invalid input parameters');
end

% Open the data file - .p extension not required.
[P N E data_file] = file_with_ext( data_file, ...
                                   {'.p' '.P' ''}, ...
                                   ['Unable to find file: ' data_file] );

% Finished parsing the data_file name, now open it.
[fdin, error_message] = fopen_odas(data_file,'r');
if ~isempty(error_message)
	error(error_message);
end
if fdin<3
	error('Unable to open file %s', data_file)
end

% Read the header - note that we ignore the 'endian' format of the file
% and deal with it as an array of 8 bit values.
header = fread(fdin, 128, '*uint8');
    
% Read 8 bits of the 'endian' flag, if 0 then use the other 8 bits
endian = header(127);
if (header(127) == 0)
	endian = header(128);
end
    
% Now test the endian format of the header so we read the correct size
% of the setup file string.
if (endian == 1)
    % The header format is in little endian....
    config_file_size = uint16(header(24)) * 256 + uint16(header(23));
else
    config_file_size = uint16(header(23)) * 256 + uint16(header(24));
end

% Read the embedded portion that is the setup file.
setup_file = fread(fdin, double(config_file_size), '*uchar');

display_result = 1;

% Write to the output file, if declared.
if exist('fdout', 'var')
	fwrite(fdout, setup_file, 'uchar');
    display_result = 0;
end

out_file = char(setup_file');
if ~ispc,
    out_file = regexprep(out_file, '\x{0D}\x{0A}', '\x{0A}');
end
if nargout == 0 && display_result,
    disp(out_file);
end

try
    fclose(fdin);
    fclose(fdout);
catch
end