%% setupstr
% Read attributes of a configuration string from a RSI raw binary data file.
%%
% <latex>\index{Functions!setupstr}</latex>
%
%%% Syntax
%
%   varargout = setupstr( varargin )
%
% * [varargin] see Usage
% * [varargout] see Usage
%
%%% Description
%
% This function will extract requested attributes from an input string. The
% string should be a configuration string that was extracted out of a RSI
% data file. This string can also be obtained by reading a configuration
% file (.cfg). 
%
%%% Usage
%
% The setupstr function is called with up to four arguments.  Depending on the
% number and type of arguments used, different results are returned.
%
% The first argument is the configuration string, either extracted from a data
% file or loaded from a configuration file.  When this is the only argument, the
% $\texttt{setupstr}$ function returns an indexed structure containing the contents of the
% configuration string. This indexed structure can be used in place of the 
% configuration string, for subsequent function calls, to greatly increase
% the speed of this function.
%
%    >> obj = SETUPSTR('config_string')
%
% Parse the string $\texttt{'config\_string'}$ and return $\texttt{obj}$, a
% data structure containing the values in $\texttt{config\_string}$. The
% returned $\texttt{obj}$ can be used in subsequent calls to the function. 
%
%    >> S = SETUPSTR( obj, 'section' )
%
% Find sections within $\texttt{obj}$ that match the value
% $\texttt{'section'}$. See Notes for additional information regarding
% search queries. The returned value $\texttt{S}$ is a cell array
% containing the matching section identifiers in string format.
%
%    >> V = SETUPSTR( obj, 'section', 'parameter' )
%
% Find parameter values within $\texttt{obj}$ that are located in
% $\texttt{'section'}$ and have the parameter name $\texttt{'parameter'}$.
% See Notes for additional information regarding search queries. The
% returned value $\texttt{V}$ is a cell array containing the matching
% string values.
%
%    >> [S, P, V] = SETUPSTR( obj, 'section', 'parameter', 'value')
%
% Find values within $\texttt{obj}$ that are located in
% $\texttt{'section'}$, have the parameter name $\texttt{'parameter'}$, and
% the value $\texttt{'value'}$. See Notes for additional information
% regarding search queries. The returned values are a tuple of cell arrays
% containing the matching sections, parameter names, and parameter values.
%
% * [config_string] Configuration string.
% * [obj]   A previously returned structure containing indexed values
%           from a configuration string. 
% * [section] Search query for a matching section identifier.
% * [parameter] Search query for a matching parameter name.
% * [value]   Search query for a matching parameter value.
% * []
% * [obj] Structure containing indexed values from a configuration string,
%         that can be used in subsequent calls to this function to speed
%         the search.
% * [S]   Matching sections as an array of cells.
% * [P]   Matching parameter names as an array of cells.
% * [V]   Matching parameter values as an array of cells.
%
%%% Notes
%
% The $\texttt{[root]}$ section in a configuration file is usually not
% declared explicitly. It exists implicitly and consists of all content
% before the first explictly declared section. To access the parameters
% within the $\texttt{root}$ section, use $\texttt{'root'}$ for the section
% identifier.
%
%
%%% Examples
%
%    >> cfg = setupstr( setupfilestr );
%    >> value = setupstr( cfg, 'P', 'coef0' )
%
% From the variable $\texttt{setupfilestr}$, which is a string containing
% the entire configuration file used for data acquisition, query the value
% of the parameter $\texttt{'coef0'}$ for the pressure channel. 
%
%    >> sections = setupstr( setupfilestr, '' )
%
% Find all sections within a configuration file.
%
%    >> cfg = setupstr( setupfilestr );
%    >> [S,P,V] = setupstr( cfg, '', 'id', '' )
%
% Find all sections that have a parameter with name $\texttt{id}$, and
% place the section names into the cell-array $\texttt{S}$. The values of
% the $\texttt{id}$-parameters are placed into $\texttt{V}$, and $\texttt{P}$
% is filled with the value $\texttt{'id'}$. $\texttt{S}$,  $\texttt{P}$, and
% $\texttt{V}$ have identical length.

% Version History:
% 
% * 2012-04-30 (WID) initial function
% * 2012-10-24 (WID) Updated documentation
% * 2012-11-22 (WID) added support for arg1 being a configuration string or a
%                    structure.
% * 2012-11-22 (WID) updated documentation to match the manual
% * 2012-02-26 (WID) updated documentation - incorrect pressure example
% * 2013-05-23 (WID) modified so that no-fast, no-slow, and num_rows are
%                    no longer required.
% * 2014-01-21 (WID) corrected bug with single row matrix
% * 2015-11-13 (WID) remove channels that are not found within the matrix.
% * 2015-11-19 (RGL) Corrected documentation.
% * 2015-11-19 (WID) Tweak channels of type "accel" into channels of type
%                    "piezo" if required.
% * 2016-06-21 (WID) Fix crash when a section has no parameter/value pairs.
% * 2016-12-29 (WID) Added fix for XMP therm channel with diff_gain.
% * 2018-08-28 (JMM) Added a warning message if info is removed from setup
%                    string
% * 2020-02-11 (WID) Allow for "recordDuration" in place/addition to
%                    "recsize".  Set to default value of 1 if not defined.
% * 2020-02-12 (WID) Fix bug introduced with last update.

function varargout = setupstr(arg1, arg2, arg3, arg4)

if nargin == 0
  display_help();
  varargout{1} = [];      % default return value prevents error during help.
end


if nargin == 1
  if isstruct(arg1), varargout{1} = arg1; return; end
  arg1 = tweak_setup_structure( parse_setup_file( arg1 ) );
  varargout{1} = arg1;
end


% Assume a section expression was provided, return matching section.
if nargin == 2
  if isa( arg1, 'char' )
    arg1 = tweak_setup_structure( parse_setup_file( arg1 ) );
  end
  out = {};
  for i = section_indexes(arg1, arg2)
    out{end+1} = arg1(i).k{1};
  end
  varargout{1} = out;
end


% Assume a section and key were provided, return all matching values.
if nargin == 3
  if isa( arg1, 'char' )
    arg1 = tweak_setup_structure( parse_setup_file( arg1 ) );
  end
  out = {};
  for i = section_indexes(arg1, arg2)
    for ii = key_indexes(arg1, i, arg3)
      out{end+1} = arg1(i).v(ii).v;
    end
  end
  varargout{1} = out;
end


% Assume a section, key, and value were provided - return all matching truples.
if nargin == 4
  if isa( arg1, 'char' )
    arg1 = tweak_setup_structure( parse_setup_file( arg1 ) );
  end
  if isempty(arg4), value = '.*'; else value = ['^' arg4 '$']; end
  sections = [];
  keys = [];
  values = [];
  for i = section_indexes(arg1, arg2)
    for ii = key_indexes(arg1, i, arg3)
      if ~isempty(regexpi(char(arg1(i).v(ii).v), value, 'match'))
        sections{end+1} = arg1(i).k{1};
        keys{end+1} = arg1(i).v(ii).k;
        values{end+1} = arg1(i).v(ii).v;
      end
    end    
  end
  varargout{1} = sections;
  varargout{2} = keys;
  varargout{3} = values;
end

end



function display_help ()
disp(' ');
disp('SETUPSTR:  parse / query ODAS setup strings');
disp(' ');
disp('  OBJ = setupstr (setup_string)');
disp('    Parse setup_string into an array of struct elements.  The returned');
disp('    array should be passed into subsequent function calls.');
disp(' ');
disp('  S = setupstr (struct, section)');
disp('    Search for matching sections in the struct array.  Standard regex');
disp('    wildcards are accepted.  A cell array of matching sections is');
disp('    returned.');
disp(' ');
disp('  V = setupstr (struct, section, parameter)');
disp('    Search for matching section/parameter pairs in the struct array and');
disp('    return a cell array of associated values.  Standard regex wildcards');
disp('    are accepted.');
disp(' ');
disp('  [S, K, V] = setupstr (struct, section, parameter, value)');
disp('    Search for section, parameter, and value entries.  Standard regex');
disp('    wildcards are accepted.  Cell arrays with the matching sections,');
disp('    parameter, and values are returned.');
end



function out = parse_setup_file (setupFile)

out = [];
section = 'root';
s_index = 1;
index = 1;
dest = 1;

lines = strtrim(regexp(setupFile, '\n', 'split'));

for i=1:length(lines)
  
  % remove comments
  [a,b] = regexp(lines{i}, '^(.*?);.*', 'tokens');
  if ~isempty(b)
	  if ~isempty(a{1}), lines{i} = a{1}{1}; end
  end
  
  % find sections
  [a,b] = regexp(lines{i}, '^\s*\[\s*(.+)\s*\]', 'tokens');
  if ~isempty(b)
    section = lower(a{1}{1});
    index = 1;                    % reset the index..
    dest = 1;
    s_index = s_index + 1;
    continue;
  end
  
  % find assignments
  expr = '^\s*(.+?)\s*=\s*(.+?)\s*$';
  [a,b] = regexp(lines{i}, expr, 'tokens');
  if ~isempty(b)
    
    % check if the key is equal to "name" - if so create a new section 
    % using the value of "name" if it doesn't match the current section.
    % Use 'dest' to ensure the first element in the array is default name.
    if strcmpi(a{1}{1}, 'name')
      out(s_index).k(1) = lower(a{1}(2));
      dest = 2;
    end
    
    % we found an entry - time to add it.
    out(s_index).k(dest) = {section};
    out(s_index).v(index).k = lower(a{1}{1});
    
    % When we add the entry value, we want to add the evaluated version of the
    % entry if available.  For example, "4 + 7" should result in the string "11"
    % being added as the value.
    try
      temp = sprintf('%d',eval_in_clean_workspace( a{1}{2} ));
    catch
      temp = a{1}{2};
    end
    out(s_index).v(index).v = temp;
    
    index = index + 1;
    continue;
  end
end
end


function out = section_indexes (input, expression)
% return list of index values of struct elements where the section matches
% the provided expression.

if isempty(expression)
  sect = '.*'; 
else
  if iscell(expression), expression = expression{1}; end
  sect = ['^' expression '$'];
end

out = [];
for i=1:length(input)
  for ii=1:length(input(i).k)
    result = regexpi(input(i).k{ii}, sect, 'match');
    if ~isempty(result), out = [out i]; break; end
  end
end

end


function out = key_indexes (input, section, expression)
% return list of index values at the given section where the key matches
% the provided expression.

if isempty(expression)
  key = '.*'; 
else
  if iscell(expression), expression = expression{1}; end
  key = ['^' expression '$'];
end

out = [];
for i=1:length(input(section).v)
  result = regexpi(input(section).v(i).k, key, 'match');
  if ~isempty(result), out = [out i]; end
end

end



% Make changes to the configuration structure right before it is returned.
function clean = tweak_setup_structure( dirty )

    % The returned value will be a slightly modified version of the
    % original value.
    clean = dirty;

    % If no root entries are available, the following code does not work 
    % correctly.  Add a blank root section if root is missing.
    root_processed = false;
    for x = 1:length(clean)
        % Find the root section.
        if isempty(clean(x).k), continue; end
        if ~strcmp(clean(x).k{1},'root'), continue; end
        root_processed = true;
    end
    if ~root_processed
        clean(end+1).k{1} = 'root';
    end

    fast = setupstr( clean, 'root', 'no-fast' );
    slow = setupstr( clean, 'root', 'no-slow' );
    rows = setupstr( clean, 'matrix', 'num_rows' );
    
    % Inject values for no-fast, no-slow, and num_rows.
    matrix = [];
    for row = setupstr( clean, 'matrix', 'row[0-9]+' )
        % Generate the matrix found with the configuration file.
        C = textscan(row{1}, '%d');
        matrix(end+1,:) = C{1}';
    end
    
    % Purge the data structure of sections with no parameter / value pairs.
    for ch = setupstr( clean, '', '', '' )
        idx = 1;
        while idx <= length(clean)
            if isempty(clean(idx).k)
                clean = [clean(1:idx-1) clean(idx+1:end)];
            end
            idx = idx + 1;
        end
    end
    
    % If required, insert the number of rows.
    if isempty(rows)
        rows = size(matrix, 1);
        for x = 1:length(clean)
            if strcmp(clean(x).k{1},'matrix')
                clean(x).v(end+1).k = 'num_rows';
                clean(x).v(end).v = num2str(rows);
                break;
            end
        end
    end
    
    % If required, insert the number of fast and slow columns.
    diffmatrix = sum(abs(diff(matrix)),1);
    newslow = length(find(diffmatrix ~= 0));
    newfast = length(find(diffmatrix == 0));

    % Find the correct section but only instert the values if they are not
    % already there.
    for x = 1:length(clean)
        % Find the root section.
        if ~strcmp(clean(x).k{1},'root'), continue; end
        clean(x).v(end+1).k = 'num_fast';
        clean(x).v(end).v = num2str(newfast);
        if isempty(fast)
            clean(x).v(end+1).k = 'no-fast';
            clean(x).v(end).v = num2str(newfast);
        end
        clean(x).v(end+1).k = 'num_slow';
        clean(x).v(end).v = num2str(newslow);
        if isempty(slow)
            clean(x).v(end+1).k = 'no-slow';
            clean(x).v(end).v = num2str(newslow);
        end
        break;
    end
    
    % Supplement id_even and id_odd with an id parameter.    
    previous = '';
    for name = setupstr( clean, '', 'id|id_even|id_odd', '' )
        ch = char(name);
        
        % Process each section only once.  They will be in sequential
        % order.
        if strcmpi(ch, previous), continue; else, previous = ch; end
        
        % Find all declared id values for the section.
        id     = char(setupstr(clean, ch, 'id'));
        ideven = char(setupstr(clean, ch, 'id_even'));
        idodd  = char(setupstr(clean, ch, 'id_odd'));        
        
        if ~isempty(id)
            % ID is declared so calculate id_even and id_odd if required.
            tmpid = eval(['[' id(1,:) ']']);
            if ~isnumeric(tmpid), continue; end
            if length(tmpid) == 1, continue; end
            tmpid = sort(tmpid);
                
            %%%% Add the id_even and id_odd values here.....
            for x = 1:length(clean)
                % Find the correct section.
                if ~strcmp(clean(x).k{1}, ch), continue; end
                clean(x).v(end+1).k = 'id_even';
                clean(x).v(end).v = num2str(tmpid(1));
                clean(x).v(end+1).k = 'id_odd';
                clean(x).v(end).v = num2str(tmpid(2));
                break;
            end
        
        elseif ~isempty(ideven) && ~isempty(idodd)
            % Both id_even and id_odd required for calculating new value of
            % "id".

            %%%% Add id value here......
            for x = 1:length(clean)
                % Find the correct section.
                if ~strcmp(clean(x).k{1}, ch), continue; end
                clean(x).v(end+1).k = 'id';
                clean(x).v(end).v = sprintf('%s,%s', ideven, idodd);
                break;
            end
        end
    end

    % Purge channels that are not found within the address matrix.
    for ch = setupstr( clean, '' )
        if isempty(setupstr( clean, ch{1}, 'name' )), continue; end
        if isempty(setupstr( clean, ch{1}, 'type' )), continue; end
        ids = setupstr( clean, ch{1}, 'id' );
        
        for id = str2num(ids{1})  % do not use str2double - it causes errors
            if ~any(id == matrix)
                warning('Removing channel "%s" (ID %d not found within address matrix).', ch{1}, id)
                for x = length(clean):-1:1
                    if ~any(strcmp(clean(x).k,ch)), continue; end
                    % Remove channels that are not in the address matrix.
                    clean = [clean(1:x-1) clean(x+1:end)];
                end
            end
        end
    end
    
    % Change type "accel" to "piezo" for configuration files made before
    % type "piezo" existed.  "coef[0|1]" used to differentiate between the
    % two different types.
    for ch = setupstr( clean, '', 'type', 'accel' )
        coef0 = setupstr( clean, ch{1}, 'coef0' );
        coef1 = setupstr( clean, ch{1}, 'coef1' );
        if isempty(coef0) || isempty(coef1) || ~strcmp(coef0, '0') || ~strcmp(coef1, '1')
            continue;
        end
        
        for x = 1:length(clean)
            if ~strcmp(clean(x).k{1},ch{1}), continue; end
            
            for xx = 1:length(clean(x).v)
                if ~strcmp(clean(x).v(xx).k,'type'), continue; end
                clean(x).v(xx).v = 'piezo';
                break;
            end
            break;
        end
    end
    
    % Correct XMP configuration file error.  diff_gain added to channels
    % without pre-emphasis.  Causes problems for odas_p2mat so this
    % parameter should be removed.
    if ~isempty(setupstr( clean, 'T1', 'type', 'xmp_therm' ))
        for x = 1:length(clean)
            if ~strcmp(clean(x).k{1}, 't1'), continue; end
            for xx = 1:length(clean(x).v)
                if ~strcmp(clean(x).v(xx).k,'diff_gain'), continue; end
                clean(x).v(xx) = [];
                break;
            end
            break;
        end        
    end
    
    % Allow for "recordDuration" in place of "recsize".  Either parameter
    % can be used in a configuration file.  Future code should make use of
    % recordDuration to avoid the ambiguity associated with recsize.  This
    % section ensures that code written for either name will continue to
    % work.
    %
    % Note, this section will also add a default value of 1 when both
    % recsize and recordDuration are missing.
    recordDuration = setupstr( clean, 'root', 'recordduration' );
    recsize        = setupstr( clean, 'root', 'recsize' );
    
    if isempty(recordDuration) && isempty(recsize)
        recordDuration = 1;
        recsize = 1;
    elseif isempty(recordDuration)
        recordDuration = recsize{1};
    elseif isempty(recsize)
        recsize = recordDuration{1};
    elseif recordDuration ~= recsize
        error(['Both "recsize" and "recordDuration" are defined in the '
            'configuration file but they differ.  Only define one '
            'value or remove both values to accept the default value '
            'of "1".']);
    end    
    
    for x = 1:length(clean)
        % Find the root section.
        if ~strcmp(clean(x).k{1},'root'), continue; end

        recsize_done = 0;
        recordDuration_done = 0;
        
        % Test each value in the root section. Note when value already
        % present.
        for i = 1:length(clean(x).v)
            if strcmp(clean(x).v(i).k, 'recordduration')
                recordDuration_done = 1;
            end
            if strcmp(clean(x).v(i).k, 'recsize')
                recsize_done = 1;
            end
        end
        % Add values that were not found.
        if (~recsize_done)
            clean(x).v(end+1).k = 'recsize';
            clean(x).v(end).v = num2str(recsize);
        end
        if (~recordDuration_done)
            clean(x).v(end+1).k = 'recordduration';
            clean(x).v(end).v = num2str(recordDuration);
        end
        break;
    end
    
end

% Evaluate string in a seperate workspace to avoid conflicts with local
% variables.
function clean = eval_in_clean_workspace( statement )
    clean = eval( statement );
end

