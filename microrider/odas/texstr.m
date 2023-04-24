%% texstr
% Format a string for use by a TeX interpreter
%%
% <latex>\index{Functions!texstr}</latex>
%
%%% Syntax
%   result = texstr( input, escChar )
%
% * [input] String, cell, or cell array to reformat.
% * [escChars] String containing characters to interpret literally.
%         Default is `_\'.
% * []
% * [result] Input with all specified characters prepended with a forward
%         slash `\'. Output type is the same as the input type -- string,
%         cell, or cell array.
%
%%% Description
% Parse the input string and escape specified characters with a
% forward slash `\'.  This allows the strings to be used with third
% party applications such as LaTeX and was originally written to 
% facilitate underscore characters within Matlab plot labels. 

% Version History
%
% * 2013-07-10 (WID) Written to replace fix_underscore.
% * 2013-09-10 (WID) Minor changes to documentation and added support for
%                    cell string inputs.
% * 2015-01-02 (WID) Facilitate cell arrays.  Return type now the same as
%                    the input type.

function result = texstr(input, escChars)

% Generate a list of escape characters by using the provided characters or
% using default values.  Force the forward slash '\' character to be first
% in the list to prevent escape characters from being escaped.
try
    chars = ['\' escChars];
catch
    % When "escChars" is not provided, default to '\' and '_'.
    chars = ['\', '_'];
end

% Convert strings into cell array so that the following algorithm will work
if ~iscell(input)
    input_cell = {input};
else
    input_cell = input;
end

% For each string in the cell array, prefix the required characters.
for r = 1:size(input_cell,1)
    for c = 1:size(input_cell,2)
        for ch = chars
            input_cell{r,c} = strrep(input_cell{r,c}, ch, ['\' ch]);
        end
    end
end

% Convert output to the same format as the input. ie, cell or string
if ~iscell(input)
    result = input_cell{1};
else
    result = input_cell;
end

