%% make_scientific
% Convert number into a string in scientific notation
%%
% <latex>\index{Functions!make\_scientific}</latex>
%
%%% Syntax
%   scString = make_scientific( N, Digits, force )
%
% * [N] Number to convert to scientific notation
% * [digits] Number of significant digits required in the output
% * [force] Force inclusion of exponent even when equal to 0
% * []
% * [scString] String representation of 'N' in scientific notation.  When 
%         'N' is a vector of length greater than 1, 'scString' is a cell 
%         array of strings.
%
%%% Description
% Converts the number $\texttt{N}$ into a string representation in
% scientific notation. The resulting string is optimized for use in Matlab
% plotting. When $\texttt{N}$ is a matrix, each value is converted, and
% returned in a cell array of strings.
%
%%% Examples
%
%   >> plotLabel = make_scientific( pi/100, 5 )
%
% Returns Pi/100 with 5 significant digits:
%
%        3.1416 \times 10^{-2}
%
% Rendered in a MATLAB plot as:
%
% <latex>\hspace{0.5in} $3.1416 \times 10^{-2}$</latex>
%

% Version History
%
% * 2012-09-03 (WID) documentation changed for Matlab publishing
% * 2012-09-03 (WID) added support for matrix input / string array output
% * 2012-11-05 (WID) documentation update
% * 2016-06-21 (WID) when not forced, remove the exponent when close to 0
% * 2019-01-08 (RGL) when not forced, remove the exponent when exponent == 0.

function scString = make_scientific(N, Digits, force)

if nargin < 3, force = false; end

if ~isnumeric(N) || ~isnumeric(Digits)
    error ('Input must be a number')
end

Digits = round(Digits); % in case of non-integer number
if Digits < 1, Digits = 1; end
scString = {};
sci_format = ['%1.' num2str(Digits-1) 'f'];

for Num = N,
    n = find(Num < 0); % find negative values
    Num(n) = -Num(n); % temporarily make negative values positive
    M = log10(Num); % take base 10 logarithm
    exponent = floor(M);
    mantissa = M - floor(M);
    mantissa = 10^(mantissa);
    mantissa(n) = -mantissa(n); % return negative values to les than zero
    if force || exponent ~= 0 
%        sc_string = [num2str(mantissa, Digits) ' \times 10^{' num2str(exponent) '}'];
        sc_string = [sprintf(sci_format, mantissa) ' \times 10^{' num2str(exponent) '}']; % retain the trailing zeros
    else
        sc_string = sprintf(sci_format, mantissa * 10^exponent); % retain the trailing zeros
%        sc_string = num2str(mantissa * 10^exponent, Digits);
    end
    scString{end + 1} = sc_string;
end

if length(scString) == 1, scString = scString{1}; end

end
