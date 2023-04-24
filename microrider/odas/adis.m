%% adis
% Convert inclinometer data into meaningful raw counts
%%
% <latex>\index{Functions!adis}</latex>
%
%%% Syntax
%   [xOut, oldFlag, errorFlag] = adis( xIn )
%
% * [xIn] Vector of words from inclinometer.
% * []
% * [xOut] Data vector extracted from words as 2s-complimented count values.
% * [oldFlag] Index vector to data elements that are not new.
% * [errorFlag] Index vector to data elements with error-flags set.
%
% Normally, errorFlag and oldFlag will be empty vectors.
%
%%% Description
% The ADIS16209 inclinometer outputs words that combine data with status flags.
% The data and status flags must be seperated before the data can be converted
% into physical units.  This function extracts the data and converts it into
% valid 16bit, 2s-compliment words.
%
% ADIS16209 words have the following format:
%
%       bit15 = new data present
%       bit14 = error flag
%   bit[13:0] = X and Y inclination data, 2s-compliment (X and Y channels)
%   bit[11:0] = temperature data, unsigned              (temperature channel)
%
%%% Examples
%
%    >> [raw, oldData, errors] = adis( ADIS16209_words )
%
% Extract and convert raw count values from ADIS16209 words.
%

% *Version History:*
%
% * 2010-06-07 (RGL) original version
% * 2012-04-30 (WID) added documentation tags for matlab publishing
% * 2012-05-08 (WID) rewrite of documentation + addition of new tags

function [X_Out, old_flag, error_flag] = adis(X_In) 

error_flag = [];
old_flag = find(X_In >= 0); % MS-bit is not set so data is not new.

n = find(X_In < -2^14); % these data points have the MS-bit set and are new.
if ~isempty(n)
    X_In(n) = X_In(n) + 2^15; % shift to clear MS-bit
end

n = find(X_In >= 2^14); % Second MS-bit is set, which indicates an error. 
if ~isempty(n)
    error_flag = n;
    X_In(n) = X_In(n) - 2^14;
end

% The data are 2s-compliment for inclination. So we find the upper half of the 
% 14-bit range. This is the data with negative values. So we shift it down
% to get a proper 2s-compliment conversion. The temperature data is 12-bit
% only and will automatically not get shifted.
n = find(X_In >= 2^13);
if ~isempty(n)
    X_In(n) = X_In(n) - 2^14;
end
X_Out = X_In;

    

