%% median_filter
% Used for removing flyers from ADV data
%%
% <latex>\index{Functions!median\_filter}</latex>
%
%%% Syntax
%   [result, points] = median_filter( dIn, threshold, filterLen, 
%                                     extraPoints, stDev, k )
%
% * [dIn] Matrix of column vectors to clean. Typically three columns -- one
%       for each velocity component.
% * [threshold] Data points are defined as bad when they differ by more than
%       this value from the median. Can be empty if using the stDev option.
% * [filterLen] Length of data segment [samples] used to calculate the median
%       and standard deviation. 
% * [extraPoints] Number of points, adjacent to bad points, to also replace
%      with interpolated values.
% * [stDev] empty string '', or 'std_dev'.
% * [k] empty, or the scale factor for the stDev option.
% * []
% * [result] Copy of dIn with bad data points replaced with interpolated values.
% * [points] Vector of indicies to the bad data samples. Length of points
%       equals the number of bad samples detected.
%
%
%%% Description
%
% Use a median filter to detect and remove flyers or outliers - points that
% are obviously wrong because they differ greatly from a local median.
% 
% This function is directed at the cleaning of ADV (acoustic doppler
% velocimeter) data using only velocity data without any information about
% correlation coefficients or other metrics. Can be used with other data. 
%
% The function detects erroneous data by calculating the median and standard
% deviation of segments of length $\texttt{filterLen}$. Data points that
% deviate by more then $\texttt{threshold}$ from the median are flagged as
% bad. The value of $\texttt{threshold}$ can be specified explictly, in
% which case the values of $\texttt{stDev}$ and $\texttt{k}$ are ignored.
% If $\texttt{threshold}$ is empty, and $\texttt{stDev = 'st\_dev'}$, the
% threshold is set to $\texttt{k}$ times the standard deviation. That is,
% the threshold may change from segement to segment. The standard deviation
% of a segment may be highly elevated due to a flyer, and it may be prudent
% to clean the data using two passes of this function. First with a
% very high value of $\texttt{k}$, say $\num{\sim8}$, to remove extreme
% flyers, and then with a lower value of $\num{\sim4}$.
%
% When dIn has multiple columns, a bad data point in one column invalidates
% the corresponding data points in all other columns. This behaviour is
% suited to cleaning ADV data, for which an error in one velocity component
% implies an error in the other components. If this behaviour is
% undesirable, use the function with only a single column vector.
%
% The input data may be from the analog output voltage of an ADV, in which
% case it is frequently oversampled. A short neighbourhood adjacent to the
% bad points may then also be contaminated by a flyer. Use the
% $\texttt{extraPoints}$ option to remove such points. A value of 3 seems
% to work well for analog ADV data sampled at 512 Hz, when the ADV itself
% was updating its output at a rate of $\SI{4}{\per\s}$.

% Version History:
%
% 2014-06-07, RGL and WID, Original version.
% 2014-12-05, RGL added option to use the standard deviation to determine a
%           threshold. A threshold of 3.5 to 4.5 times the standard
%           deviations appears to work well.
% 2015-03-31, RGL, forced extra_points to extra_points+1, so that the loops
%           work properly, when extra_points=0. The functional result is
%           unchanged.
% 2015-04-15, WID, Updated documentation for publishing.  Modified some 
%           variable names to avoid naming conflicts.
% 2015-08-25, WID, Modified algorithm to improve speed when excessive bad
%           points are detected.
% 2015-08-29, RGL, Change documentation.
% 2015-11-18 RGL, Updated documentation.
% 2015-12-23 RGL, Added check for 5 input arguments only.
% 2017-01-26 RGL, Replaced interpolation with linspace function. I was
%           getting size mis-match errors. 

function [result, points] = median_filter( dIn, threshold, filterLen, ...
                                           extraPoints, stDev, k )

if nargin < 4,
    extraPoints = 1;
    stDev = false;
    k = [];
elseif nargin < 5
    stDev = false;
    k = [];
elseif nargin == 6 && strcmpi(stDev,'st_dev') && isscalar(k)
    stDev = true;
else
    error('Input arguments 5 and 6 are not correct')
end

extraPoints = extraPoints + 1; %else loops will fail
result = dIn;
bad_points = zeros(length(dIn),1);

for col = dIn
    bad_points = bad_points + find_bad_points( col, threshold, ...
        filterLen, stDev, k);
end

points = find(bad_points);

last = 0;
while true
    p = find(bad_points(last+1:end),1) + last;

    % No more bad points - exit.
    if isempty(p) ||  p == last, break; end
    
    % Find the start point.
    start = 1;
    for i = p-extraPoints:-extraPoints:1
        if ~bad_points(i), start = i; break; end
    end
    start_value = dIn(start,:);
    
    % Find the end point.
    last = length(dIn);
    for i = p+extraPoints:extraPoints:length(dIn)
        if ~bad_points(i), last = i; break; end
    end
    last_value = dIn(last,:);
    
    % If the bad point is at the start or end... fudge values.
    if p == 1
        start_value = last_value;
    elseif p == length(dIn)
        last_value = start_value;
    end
    
    % Fill the gap
    d = ones(last-start+1,size(result,2));
%    size(d)
    for i = 1:size(result,2)
        if start_value(i) == last_value(i)
            d(:,i) = ones(last-start+1,1) * start_value(i);
        else
%            [last start]
%            [start_value(i) last_value(i)]
            d(:,i) = linspace(start_value(i),last_value(i), size(d,1))';
%            d(:,i) = (start_value(i):((last_value(i) - start_value(i)) / (last - start)):last_value(i))';
        end
    end
    
    result(start+1:last-1,:) = d(2:end-1,:);
    
    % We reached the end of the vector - exit.
    if length(result) <= last, break; end
end
end



function points = find_bad_points( dIn, threshold, filter_len, st_dev, k)

    points = zeros(size(dIn));

    % Smooth the input vector
    smoothed = fast_median_filter( dIn, filter_len );
    if st_dev,
        std_of_segments = fast_std_filter( dIn, filter_len );
        threshold = k*std_of_segments;
    end

    % Mark points when input > smoothed vector by the threshold value.
    points(abs(dIn - smoothed) > threshold) = 1;
end

function result = moving_average_filter( dIn, len )
% Rolf sees no value in using a moving average. Currently not implemented.
    result = zeros(size(dIn));
    w = zeros(len,1);
    wi = 0;
    s = 0;
    c = 0;

    for i = 1:length(dIn)

        if c == len, s = s - w(wi+1); c = c - 1; end

        s = s + dIn(i);
        c = c + 1;

        w(wi+1) = dIn(i);
        wi = mod(wi+1, len);

        result(i) = s / c;
    end

end

function result = moving_median_filter( dIn, len )
% Rolf sees no value in using a miving median. It increases the computation
% by a factor of filter_len without any significant benefit. Not
% implemented.
    result = zeros(size(dIn));

    for i = 1:length(dIn)
        a = max(1,round(i - len/2));
        b = min(round(i + len/2), length(dIn));
        result(i) = median(dIn(a:b));
    end
end

function result = fast_median_filter( dIn, len )
% This median filter calculates the mediam of sections of length filter_len
% and takes care of a possible runt segmens at the end of input.
    in = dIn(1:end-mod(length(dIn),len));
    in = reshape(in, len, []);

    result = median(in);
    result = repmat(result,len,1);
    result = result(:);

    runt_points = [];
    if length(dIn) ~= length(result)
        m = median(dIn(end-mod(length(dIn),len)+1:end));
        runt_points = ones(mod(length(dIn), len),1) .* m;
    end

    result = [result;runt_points];
end

function result = fast_std_filter( dIn, len)

    in = dIn(1:end-mod(length(dIn),len));
    in = reshape(in, len, []);

    result = std(in);
    result = repmat(result,len,1);
    result = result(:);

    runt_points = [];
    if length(dIn) ~= length(result)
        m = median(dIn(end-mod(length(dIn),len)+1:end));
        runt_points = ones(mod(length(dIn), len),1) .* m;
    end

    result = [result;runt_points];
end
