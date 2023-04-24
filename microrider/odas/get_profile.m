%% get_profile
% Extract indices to where an instrument is moving up or down
%%
% <latex>\index{Functions!get\_profile}</latex>
%
%%% Syntax
%   profile = get_profile( P, W, pMin, wMin, direction, minDuration, fs )
%
% * [P] Pressure vector
% * [W] Vector of the rate-of-change of pressure.
% * [pMin] Minimum pressure, for a valid profile.
% * [wMin] Minimum magnitude of the rate-of-change of pressure, for a valid profile.
% * [direction] Direction of profile, either 'up' or 'down'.
% * [minDuration] Minimum duration, for a valid profile [s].
% * [fs] Sampling rate of P and W in samples per second.
% * []
% * [profile] 2 X N matrix where each column contains the start- and
%      end-indices of a profile, according to the input definition. Empty
%      if no profiles were detected.
%
%%% Description
% Extract the sample indices to the sections in a data file where the 
% instrument is steadily moving in $\texttt{direction}$ (up or down) for at
% least $\texttt{minDuration}$ seconds, faster than $\texttt{wMin}$ and at
% a pressure greater than $\texttt{pMin}$. 
%
% Call this function separately for ascents and descents.
%
% Can be used with data files collected with any vertical profiler or
% glider. 

% *Version History:*
%
% * 2011-12-26 (RGL) initial
% * 2012-04-30 (WID) comments added to allow for Matlab publishing
% * 2015-07-27 (WID) Documentation update
% * 2015-10-29 (RGL) Documentation update
% * 2015-11-18 RGL Documentation updates.

function profile = get_profile(P, W, P_min, W_min, direction, min_duration, fs)

min_samples = min_duration*fs; %The minimum number of contiguous samples to be declared a profile.
if strcmpi(direction,'up'); W = -W; end

n = find((P > P_min) & (W >= W_min));
if (length(n) < min_samples)
    profile = [];
    return
end
    
diff_n = diff(n);
diff_n = [diff_n(1); diff_n]; % use diff function to find gaps in profiles

m = find(diff_n > 1);% This locates the breaks bewtween profiles
m = m(m ~= 1); % Reject m == 1

if isempty(m)
    profile = [n(1) ; n(end)]; % There is only a single profile
else
    profile = zeros(2,length(m)+1);
    profile(:,1) = [n(1); n(m(1)-1)]; % This is the first profile
    
    for index = 2:length(m)
        profile(:,index) = [n(m(index-1));n(m(index)-1)];
    end
    profile(:,end) = [n(m(end));n(end)]; % This is the last profile
end
%Now check the length of each profile
profile_length = profile(2,:) - profile(1,:); % the length of each profile

mm = find(profile_length >= min_samples); % find the profiles longer than min_samples
profile = profile(:,mm); % keep only those profiles
