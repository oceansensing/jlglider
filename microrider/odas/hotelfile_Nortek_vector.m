%% hotelfile_Nortek_vector
% Generate a hotel-file from Nortek Vector dat- and sen-files
%%
% <latex>\index{Functions!hotelfile\_Nortek\_vector}</latex>
%
%%% Syntax
%    dat = hotelfile_Nortek_vector(file_name, sampling_rate )
%
% * [file_name] Name of data file to read.
% * [sampling_rate] The rate of sampling. If absent, the *.hdr file in the
%      local directory is used to determine the sampling rate.
% * []
% * [dat] Optional. When specified, the vector data are return in this
%      structure. Otherwise, the data are written into a hotel file. 
%
%%% Description
% Generate a hotel-file containing the data from a Nortek Vector ADV
% $\texttt{*.dat}$ and $\texttt{*.sen}$ files. If a $\texttt{*.hdr}$ file
% exists, it will be used to determine the sampling rate of the data in the
% $\texttt{*.dat}$ file. 
%
% The hotel file is a mat-file with the same base name as the data files.
% The hotel file can be used with $\texttt{odas\_p2mat}$ and
% $\texttt{quick\_look}$ to incoprporate the ADV data into an RSI data
% mat-file.
%
%%% Notes
%
% # Nortek Vector instruments with outdated firmware can generate data
% files containing $``(null)"$ entries. This function can not process such
% files. Replacing all $``(null)"$ entries with $``0"$ solves this
% problem. Ensure the Nortek Vector firmware is properly updated to prevent
% this problem from occuring.
%

% "*.sen" files contain the following Vector data.
%
%  1   Month                            (1-12)
%  2   Day                              (1-31)
%  3   Year
%  4   Hour                             (0-23)
%  5   Minute                           (0-59)
%  6   Second                           (0-59)
%  7   Error code
%  8   Status code
%  9   Battery voltage                  (V)
% 10   Soundspeed                       (m/s)
% 11   Heading                          (degrees)
% 12   Pitch                            (degrees)
% 13   Roll                             (degrees)
% 14   Temperature                      (degrees C)
% 15   Analog input
% 16   Checksum                         (1=failed)
%
%
% "*.dat":
%  1   Burst counter
%  2   Ensemble counter                 (1-65536)
%  3   Velocity (Beam1|X|East)          (m/s)
%  4   Velocity (Beam2|Y|North)         (m/s)
%  5   Velocity (Beam3|Z|Up)            (m/s)
%  6   Amplitude (Beam1)                (counts)
%  7   Amplitude (Beam2)                (counts)
%  8   Amplitude (Beam3)                (counts)
%  9   SNR (Beam1)                      (dB)
% 10   SNR (Beam2)                      (dB)
% 11   SNR (Beam3)                      (dB)
% 12   Correlation (Beam1)              (%)
% 13   Correlation (Beam2)              (%)
% 14   Correlation (Beam3)              (%)
% 15   Pressure                         (dbar)
% 16   Analog input 1
% 17   Analog input 2
% 18   Checksum                         (1=failed)
%
% 2010-03-01, Rolf Lueck, RSI, origninal version
% 2011-06-24, Rolf Lueck, RSI, changed to using a structure.
% 2014-06-04, WID, Fixed error with certain data files.  Use textscan.
%                  Modified file opening / error reporting.
% 2014-11-10, RGL, reverted back to using importdata function
% 2015-11-13, RGL, Turned into a function consistent with hotel file
%     structures.

function dat = hotelfile_Nortek_vector(file_name, sampling_rate )


% Remove extension if provided.
[P,N,E] = fileparts(file_name);
if isempty(P), P = '.'; end
file_name = [P filesep N];


[P,N,E,fullName] = file_with_ext( file_name, {'.hdr','.HDR'} );
if ~isempty(fullName)
    % The header file is present. Find the Sampling rate from the file.
    header = fileread( fullName );
    exp = 'Sampling rate\W+(\d+.?\d*) Hz';
    [tokens, matches] = regexp( header, exp, 'tokens', 'match' );
    
    if ~isempty(tokens)
        rate = str2double(tokens{1});
    end
end

if nargin >= 2, rate = sampling_rate; end


if ~exist('rate', 'var') || isempty(rate)
    error('Rate required.  Please provide a header file or specify a sampling rate.');
end

[P,N,E,fullName] = file_with_ext( file_name, ...
                                  {'.sen','.SEN'}, ...
                                  'Unable to find .SEN file.  This file is required.');

% Process the sense file only if it was found.
disp(['Loading "' N E '"']);

try
    filedata = importdata(fullName);% This is a built in Matlab function
catch
    error(['Unable to load .sen file. "(null)" elements might be ' ...
           'present.  Replace with "0" and try again.']);
end

fprintf('Data file contains %d vectors or length %d.\n', size(filedata,2), size(filedata,1));
fprintf('Finished loading sensor file "%s". Converting into doubles...\n\n', [N E]);


dat.Year          = filedata(:,3);
dat.Month         = filedata(:,1);
dat.Day           = filedata(:,2);
dat.Hour          = filedata(:,4);
dat.Minute        = filedata(:,5);
dat.Second        = filedata(:,6);

senstime = datenum(dat.Year, dat.Month, dat.Day, ...
                   dat.Hour, dat.Minute, dat.Second);

dat.Error.data    = filedata(:,7);
dat.Error.time    = senstime;
dat.Status.data   = filedata(:,8);
dat.Status.time   = senstime;
dat.BV.data       = filedata(:,9);
dat.BV.time       = senstime;
dat.SoundVel.data = filedata(:,10);
dat.SoundVel.time = senstime;
dat.Heading.data  = filedata(:,11);
dat.Heading.time  = senstime;
dat.Pitch.data    = filedata(:,12);
dat.Pitch.time    = senstime;
dat.Roll.data     = filedata(:,13);
dat.Roll.time     = senstime;
dat.T_vector.data = filedata(:,14);
dat.T_vector.time = senstime;
dat.AnaIn.data    = filedata(:,15);
dat.AnaIn.time    = senstime;
dat.Checksum.data = filedata(:,16);
dat.Checksum.time = senstime;


[P,N,E,fullName] = file_with_ext( file_name, {'.dat','.DAT'} );

if ~isempty(fullName)
    disp(['Loading ""' N E '"']);

    filedata = importdata(fullName);% This is a built in Matlab function
    fprintf('Data file contains %d vectors or length %d.\n', size(filedata,2), size(filedata,1));
    fprintf('Finished loading file "%s". Converting into doubles...\n\n', [N E]);

    dat.Burst    = filedata(:,1);
    dat.Ensemble = filedata(:,2);

    %%%%%
    % Generate a time vector for the contents of the .DAT file.  The .DAT
    % file is broken into "Bursts" that roughly correspond to any time gaps
    % within the .SEN file.  I say roughly because if they do not this
    % algorithm will not work.  I have not observed a case where this
    % assumption does not hold true.
    %%%%%
    
    % The .SEN sample duration in units of MATLAB time (= 1 second)
    t_sen = addtodate(0, 1, 'second');
    
    %%%% First step, generate indices to changes within the senstime vector
    d_sentime = [inf;diff(senstime);inf];
    sen_breaks = find (d_sentime > t_sen);
    
    % The .DAT sample duration in units of MATLAB time.
    t_dat = t_sen / rate;

    % 1 is prepended and postpended to ensure matches to the extents of the
    % bursts.
    burst_breaks = find([1;diff(dat.Burst);1]);
    
    % Generate the time vectors within "burst".
    for i = 1:length(burst_breaks)-1
        % .DAT files start at the same time as .SEN files.  But the time
        % that .SEN files record is the time AFTER capturing a 1 second
        % record of data. So the .DAT file starts (1 - 1/rate) seconds
        % before the first time in .SEN.
        starttime = senstime(sen_breaks(i)) - t_sen + t_dat ;
        burst{i} = ((burst_breaks(i):burst_breaks(i+1)-1) - burst_breaks(i)) * t_dat + starttime;
    end
    
    % Combine time vectors into a single vector.
    dat_time = [];
    for i = 1:length(burst)
        dat_time = [dat_time burst{i}];
    end
    dat_time = dat_time(:);
    
    
    dat.Ux.data       = filedata(:,3);
    dat.Ux.time       = dat_time;
    dat.Uy.data       = filedata(:,4);
    dat.Uy.time       = dat_time;
    dat.Uz.data       = filedata(:,5);
    dat.Uz.time       = dat_time;
    dat.Amp1.data     = filedata(:,6);
    dat.Amp1.time     = dat_time;
    dat.Amp2.data     = filedata(:,7);
    dat.Amp2.time     = dat_time;
    dat.Amp3.data     = filedata(:,8);
    dat.Amp3.time     = dat_time;
    dat.SNR1.data     = filedata(:,9);
    dat.SNR1.time     = dat_time;
    dat.SNR2.data     = filedata(:,10);
    dat.SNR2.time     = dat_time;
    dat.SNR3.data     = filedata(:,11);
    dat.SNR3.time     = dat_time;
    dat.Corr1.data    = filedata(:,12);
    dat.Corr1.time    = dat_time;
    dat.Corr2.data    = filedata(:,13);
    dat.Corr2.time    = dat_time;
    dat.Corr3.data    = filedata(:,14);
    dat.Corr3.time    = dat_time;
    dat.Pres.data     = filedata(:,15);
    dat.Pres.time     = dat_time;
    dat.Ana1.data     = filedata(:,16);
    dat.Ana1.time     = dat_time;
    dat.Ana2.data     = filedata(:,17);
    dat.Ana2.time     = dat_time;
    dat.Checksum.data = filedata(:,18);
    dat.Checksum.time = dat_time;
end


if nargout == 0
    fprintf('Complete, saving result into "%s".\n', [N '.mat']);
    save( [N '.mat'], '-struct', 'dat', '-v6' );
end

if ~exist('dat', 'var')
    fprintf('Unable to find the requested file: %s\n', file_name);
end
