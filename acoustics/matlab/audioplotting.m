[electa_wav, Fs] = audioread('/Users/jack/oceansensing Dropbox/C2PO/glider/gliderData/electa-20231127-norse-complete/loggerhead/20231116T061756_9040035954096020_2.0dB_3.5V_ver3.00.wav');
[signal, Fs_sig] = audioread('up4.1s_4.2s_610_890_48k.wav');

signal = signal/max(signal);
normmf=sum(abs(signal.^2));

corrout0=xcorr(electa_wav(:,1)-mean(electa_wav(:,1)),signal)./normmf;
corrout= corrout0(fix((length(corrout0)+1)/2) : end );
mf_out1 = abs(hilbert(corrout));
 
%mf_out=mf_out1(1:2:end);   % half the sample rate...
%fsmf=Fs./2;

duration = length(mf_out1) / Fs;
t = linspace(0, duration, length(mf_out1));

plot(t, mf_out1)

% Find local maxima
[peaks, locs] = findpeaks(mf_out1);

% Define the range for peak values
min_peak_value = 0.5e-4;  % Minimum peak value
max_peak_value = 2e-3;  % Maximum peak value

% Filter indices based on peak values within the range
valid_indices = locs(peaks >= min_peak_value & peaks <= max_peak_value);

detection_times = t(valid_indices);
detections = mf_out1(valid_indices);

%scatter(detection_times, detections)