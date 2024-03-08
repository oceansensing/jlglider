using GLMakie
using WAV
using Statistics
using DSP

file_path = "/Users/jack/oceansensing Dropbox/C2PO/glider/gliderData/electa-20231127-norse-complete/loggerhead/20231116T063756_9040035954096020_2.0dB_3.5V_ver3.00.wav";
signal_path = "/Users/jack/Documents/acoustic_code/julia/up4.1s_4.2s_610_890_48k.wav"

wav_data, Fs = wavread(file_path);
signal_data, Fs_sig = wavread(signal_path);

duration = (length(wav_data)) / Fs;
t = LinRange(0, duration, length(wav_data));

signal_data = signal_data[:,1]/maximum(signal_data[:,1]);
normmf = sum(abs.(signal_data.^2));

wav_mean_sub = wav_data[:,1] .- mean(wav_data[:,1]);

corrout0 = xcorr(wav_mean_sub, signal_data; padmode = :longest)./normmf;
corrout = corrout0[floor(Int, (length(corrout0) + 1) / 2):end];
mf_out = abs.(hilbert(corrout));

duration = (length(mf_out)) / Fs;
t = LinRange(0, duration, length(mf_out));

fig = Figure(size=(1000,800), fontsize = 20);
ax = Axis(fig[1,1],
    title = "Matched Filter Processing",
    xlabel = "Time (s)",
    ylabel = "Similarity"
);

GLMakie.lines!(t,mf_out);
fig