using FFTW, DSP

ii = 1

mrp = norse23mr[ii].mr;
mrpz = norse23mr[ii].z[:];

minind = findall(mrpz .== minimum(mrpz))
tbott = round(mrp.t_fast[minind][1], RoundDown);
tinddn = findall(tbott .>= mrp.t_fast[:] .>= 0);

# Generate a sample signal (replace this with your actual signal)
fs = 512  # Sample rate (Hz)
t = 0:1/fs:tbott  # Time vector
#signal = sin.(2π * 5t) + sin.(2π * 50t)  # A signal with two frequencies: 5Hz and 50Hz

# Design a high-pass filter
cutoff_frequency_hp = 0.1  # Cutoff frequency (Hz)
nyquist_rate = fs / 2
normalized_cutoff_hp = cutoff_frequency_hp / nyquist_rate
digital_filter_hp = digitalfilter(DSP.Highpass(normalized_cutoff_hp), DSP.Butterworth(1))

# Design a low-pass filter
cutoff_frequency_lp = 0.5 # Cuttoff frequency (Hz)
normalized_cutoff_lp = cutoff_frequency_lp / nyquist_rate
digital_filter_lp = digitalfilter(DSP.Lowpass(normalized_cutoff_lp), DSP.Butterworth(4))


# Apply the highpass filter
sh1_hp_raw = DSP.filtfilt(digital_filter_hp, mrp.sh1[tinddn]);

# rectifying the signal
sh1_hp = abs.(sh1_hp_raw);

# create a copy of the rectified signal
sh1_hp_cp = deepcopy(sh1_hp);

# Applying the lowpass filter
sh1_lp = DSP.filtfilt(digital_filter_lp, sh1_hp_cp);

spikethreshold = 13.0;
bind = findall(sh1_hp ./ sh1_lp .> spikethreshold);
gind = findall(sh1_hp ./ sh1_lp .<= spikethreshold);

Plots.plot(sh1_hp[gind])
