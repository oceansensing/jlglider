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
cutoff_frequency = 0.1  # Cutoff frequency (Hz)
nyquist_rate = fs / 2
normalized_cutoff = cutoff_frequency / nyquist_rate
hp_filter = DSP.Butterworth(4)  # Order 4 filter
digital_filter = digitalfilter(DSP.Highpass(normalized_cutoff), hp_filter)

# Apply the filter
hp_sh1 = DSP.filtfilt(digital_filter, mrp.sh1[tinddn]);
