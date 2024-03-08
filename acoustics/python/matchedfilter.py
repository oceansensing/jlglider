import matplotlib.pyplot as plt
from scipy.io import wavfile
import scipy.signal as spsig
import numpy as np

import acoustic_functions as acfunc

electa_wav_path = '/Users/jack/oceansensing Dropbox/C2PO/glider/gliderData/electa-20231127-norse-complete/loggerhead/20231116T063756_9040035954096020_2.0dB_3.5V_ver3.00.wav'
signal_path = 'up4.1s_4.2s_610_890_48k.wav'

Fs, wav = wavfile.read(electa_wav_path)
Fs_sig, signal = wavfile.read(signal_path)

if Fs != Fs_sig:
    raise ValueError("Sample rate of waveform and reference signal do not match")

mfout = acfunc.matched_filter(wav, signal)

duration = int(len(mfout)/Fs)
t = np.linspace(0,duration,len(mfout))
plt.plot(t,mfout)
plt.show()