import matplotlib.pyplot as plt
from scipy.io import wavfile
import scipy.signal as spsig
import numpy as np

def matched_filter(wav, signal):

    signal = signal/np.max(signal)
    normmf = np.sum(np.abs(signal**2))

    corrout = spsig.correlate(wav-np.mean(wav),signal)/normmf

    mfout = np.abs(spsig.hilbert(corrout))

    return mfout