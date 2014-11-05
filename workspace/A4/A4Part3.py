import os
import sys
import numpy as np
from scipy.signal import get_window
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../software/models/'))
import stft
import utilFunctions as UF

eps = np.finfo(float).eps

"""
A4-Part-3: Computing band-wise energy envelopes of a signal

Write a function that computes band-wise energy envelopes of a given audio signal by using the STFT.
Consider two frequency bands for this question, low and high. The low frequency band is the set of all the 
frequencies from 0 - 3000 Hz and the high frequency band is the set of all the frequencies from 
3000 - 10000 Hz. At a given frame, the value of the energy envelope of a band can be computed as the 
sum of squared values of all the frequency coefficients in that band. Compute the energy envelopes in 
decibels. 

Refer to "A4-STFT.pdf" document for further details on computing bandwise energy.

The input arguments to the function are the wav file name including the path (inputFile), window 
type (window), window length (M), FFT size (N) and hop size (H). The function should return a numpy 
array with two columns, where the first column is the energy envelope of the low frequency band and 
the second column is that of the high frequency band.

Use stft.stftAnal() to obtain the STFT magnitude spectrum for all the audio frames. Then compute two 
energy values for each frequency band specified. While calculating frequency bins for each frequency 
band, consider only the bins that are within the specified frequency range. For example, for the low 
frequency band consider only the bins with frequency > 0 Hz and < 3000 Hz. This way we also remove the 
DC offset in the signal in energy envelope computation.

To get a better understanding of the energy envelope and its characteristics you can plot the envelopes 
together with the spectrogram of the signal. You can use matplotlib plotting library for this purpose. 
To visualize the spectrogram of a signal, a good option is to use colormesh. You can reuse the code in
sms-tools/lectures/4-STFT/plots-code/spectrogram.py. Either overlay the envelopes on the spectrogram 
or plot them in a different subplot. Make sure you use the same range of the x-axis for both the 
spectrogram and the energy envelopes.

EXAMPLE: Running your code on piano.wav file with window = 'blackman', M = 513, N = 1024, H = 128, in 
the plots you can clearly notice the sharp attacks and decay of the piano notes (See figure in the 
accompanying pdf). In addition, you can also visually analyse which of the two energy envelopes is better 
for detecting onsets of the piano notes.
"""
def computeEngEnv(inputFile, window, M, N, H):
    """
    Inputs:
            inputFile (string): input sound file (monophonic with sampling rate of 44100)
            window (string): analysis window type (choice of rectangular, triangular, hanning, hamming, 
                blackman, blackmanharris)
            M (integer): analysis window size (odd positive integer)
            N (integer): FFT size (power of 2, such that N > M)
            H (integer): hop size for the stft computation
    Output:
            The function should return a numpy array engEnv with shape Kx2, K = Number of frames
            containing energy envelop of the signal in decibles (dB) scale
            engEnv[:,0]: Energy envelope in band 0 < f < 3000 Hz (in dB)
            engEnv[:,1]: Energy envelope in band 3000 < f < 10000 Hz (in dB)
    """
    
    ### your code here
    F1 = 3000
    F2 = 10000

    (fs, x) = UF.wavread(inputFile)
    w = get_window(window, M)
    mX, pX = stft.stftAnal(x, fs, w, N, H)
    y = stft.stftSynth(mX, pX, M, H)
    print "length x = " + str(x.size) + " Length y = " + str(y.size)
    mX = np.power(mX,10)/20 # de-convert from decibel
    # Note that mX is an array of DFT frames corresponding to each frame in the entire signal.
    # Each ROW of mX is a DFT of the particular frame. Toral numebr of rows is the total number of frames.
    
    # Find number of bins in the region 0 - 3k, f = k*Fs/N
    k0 = 0
    k1 = F1*N/fs # Bin number/index of first limit
    k2 = F2*N/fs
    R1 = [k0,k1]
    R2 = [k1,k2]
    numFrames = int(mX[:,0].size)
    print "Number of Frames = " + str(numFrames)
    # Find the energy in the spectrum f > 0 and <3k
    e1 = np.zeros((numFrames,k1-k0-1)) 
    e2 = np.zeros((numFrames,k2-k1-1))
    for i in range(0,numFrames):
        e1[i] = np.sum( mX[i,k0+1:k1] * mX[i,k0+1:k1])
        e2[i] = np.sum( mX[i,k1+1:k2] * mX[i,k1+1:k2])
    
    # Convert back to dB
    e1 = 10*np.log10(e1) 
    e2 = 10*np.log10(e2)
    
    

    
    
    
    
    
#You can put the code that calls the above functions down here
if __name__ == "__main__":
    computeEngEnv('../../sounds/piano.wav','blackman',512,1024,128)