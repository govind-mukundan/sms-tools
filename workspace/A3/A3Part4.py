import sys
sys.path.append('../../software/models/')
from dftModel import dftAnal, dftSynth
from scipy.signal import get_window
import matplotlib.pyplot as plt
import numpy as np
"""
A3-Part-4: Suppressing frequency components using DFT model

Given a frame of the signal, write a function that uses the dftModel functions to suppress all the frequency 
components till the first bin >= 70Hz in the signal and returns the output of the dftModel with and without filtering. 

You will use the DFT model to implement a very basic form of filtering to suppress frequency components. 
When working close to mains power lines, there is a 50/60 Hz hum that can get introduced into the 
audio signal. You will try to remove that using a basic DFT model based filter. You will work on just 
one frame of a synthetic audio signal to see the effect of filtering. 

You can use the functions dftAnal and dftSynth provided by the dftModel file of sms-tools. Use dftAnal 
to obtain the magnitude spectrum (in dB) and phase spectrum of the audio signal. Set the values of 
the magnitude spectrum that correspond to frequencies <= 70 Hz to -120dB (there may not be a bin 
corresponding exactly to 70 Hz, choose the nearest bin of equal or higher frequency). Use dftSynth 
to synthesize the filtered output signal and return the output. The function should also return the 
output of dftSynth without any filtering (without altering the magnitude spectrum in any way). 
You will use a hamming window to smooth the signal. Hence, do not forget to scale the output signals 
by the sum of the window values (as done in sms-tools/software/models_interface/dftModel_function.py).
To understand the effect of filtering, you can plot both the filtered output and non-filtered output of the dftModel. 

Please note that this question is just for illustrative purposes and filtering is not usually done 
this way - such sharp cutoffs introduce artifacts in the output. 

The input is a M length input signal x that contains undesired frequencies below 70 Hz, sampling 
frequency fs and the FFT size N. The output is a tuple with two elements (y, yfilt), where y is the 
output of dftModel with the unaltered original signal and yfilt is the filtered output of the dftModel.

EXAMPLE: For an input signal with 40 Hz, 100 Hz, 200 Hz, 1000 Hz components, yfilt will only contain
100 Hz, 200 Hz and 1000 Hz components. 

"""
def suppressFreqDFTmodel(x, fs, N):
    """
    Inputs:
        x (numpy array) = input signal of length M (odd)
        fs (float) = sampling frequency (Hz)
        N (positive integer) = FFT size
    Outputs:
        The function should return a tuple (y, yfilt)
        y (numpy array) = Output of the dftSynth() without filtering (M samples long)
        yfilt (numpy array) = Output of the dftSynth() with filtering (M samples long)
    The first few lines of the code have been written for you, do not modify it. 
    """
    print str(fs) + "---" +  str(N) + "===" + str(len(x))
    M = len(x)
    w = get_window('hamming', M)
    outputScaleFactor = sum(w)
    ## Your code here
    # compute the dft of the sound fragment
    mX, pX = dftAnal(x, w, N)
    y = dftSynth(mX, pX, M)*outputScaleFactor
    # Filter - set the bins < 70Hz = 0 or -120dB
    # Note that DFT.Anal() returns the +ve side of the spectrum only
    nbin=np.ceil(N*70.0/fs)
    mX[:nbin+1]=-120
    print "Number of Bins = " + str(nbin+1)
    nbin = 0
    for k in range(0,len(mX)):
        f = np.ceil(k*fs/N)
        if(f <= 80.0):
            print "freq: " + str(f) + " Bin Index = " + str(k)
            nbin = nbin + 1
    print "Number of Bins 2x = " + str(nbin)
    # compute the inverse dft of the spectrum
    yfilt = dftSynth(mX, pX, M)*outputScaleFactor

    return(y,yfilt)


##########################
#You can put the code that calls the above functions down here
if __name__ == "__main__":
    
    f1 = 71
    f2 = 200
    f3 = 300
    f4 = 500
    N = 1024
    Fs = 10000
    n = np.arange(0,N)
    x = np.cos(2*np.pi*f1*n/Fs) + np.cos(2*np.pi*f2*n/Fs) + np.cos(2*np.pi*f3*n/Fs) + np.cos(2*np.pi*f4*n/Fs) 
    y,yfilt = suppressFreqDFTmodel(x,Fs,N)
    plt.subplot(311)
    plt.plot(y)
    plt.subplot(312)
    plt.plot(yfilt)
    
    M = len(x)
    w = get_window('hamming', M)
    outputScaleFactor = sum(w)
    mX, pX = dftAnal(yfilt, w, N)
    plt.subplot(313)
    plt.plot(float(Fs)*np.arange(mX.size)/float(N), mX, 'r')
    plt.axis([0, Fs/2.0, min(mX), max(mX)])
    plt.title ('magnitude spectrum: mX')
    plt.ylabel('amplitude (dB)')
    plt.xlabel('frequency (Hz)')