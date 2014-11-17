import numpy as np
from scipy.signal import get_window
import math
import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../software/models/'))
import dftModel as DFT
import utilFunctions as UF

""" 
A5-Part-1: Minimizing the frequency estimation error of a sinusoid

Write a function that estimates the frequency of a sinusoidal signal at a given time instant. The 
function should return the estimated frequency in Hz, together with the window size and the FFT 
size used in the analysis.  

The input arguments to the function are the wav file name including the path (inputFile) containing 
the sinusoidal signal, and the frequency of the sinusoid in Hz (f). The frequency of the input sinusoid  
can range between 100Hz and 2000Hz. The function should return a three element tuple of the estimated 
frequency of the sinusoid (fEst), the window size (M) and the FFT size (N) used.

The input wav file is a stationary audio signal consisting of a single sinusoid of length 1 second. 
Since the signal is stationary you can just perform the analysis in a single frame, for example in 
the middle of the sound file (time equal to .5 seconds). The analysis process would be to first select 
a fragment of the signal equal to the window size, M, centered at .5 seconds, then compute the DFT 
using the dftAnal function, and finally use the peakDetection and peakInterp functions to obtain the 
frequency value of the sinusoid.

Use a Blackman window for analysis and a magnitude threshold t = -40 dB for peak picking. The window
size and FFT size should be chosen such that the difference between the true frequency (f) and the 
estimated frequency (fEst) is less than 0.05 Hz for the entire allowed frequency range of the input 
sinusoid. The window size should be the minimum positive integer of the form 100*k + 1 (where k is a 
positive integer) for which the frequency estimation error is < 0.05 Hz. For a window size M, take the
FFT size (N) to be the smallest power of 2 larger than M. 

HINT: If the specified frequency range would have been 440-8000 Hz, the parameter values that satisfy 
the required conditions would be M = 1101, N = 2048. Note that for a different frequency range, like 
the one specified in the question, this value of M and N might not work. 

"""
def minFreqEstErr(inputFile, f):
    """
    Inputs:
            inputFile (string) = wav file including the path
            f (float) = frequency of the sinusoid present in the input audio signal (Hz)
    Output:
            fEst (float) = Estimated frequency of the sinusoid (Hz)
            M (int) = Window size
            N (int) = FFT size
    """
    # analysis parameters:
    window = 'blackman'
    t = -40
    
    ### Your code here
    (fs, ip) = UF.wavread(inputFile)

    # M >= Bs * Fs / Delta_F where Bs is the width of the main lobe in bins
    Bs = 6 #for blackmann
    
    W = 6 * fs/0.05
    k = 21 # Found this by running hte main code
    M = 100 * k + 1
    N = findNFromM(100*k+1)
    x = ip[fs/2 - M/2 : fs/2 + M/2 + 1]
    w = get_window(window, M)
    mX, pX = DFT.dftAnal(x, w, N)
    ploc = UF.peakDetection(mX, t)
    iploc, ipmag, ipphase = UF.peakInterp(mX, pX, ploc) 
    ipfreq = fs*iploc/float(N)                            # convert peak locations to Hertz
    print "FEst = " + str(ipfreq) + "M = " + str(M) + "N = " + str(N)
    
    return (float(ipfreq),int(M),int(N))

def findNFromM(M):
    
    for i in np.arange(0,M/2):
        if(pow(2,i) > M):
            break
        
    return pow(2,i)
    
def findErrForK(ip, k,f):
    """
    Inputs:
            inputFile (string) = wav file including the path
            f (float) = frequency of the sinusoid present in the input audio signal (Hz)
    Output:
            fEst (float) = Estimated frequency of the sinusoid (Hz)
            M (int) = Window size
            N (int) = FFT size
    """
    # analysis parameters:
    window = 'blackman'
    t = -40
    
    ### Your code here
    
    fs = 44100
    
    # M >= Bs * Fs / Delta_F where Bs is the width of the main lobe in bins
    #Bs = 6 #for blackmann
    
    #W = 6 * fs/0.05
    N = findNFromM(100*k+1)
    t = -40

    M = 100 * k + 1
    x = ip[fs/2 - M/2 : fs/2 + M/2 + 1]
    w = get_window(window, M)
    mX, pX = DFT.dftAnal(x, w, N)
    ploc = UF.peakDetection(mX, t)
    iploc, ipmag, ipphase = UF.peakInterp(mX, pX, ploc) 
    ipfreq = fs*iploc/float(N)                            # convert peak locations to Hertz
        
    return f - ipfreq
    
#You can put the code that calls the above functions down here
if __name__ == "__main__":
    #minFreqEstErr('../../sounds/sine-440.wav',150)
    minFreqEstErr('../../sounds/sine-440.wav',150)
    Fs = 44100
    n = np.arange(0,Fs)
    k_old = 0
    # For each value of K find the error over the entire spectrum.
    # Need to find the smallest K for which the error over the entire specrum is < 0.05Hz
    found = False
    for k in np.arange(1,30):
        found = True
        for i in np.arange(100,2000):
            x = np.cos(2*np.pi*i*n/Fs)
            err = findErrForK(x,k,i)
            if(err > 0.05):
                print "k not = " +str(k)
                found = False
                break
        if(found == True):
            print "k = " +str(k)
            break
    
    print "M = " + str(100*k + 1)