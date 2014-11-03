import numpy as np
from scipy.signal import get_window
from scipy.fftpack import fft
import math
import matplotlib.pyplot as plt
eps = np.finfo(float).eps

""" 
A4-Part-1: Extracting the main lobe of the spectrum of a window

Write a function that extracts the main lobe of the magnitude spectrum of a window given a window type 
and its length (M). The function should return the samples corresponding to the main lobe in decibels (dB).

To compute the spectrum, take the FFT size (N) to be 8 times the window length (N = 8*M) 
(For this part, N need not be a power of 2). 

The input arguments to the function are the window type (window) and the length of the window (M). 
The function should return a numpy array containing the samples corresponding to the main lobe of 
the window. 

In the returned numpy array you should include the samples corresponding to both the local minimas
across the main lobe. 

The possible window types that you can expect as input are rectangular ('boxcar'), 'hamming' or
'blackmanharris'.

EXAMPLE: If you run your code using window = 'blackmanharris' and M = 100, the output numpy array 
should contain 65 samples.

NOTE: You can approach this question in two ways: 1) You can write code to find the indices of the 
local minimas across the main lobe. 2) You can manually note down the indices of these local minimas 
by plotting and a visual inspection of the spectrum of the window. If done manually, the indices have 
to be obtained for each possible window types separately (as they differ across different window types).

Tip: log(0) is not well defined, so its a common practice to add a small value such as eps = 1e-16 to the 
magnitude spectrum before computing it in dB. This is optional and will not affect your answers.
"""
def extractMainLobe(window, M):
    """
    Input:
            window (string): Window type to be used (Either rectangular ('boxcar'), 'hamming' or '
                blackmanharris')
            M (integer): length of the window to be used
    Output:
            The function should return a numpy array containing the main lobe of the magnitude spectrum 
            of the window in decibels (dB).
    """

    ### Your code here
    print window + " len = " + str(M)
    w = get_window(window, M)
    N = 8*M
    hN = N/2
    hM = M/2
    hM1 = int(math.floor((w.size+1)/2))                     # half analysis window size by rounding
    hM2 = int(math.floor(w.size/2))                         # half analysis window size by floor
    # zero pad on both sides of the window
    fftBuffer = np.zeros(N)
    mX1 = np.zeros(N)
    fftBuffer[hN-hM1:hN+hM2] = w
    # Find DFT of window
    W = fft(fftBuffer)   
    absX = abs(W)
    absX[absX < eps] = eps    # if zeros add epsilon to handle log
    mX = 20*np.log10(absX)
    print mX
    # Find the local minimum from the center
    Min = 0
    for i in np.arange(0,hN):
        if(mX[i] - mX[i+1] < 0.0):
            print "Found Minimum at: " + str(i) + " Total Len = " + str(2*i + 1)
            break        
    
    op = np.zeros(i*2+1)
    op[i:len(op)] = mX[0:i+1]
    op[0:i] = mX[N-i:N]
    print op
    plt.plot(mX)
    mX1[:hN] = mX[hN:]
    mX1[N-hN:] = mX[:hN]  
    # Here we normalize the window so that the peak has value of 1dB
    # It's also OK to plot just mX1, in which case you need to extend the Y axes to say 50dB to -120 dB
    #plt.plot(np.arange(-hN, hN), mX1-max(mX), 'r', lw=1.5) 
    #plt.axis([-hN,hN,-40,0])
    #plt.title('mW, N = 1024')
    
    return(op)
    
    
    
    
    
#You can put the code that calls the above functions down here
if __name__ == "__main__":
    extractMainLobe('blackmanharris',100)
    
    

            