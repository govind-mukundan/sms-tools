from scipy.fftpack import fft
import numpy as np
from fractions import gcd
import matplotlib.pyplot as plt

"""
A3-Part-1: Minimize energy spread in DFT of sinusoids
Given a signal consisting of two sinusoids, write a function that selects the first M samples from 
the signal and returns the positive half of the DFT magnitude spectrum (in dB), such that it has 
only two non-zero values. 

M is to be calculated as the smallest positive integer for which the positive half of the DFT magnitude 
spectrum has only two non-zero values. To get the positive half of the spectrum, first compute the 
M point DFT of the input signal (for this you can use the fft function of scipy.fftpack, which is 
already imported in this script). Consider only the first (M/2)+1 samples of the DFT and compute the
magnitude spectrum of the positive half (in dB) as mX = 20*log10(abs(X[:M/2+1])), where X is the DFT 
of the input.

The input arguments to this function are the input signal x (of length W >= M) consisting of two 
sinusoids of frequency f1 and f2, the sampling frequency fs and the value of frequencies f1 and f2. 
The function should return the positive half of the magnitude spectrum mX. For this question, 
you can assume the input frequencies f1 and f2 to be positive integers and factors of fs, and 
that M is even. 

Due to the precision of the FFT computation, the zero values of the DFT are not zero but very small
values < 1e-12 (or -240 dB) in magnitude. For practical purposes, all values with absolute value less 
than 1e-6 (or -120 dB) can be considered to be zero. 

EXAMPLE: For an input x sampled at 10000 Hz and composed of sinusoids of frequency 80 Hz and 200 Hz, 
the DFT size for the required condition is M = 250 samples and the non-zero values in the DFT spectrum 
are at bin indices 2 and 5 (corresponding to the frequency values of 80 and 200 Hz, respectively). 
The output mX that the function returns is 126 samples in length. 

"""
def minimizeEnergySpreadDFT(x, fs, f1, f2):
    """
    Inputs:
        x (numpy array) = input signal 
        fs (float) = sampling frequency in Hz
        f1 (float) = frequency of the first sinusoid component in Hz
        f2 (float) = frequency of the second sinusoid component in Hz
    Output:
        The function should return 
        mX (numpy array) = The positive half of the DFT spectrum of the M sample segment of x. 
                           mX is (M/2)+1 samples long (M is to be computed)
    """
    ## Your code here
    # If the DFT is to be an impulse, the bin frequency must be aligned with the sinusoid frequency
    # Bin frequency = multiples of Fs/N; f1 and f2 must be a multiple of the bin frequency
    # So lets set the bin frequency to be the GCD of f1 and f2, this will give you the smallest N
    # that ensures that f1 and f2 are multiples of the bin frequency 
    bin_freq = gcd(f1,f2)
    print 'GCD of ' + str(f1) + ' and ' + str(f2) + ' is ' + str(bin_freq)
    M = fs/bin_freq    
    print 'M = ' +str(M)
    
    # Now lets compute the DFT using FFT
    x1 = x[0:M]
    X = fft(x1)
    # Now return the +ve half of the magnitude spectrum in decibels    
    mX = 20*np.log10(abs(X[:M/2+1]))
    
    return mX
    

def gcd(a,b):
    """Calculate the greatest common divisor of a and b"""
    while b:
        a, b = b, a%b
        
    return a
    
    
    
##########################
#You can put the code that calls the above functions down here    
if __name__ == "__main__":
    f1 = 80
    f2 = 200
    N = 20000
    Fs = 10000
    n = np.arange(0,N)
    x = np.cos(2*np.pi*f1*n/Fs) + np.cos(2*np.pi*f2*n/Fs) 
    
    #x = np.cos(2*np.pi*f1*n/Fs + f1) + np.cos(2*np.pi*f2*n/Fs + f2)+ np.cos(2*np.pi*300*n/Fs + 300) + np.cos(2*np.pi*400*n/Fs + 400)
    #plt.plot(x)
    mX = minimizeEnergySpreadDFT(x, Fs, f1, f2)
    plt.plot(mX,'o')
    print mX