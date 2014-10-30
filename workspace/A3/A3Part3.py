import numpy as np
from scipy.fftpack import fft, fftshift
import math

"""
A3-Part-3: Symmetry properties of the DFT

Write a function to check if the input signal is real and even using the symmetry properties of its DFT. 
The function will return the result of this test, the zerophase windowed version of the input signal 
(dftbuffer), and the DFT of the dftbuffer. 

Given an input signal x of length M, do a zero phase windowing of x without any zero-padding (a dftbuffer, 
on the same lines as the fftbuffer in sms-tools). Then compute the M point DFT of the zero phase windowed 
signal and use the symmetry of the computed DFT to test if the input signal x is real and even. Return 
the result of the test, the dftbuffer computed, and the DFT of the dftbuffer. 

The input argument is a signal x of length M. The output is a tuple with three elements 
(isRealEven, dftbuffer, X), where 'isRealEven' is a boolean variable which is True if x is real 
and even, else False. dftbuffer is the M length zero phase windowed version of x. X is the M point 
DFT of the dftbuffer. 

To make the problem easier, we will use odd length input sequence in this question (M is odd). 

Due to the precision of the FFT computation, the zero values of the DFT are not zero but very small
values < 1e-12 (or -240 dB) in magnitude. For practical purposes, all values with absolute value less 
than 1e-6 (or -120 dB) can be considered to be zero. Use an error tolerance of 1e-6 on the linear 
scale to compare if two floating point arrays are equal. 

EXAMPLE: If x = np.array([ 2, 3, 4, 3, 2 ]), which is a real and even signal, the function returns 
(True, array([ 4., 3., 2., 2., 3.]), array([14.0000+0.j, 2.6180+0.j, 0.3820+0.j, 0.3820+0.j, 2.6180+0.j]))
(values are approximate)

"""
def testRealEven(x):
    """
    Inputs:
        x (numpy array)= input signal of length M (M is odd)
    Output:
        The function should return a tuple (isRealEven, dftbuffer, X)
        isRealEven (boolean) = True if the input x is real and even, and False otherwise
        dftbuffer (numpy array, possibly complex) = The M point zero phase windowed version of x 
        X (numpy array, possibly complex) = The M point DFT of dftbuffer 
    """
    ## Your code here
    # IF both real ane even, mag(F(k)) = mag(F(N-k)) AND phase(F(k)) = 0
    # IF only real, mag(F(k)) = mag(F(N-k)) AND phase(F(k)) = -phase(F(N-k))
    
    # 1. Zero phase window --> cyclic shift to align the midpoint of the signal to the start fo the buffer
    #x = x * np.hamming(len(x))    
    x = np.roll(x, (x.size+1)/2)
    #print x
    # 2. compute DFT
    X = fft(x)
    M = len(x)
    mX = 20*np.log10(abs(X))

    # Beaautify the output to give zeros for values less than -120dB        
    for i in range(0,len(mX)):
        #print mX[i]
        if(mX[i] < -120):
            mX[i] = 0
    mag_flag = False
    ph_flag = False
    print(mX)
    # 3. Test for real + even
    for i in range(1,M/2+1):
        print 'M1 = ' + str(abs(mX[i])) + ' M2 = ' + str( abs(mX[M-i]))
        if (abs(mX[i]) != abs(mX[M-i])):
            mag_flag = True;
            
    # Note that phase of Pi is also same as 0
    for i in range(0,M):
        print 'PH = ' + str(np.angle(X[i])) 
        if ( (abs(np.angle(X[i])) < 1e-6) or (abs(np.angle(X[i])) > 3.140 and (abs(np.angle(X[i])) < 3.143)) ):
            print ""
        else:
            ph_flag = True
            print "Phase condition failed!"
    
    if((ph_flag == False) and (mag_flag == False)):
        print "Real Even Condition True"
        RealEven = True
    else:
        print "Real Even Condition False"
        RealEven = False
        
        
    return RealEven,x,X
        
        
##########################
#You can put the code that calls the above functions down here
if __name__ == "__main__":
    
    x = np.array([ 2, 3, 4, 3, 2 ]);
    testRealEven(x)
    
        
    