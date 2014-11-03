# -*- coding: utf-8 -*-
"""
Created on Sun Nov  2 15:10:20 2014

@author: govind
"""

import sys
sys.path.append('../software/models/')
from dftModel import dftAnal, dftSynth
from scipy.signal import get_window
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft, ifft


def getDFT(x,N):
    
    X = np.array([])
    for k in range(N):
	s = np.exp(1j*2*np.pi*k/N*np.arange(N))
	X = np.append(X, sum(x*np.conjugate(s)))
 
    return(X)
    
def Test():
    M = 64
    N = 1024
    hN = N/2     
    hM = M/2
    fftbuffer = np.zeros(N)
    mX1 = np.zeros(N)
    
    plt.figure(1, figsize=(9.5, 6))
    fftbuffer[hN-hM:hN+hM]=np.ones(M)
    plt.subplot(2,1,1)
    plt.plot(np.arange(-hN, hN), fftbuffer, 'b', lw=1.5)
    plt.axis([-hN, hN, 0, 1.1])
    plt.title('w (rectangular window), M = 64')
    
    
    X = fft(fftbuffer)
    mX = 20*np.log10(abs(X)) 
    mX1[:hN] = mX[hN:]
    mX1[N-hN:] = mX[:hN]      
    
    plt.subplot(2,1,2)
    plt.plot(np.arange(-hN, hN), mX1-max(mX), 'r', lw=1.5)
    plt.axis([-hN,hN,-40,0])
    plt.title('mW, N = 1024')
    plt.annotate('main-lobe', xy=(0,-10), xytext=(-200, -5), fontsize=16, arrowprops=(dict(facecolor='black', width=2, headwidth=6, shrink=0.01)))
    plt.annotate('highest side-lobe', xy=(32,-13), xytext=(100, -10), fontsize=16, arrowprops=(dict(facecolor='black', width=2, headwidth=6, shrink=0.01)))
    
    
    M = 64
    N = 512
    hN = N/2     
    hM = M/2
    fftbuffer = np.zeros(N)
    mX1 = np.zeros(N)
    
    plt.figure(1, figsize=(7.5, 3.5))
    fftbuffer[hN-hM:hN+hM]=np.hanning(M)
    plt.subplot(2,1,1)
    plt.plot(np.arange(-hN, hN), fftbuffer, 'b', lw=1.5)
    plt.axis([-hN, hN, 0, 1.1])
    
    
    X = fft(fftbuffer)
    mX = 20*np.log10(abs(X)) 
    mX1[:hN] = mX[hN:]
    mX1[N-hN:] = mX[:hN]      
    
    plt.subplot(2,1,2)
    plt.plot(np.arange(-hN, hN), mX1-max(mX), 'r', lw=1.5)
    plt.axis([-hN,hN,-80,0])
    
    plt.tight_layout()
    plt.savefig('hanning.png')
    plt.show()

##########################
#You can put the code that calls the above functions down here
if __name__ == "__main__":
    
    f1 = 50
    f2 = 200
    f3 = 300
    f4 = 500
    N = 900
    N_Zero = 124
    Fs = 1000
    n = np.arange(0,N)
   # x = np.cos(2*np.pi*f1*n/Fs) + np.cos(2*np.pi*f2*n/Fs) + np.cos(2*np.pi*f3*n/Fs) + np.cos(2*np.pi*f4*n/Fs) 
    x1 = np.cos(2*np.pi*f1*n/Fs + 0.3 * np.pi) +  np.cos(2*np.pi*f2*n/Fs)
    #x1 = np.random.rand(N) * np.cos(2*np.pi*f1*n/Fs)
    #x1 = np.ones(len(x1) + N_Zero)
    # Zero padding
    x = np.zeros(len(x1) + N_Zero)
    x[0:len(x1)] = x1
    print len(x)
    
    M = len(x)
    w = get_window('hanning', M)
    fftbuffer = x*w
   #fftbuffer = np.roll(fftbuffer, (fftbuffer.size+1)/2)
    #X = fft(fftbuffer)                                      # compute FFT
    X = getDFT(fftbuffer,len(x))
    absX = abs(X)                                      # compute ansolute value of positive side
    absX[absX<np.finfo(float).eps] = np.finfo(float).eps    # if zeros add epsilon to handle log
    mX = 20 * np.log10(absX)                                # magnitude spectrum of positive frequencies in dB   
    mX[mX<-120] = 0
    pX = np.unwrap(np.angle(X))                        # unwrapped phase spectrum of positive frequencies [:hN]
    #pX = np.angle(X)  
    N = len(fftbuffer)
    
    
    plt.figure()
    plt.subplot(311)
    plt.plot(fftbuffer)
    plt.axis([0, len(fftbuffer), min(fftbuffer), max(fftbuffer)])
    plt.subplot(312)
    #plt.plot(mX)
    plt.plot(float(Fs)*np.arange(mX.size/2+1)/float(N), mX[:N/2+1], 'r')
    plt.axis([0, Fs/2.0, min(mX), max(mX)])
    plt.title ('magnitude spectrum: mX')
    plt.ylabel('amplitude (dB)')
    plt.xlabel('frequency (Hz)')

    #mX, pX = dftAnal(yfilt, w, N)
    print str(N/2+1)
    plt.subplot(313)
    plt.plot(float(Fs)*np.arange(pX.size/2+1)/float(N), pX[:N/2+1], 'c')
    plt.axis([0, Fs/2.0, min(pX), max(pX)])
    plt.title ('phase spectrum: pX')
    plt.ylabel('phase (radians)')
    plt.xlabel('frequency (Hz)')


    fftbuffer = np.roll(fftbuffer, (fftbuffer.size+1)/2)
    X = getDFT(fftbuffer,len(x))
    absX = abs(X)                                      # compute ansolute value of positive side
    absX[absX<np.finfo(float).eps] = np.finfo(float).eps    # if zeros add epsilon to handle log
    mX = 20 * np.log10(absX)                                # magnitude spectrum of positive frequencies in dB   
    mX[mX<-120] = 0
    pX = np.unwrap(np.angle(X))                        # unwrapped phase spectrum of positive frequencies [:hN]

    plt.figure()
    plt.subplot(311)
    plt.plot(fftbuffer)
    plt.axis([0, len(fftbuffer), min(fftbuffer), max(fftbuffer)])
    plt.subplot(312)
    #plt.plot(mX)
    plt.plot(float(Fs)*np.arange(mX.size/2+1)/float(N), mX[:N/2+1], 'r')
    plt.axis([0, Fs/2.0, min(mX), max(mX)])
    plt.title ('magnitude spectrum: mX')
    plt.ylabel('amplitude (dB)')
    plt.xlabel('frequency (Hz)')

    #mX, pX = dftAnal(yfilt, w, N)
    print str(N/2+1)
    plt.subplot(313)
    plt.plot(float(Fs)*np.arange(pX.size/2+1)/float(N), pX[:N/2+1], 'c')
    plt.axis([0, Fs/2.0, min(pX), max(pX)])
    plt.title ('phase spectrum: pX')
    plt.ylabel('phase (radians)')
    plt.xlabel('frequency (Hz)')
    
    
    
    M = len(x1)
    w = get_window('hanning', M)
    fftbuffer = x1*w
    X = getDFT(x1,len(x1))
    absX = abs(X)                                      # compute ansolute value of positive side
    absX[absX<np.finfo(float).eps] = np.finfo(float).eps    # if zeros add epsilon to handle log
    mX = 20 * np.log10(absX)                                # magnitude spectrum of positive frequencies in dB   
    mX[mX<-120] = 0
    pX = np.unwrap(np.angle(X))                        # unwrapped phase spectrum of positive frequencies [:hN]
    N = len(fftbuffer)
    
    plt.figure()
    plt.subplot(311)
    plt.plot(fftbuffer)
    plt.axis([0, len(fftbuffer), min(fftbuffer), max(fftbuffer)])
    plt.subplot(312)
    #plt.plot(mX)
    plt.plot(float(Fs)*np.arange(mX.size/2+1)/float(N), mX[:N/2+1], 'r')
    plt.axis([0, Fs/2.0, min(mX), max(mX)])
    plt.title ('magnitude spectrum: mX')
    plt.ylabel('amplitude (dB)')
    plt.xlabel('frequency (Hz)')

    #mX, pX = dftAnal(yfilt, w, N)
    print str(N/2+1)
    plt.subplot(313)
    plt.plot(float(Fs)*np.arange(pX.size/2+1)/float(N), pX[:N/2+1], 'c')
    plt.axis([0, Fs/2.0, min(pX), max(pX)])
    plt.title ('phase spectrum: pX')
    plt.ylabel('phase (radians)')
    plt.xlabel('frequency (Hz)')

