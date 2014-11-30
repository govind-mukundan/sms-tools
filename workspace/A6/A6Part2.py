import os
import sys
import numpy as np
import math
from scipy.signal import get_window
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../software/models/'))
import utilFunctions as UF
import harmonicModel as HM
import stft

eps = np.finfo(float).eps

"""
A6Part2 - Segmentation of stable note regions in an audio signal

Complete the function segmentStableNotesRegions() to identify the stable regions of notes in a specific 
monophonic audio signal. The function returns an array of segments where each segment contains the 
starting and the ending frame index of a stable note.

The input argument to the function are the wav file name including the path (inputFile), threshold to 
be used for deciding stable notes (stdThsld), minimum allowed duration of a stable note (minNoteDur), 
number of samples to be considered for computing standard deviation (winStable), analysis window (window), 
window size (M), FFT size (N), hop size (H), error threshold used in the f0 detection (f0et), magnitude 
threshold for spectral peak picking (t), minimum allowed f0 (minf0) and maximum allowed f0 (maxf0). 
The function returns a numpy array of shape (k,2), where k is the total number of detected segments. 
The two columns in each row contains the starting and the ending frame indices of a stable note segment. 
The segments must be returned in the increasing order of their start times. 

In order to facilitate the assignment we have configured the input parameters to work with a particular 
sound, '../../sounds/sax-phrase-short.wav'. The code and parameters to estimate the fundamental frequency 
is completed. Thus you start from an f0 curve obtained using the f0Detection() function and you will use 
that to obtain the note segments. 

All the steps to be implemented in order to solve this question are indicated in segmentStableNotesRegions() 
as comments. These are the steps:

1. In order to make the processing musically relevant, the f0 values should be converted first from 
Hertz to Cents, which is a logarithmic scale. 
2. At each time frame (for each f0 value) you should compute the standard deviation of the past winStable 
number of f0 samples (including the f0 sample at the current audio frame). 
3. You should then apply a deviation threshold, stdThsld, to determine if the current frame belongs 
to a stable note region or not. Since we are interested in the stable note regions, the standard 
deviation of the previous winStable number of f0 samples (including the current sample) should be less 
than stdThsld. 
4. All the consecutive frames belonging to the stable note regions should be grouped together into 
segments. For example, if the indices of the frames corresponding to the stable note regions are 
3,4,5,6,12,13,14, we get two segments, first 3-6 and second 12-14. 
5. After grouping frame indices into segments filter/remove the segments which are smaller in duration 
than minNoteDur. Return the segment indexes.

Using '../../sounds/sax-phrase-short.wav', if you use all the default values of the input arguments in 
segmentStableNotesRegions(), then the resulting segment array should be 
array([[  8,  83],
       [ 93, 119],
       [128, 202],
       [207, 326],
       [360, 464],
       [471, 502]])

if you use all the default values except stdThsld = 30.0, then the resulting segment array should be 
array([[  8,  85],
       [ 88, 120],
       [123, 203],
       [206, 326],
       [360, 466],
       [469, 506]])
If you use all the default values except minNoteDur = 0.3, then the resulting segment array should be
array([[  8,  83],
       [128, 202],
       [207, 326],
       [360, 464]])

We also provide the function plotSpectogramF0Segments() to plot the f0 contour and the detected segments 
on the top of the spectrogram of the audio signal in order to visually analyse the outcome of your function. 

"""

def segmentStableNotesRegions(inputFile = '../../sounds/sax-phrase-short.wav', stdThsld=10, minNoteDur=0.1, 
                              winStable = 3, window='hamming', M=1024, N=2048, H=256, f0et=5.0, t=-100, 
                              minf0=310, maxf0=650):
    """
    Function to segment the stable note regions in an audio signal
    Input:
        inputFile (string): wav file including the path
        stdThsld (float): threshold for detecting stable regions in the f0 contour
        minNoteDur (float): minimum allowed segment length (note duration)  
        winStable (integer): number of samples used for computing standard deviation
        window (string): analysis window
        M (integer): window size used for computing f0 contour
        N (integer): FFT size used for computing f0 contour
        H (integer): Hop size used for computing f0 contour
        f0et (float): error threshold used for the f0 computation
        t (float): magnitude threshold in dB used in spectral peak picking
        minf0 (float): minimum fundamental frequency in Hz
        maxf0 (float): maximum fundamental frequency in Hz
    Output:
        segments (np.ndarray): Numpy array containing starting and ending frame indices of every 
                               segment.
    """
    fs, x = UF.wavread(inputFile)                               #reading inputFile
    w  = get_window(window, M)                                  #obtaining analysis window    
    f0 = HM.f0Detection(x, fs, w, N, H, t, minf0, maxf0, f0et)  #estimating F0
    # f0Detection splits the signal into "frames" of length H (hop size) and returns an f0 for each of these frames
    # So the number of frames in the signal is simply the length of f0

    ### Your code here

    # 1. convert f0 values from Hz to Cents (as described in pdf document)
    f0 = 1200 * np.log2(f0/55.0)
    numFrames = len(f0)
    
    #f0_sd = []
    # 2. create an array containing standard deviation of last winStable samples
#    for i in np.arange(winStable,numFrames - winStable -1):
#        print str(winStable) + str(i)
#        f0_sd[i - winStable] = np.std(f0[i - winStable : i])
    frameIndex = np.arange(winStable - 1, numFrames) # New index from winStable - 1 to num of frames
    print frameIndex
    # It's winStable - 1 because the SD includes current sample also. i.e. for winStable = 3, you need to use indexes 0 - 2 for first sample

    # Find a new array of standard deviations where each sample is the SD with the last winStable frames
    f0_sd = np.array(map(lambda i: np.std(f0[i - winStable + 1:i+1]),frameIndex))
    
    # 3. apply threshold on standard deviation values to find indices of the stable points in melody
    # To get the IDs of stable frames, you need to add the 
    stidx = np.where(f0_sd < stdThsld)[0] + winStable - 1 
    
    # 4. create segments of continuous stable points such that concequtive stable points belong to 
    #    same segment
    segments = groupConsecutiveRuns(stidx)

    # 5. apply segment filtering, i.e. remove segments with are < minNoteDur in length
    # Now segments is a list of arrays of the form [[startidx_1, endidx_1],[startidx_2, endidx_2]...]
    # We need to remove all segments that are not long enough
    minFrames = int(minNoteDur * fs / H)
    opseg = []
    for item in segments:
        if(filterShortSegments(item[0],item[1],minFrames) == True):
            print "Got segment: " + str(item)
            opseg.append(item)
    
    segments = np.array(opseg)
    # plotSpectogramF0Segments(x, fs, w, N, H, f0, segments)
    plotSpectogramF0Segments(x, fs, w, N, H, f0, segments)
    print str(segments)
    # return segments
    return segments

def groupConsecutiveRuns(xs):
    # Array_split splits the array into sub arrays at the given indexes
    # np.diff(a) - out[n] = a[n+1] - a[n]
    segments = np.array_split(xs, np.where(np.diff(xs) != 1)[0] + 1)
    print np.shape(segments)
    segsmall = []
    i = 0;
    # Now just get the start and end values of each segment
    for item in segments:
        print item
        segsmall.append([item[0],item[len(item) - 1]])
        i = i + 1
        
    return segsmall    

def filterShortSegments(i_start, i_end, minFrames):
    if(minFrames <= (i_end - i_start)):
        return True
    else:
        return False

## Use this function to plot the f0 contour and the estimated segments on the spectrogram
def plotSpectogramF0Segments(x, fs, w, N, H, f0, segments):
    """
    Code for plotting the f0 contour on top of the spectrogram
    """
    # frequency range to plot
    maxplotfreq = 1000.0    
    fontSize = 16

    fig = plt.figure()
    ax = fig.add_subplot(111)

    mX, pX = stft.stftAnal(x, fs, w, N, H)                      #using same params as used for analysis
    mX = np.transpose(mX[:,:int(N*(maxplotfreq/fs))+1])
    
    timeStamps = np.arange(mX.shape[1])*H/float(fs)                             
    binFreqs = np.arange(mX.shape[0])*fs/float(N)
    
    plt.pcolormesh(timeStamps, binFreqs, mX)
    plt.plot(timeStamps, f0, color = 'k', linewidth=5)

    for ii in range(segments.shape[0]):
        plt.plot(timeStamps[segments[ii,0]:segments[ii,1]], f0[segments[ii,0]:segments[ii,1]], color = '#A9E2F3', linewidth=1.5)        
    
    plt.autoscale(tight=True)
    plt.ylabel('Frequency (Hz)', fontsize = fontSize)
    plt.xlabel('Time (s)', fontsize = fontSize)
    plt.legend(('f0','segments'))
    
    xLim = ax.get_xlim()
    yLim = ax.get_ylim()
    ax.set_aspect((xLim[1]-xLim[0])/(2.0*(yLim[1]-yLim[0])))    
    plt.autoscale(tight=True) 
    plt.show()
    
##########################
#You can put the code that calls the above functions down here    
if __name__ == "__main__":
    segmentStableNotesRegions()