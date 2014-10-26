import numpy as np
import matplotlib.pyplot as plt

plt.figure(1, figsize=(9.5, 7))
N = 250
k0 = 80
k1 = 200
X = np.array([])
nv = np.arange(-N/2, N/2)
kv = np.arange(-N/2, N/2)
x =  np.cos(2*np.pi*k0/N*nv) #+ np.cos(2*np.pi*k1/N*nv) # Sampling rate / Frequency of signal = k0/N = Fs

plt.subplot(311)
plt.title('x; k = ' + str(k0) + '; N = ' + str(N))
plt.plot(nv, np.real(x),'b', lw=1.5)
plt.axis([-N/2,N/2-1,-1,1])
for k in kv:
	s = np.exp(1j*2*np.pi*k/N*nv)
	X = np.append(X, sum(x*np.conjugate(s)))

plt.subplot(312)
plt.title('magnitude spectrum: abs(X)')
plt.plot(kv, abs(X), 'ro', lw=1.5)
plt.axis([-N/2,N/2-1,0,N])

plt.subplot(313)
plt.title('phase spectrum: angle(X)')
plt.plot(kv, np.angle(X),'c', lw=1.5)
plt.axis([-N/2,N/2-1,-np.pi,np.pi])

plt.tight_layout()
plt.show()
