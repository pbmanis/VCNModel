#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

dt = 0.05  # msec


def generator(i0=0, dt=0.05, sigma0=2.0, fmod=0.2, tau=3.0, dur=10.0):
    
    # compute all times
    tb = np.linspace(0, dur*1000., int(dur*1000./dt))
    It = np.zeros_like(tb)
    # compute how    sigma varies over time:
    sigma = sigma0*(1.0+0.5*(1-np.cos(2.0*np.pi*(fmod/1000)*tb)))
    
    nrand = np.random.rand(tb.shape[0])-0.5
    sig = nrand*np.sqrt(2.0*(sigma*sigma)*dt/tau)
    It[0] = i0
    # brute now:
    for j in range(len(tb)-1):
        i = j + 1
        It[i] = It[i-1] + dt*(i0 - It[i-1])/tau + sig[i]
    return(tb, It)

if __name__ == '__main__':
    tb, It = generator()

    plt.plot(tb, It)
    plt.show()
