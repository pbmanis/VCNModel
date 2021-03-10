import numpy as np
import scipy.stats
import matplotlib.pyplot as mpl

class VS():
    def __init__(self):
        pass

    def vector_strength(self, spikes, freq):
        """
        Calculate vector strength and related parameters from a spike train, for the specified frequency
        :param spikes: Spike train, in sec.
        :param freq: Stimulus frequency in Hz
        :return: a dictionary containing:
    
            r: vector strength
            n: number of spikes
            R: Rayleigh coefficient
            p: p value (is distribution not flat?)
            ph: the circularized spike train over period of the stimulus freq, freq, in radians
            d: the "dispersion" computed according to Ashida et al., 2010, etc.
        """
    
        per = 1/freq # c
        ph = 2*np.pi*np.fmod(spikes, per)/(per) # convert to radians within a cycle
        sumcos = np.sum(np.cos(ph))
        sumsin = np.sum(np.sin(ph))
        mean_phase = np.arctan2(sumsin,sumcos)
        sumc2 = sumcos**2
        sums2 = sumsin**2
        n = len(spikes)
        vs = (1./n)*np.sqrt(sumc2+sums2)  # standard vector strength computation
        R = 2*n*vs*vs  # Raleigh coefficient
        Rp = np.exp(-n*vs*vs)  # p value for n > 50 (see Ashida et al. 2010).
        d = np.sqrt(2.*(1-vs))/(2*np.pi*freq)
        self.spikes = spikes
        return{'r': vs, 'n': n, 'R': R, 'p': Rp, 'ph': ph, 'd': d}
    
    def vs_print(self, d):
        print(f"{'Vector Strength':>24s}: {d['r']:12.4f}")
        print(f"{'Spike Count':>24s}: {d['n']:12d}")
        print(f"{'mean phase (radians)':>24s}: {scipy.stats.circmean(d['ph'])/(2*np.pi):12.4f}")
        print(f"{'SD phase':>24s}: {scipy.stats.circstd(d['ph'])/(2.*np.pi):12.4f}")
        print(f"{'Rayleigh':>24s}: {d['R']:12.4f}")
        print(f"{'p value':>24s}: {d['p']:12.4e}")
        print(f"{'d':>24s}: {d['d']:12.4f}")
  #

if __name__ == '__main__':
    
    V = VS()
    # perfect vs
    freq = 100. # Hz
    spikes = np.arange(0.1, 10., 0.010)
    print("'perfect'")
    d = V.vector_strength(spikes, freq)
    V.vs_print(d)
    # mpl.plot(1.-np.mod(spikes, 1./freq))
    # mpl.show()
    print('-'*80)
    print("Uniform random")
    spikes = np.random.uniform(0.1, 10., 5000)
    d = V.vector_strength(spikes, freq)
    V.vs_print(d)
    print('-'*80)
    print(f"5% doublets 5 msec")
    spikes = np.arange(0.1, 10., 0.010)
    delr = int(len(spikes)*30./100.)
    for n in range(delr):
        ni = int(n*len(spikes)/delr)
        spikes = np.insert(spikes, ni+1, spikes[ni]+0.005+np.random.normal(scale=0.0005))
    d = V.vector_strength(spikes, freq)
    mpl.hist(d['ph'])
    mpl.show()
    V.vs_print(d)
    print('-'*80)
    