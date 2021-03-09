import src.vcnmodel.analyzers.vector_strength as vector_strength
import numpy as np

plot_flag = False
VS = vector_strength.VectorStrength()
def compute_vs(freq=100., nsp=1000., sd=0.0, rg=None):
    if sd <= 1:
        sdn = (2.0/freq)+sd*rg.standard_normal(size=nsp)  +np.pi# just express in temrs of time
    else:  
        sdn = (1./freq) * rg.uniform(size=nsp)  # uniform across the interval
    spikes = np.cumsum(np.ones(nsp)*(1./freq))+sdn  # locked spike train with jitter
    phsp = 2*np.pi*freq*np.fmod(spikes, 1./freq)
    return(spikes, phsp)

def test_vs():
    from numpy.random import default_rng
    rg = default_rng(12345)
    freq = 100.
    nsp = 1000
    sd = 0.002 # in seconds, unles > 1 then is uniform
    x = np.array([0.,0.00005, 0.0001, 0.0002, 0.0003, 0.0005, 0.00075, 0.001, 0.002, 0.003, 0.004, 0.0050, 0.0075, 2])
    y = np.zeros_like(x)
    ph = np.zeros_like(x)
    vs = np.zeros_like(x)
    for i, sd in enumerate(x):
        spikes, phsp = compute_vs(freq, nsp, sd, rg)
        vsd = VS.vector_strength(spikes, freq)
        y[i] = vsd.circ_timeSD
        vs[i] = vsd.vs
        ph[i] = vsd.circ_phaseSD
    x[x==2] = 0.020

    if plot_flag:
        import matplotlib.pyplot as mpl
        fig, ax = mpl.subplots(3, 1)
        ax = ax.ravel()
        ax[0].plot(x, y, 'o')
        ax[0].plot([0, np.max(x)], [0, np.max(x)], '--k', alpha=0.5)
        ax[1].plot(x, vs, 'x-')
        ax[2].plot(x, ph, 's-')
        mpl.show()
        exit()

        fig, ax = mpl.subplots(2, 1)
        ax = ax.ravel()
        plot_ticks(ax[0], phsp, 'r-')
        plot_ticks(ax[1], spikes, 'b-')
        mpl.show()
        print('mean isi: ', np.mean(np.diff(spikes)), 'stim period: ', 1./freq)
        vsd = VS.vector_strength(spikes, freq)
        print('VS: ', vsd.vs, 'SD sec: ', vsd.circ_timeSD)
        print('circ_phase min max: ', np.min(vsd.circ_phase), np.max(vsd.circ_phase))
        bins = np.linspace(0, 2*np.pi, 36)
        mpl.hist(vsd.circ_phase, bins)
        mpl.show()
    
test_vs()
