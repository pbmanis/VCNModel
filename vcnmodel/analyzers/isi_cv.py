import numpy as np
from dataclasses import dataclass

@dataclass
class Regularity():
    pass
    
def isi_cv(splist:list, binwidth:float=1.0, t0:float=0.0, t1:float=300.0, tgrace:float=0.0):
    """ compute the cv and regularity according to Young et al., J. Neurophys, 60: 1, 1988.
        Analysis is limited to isi's starting at or after t0 but before t1, and ending completely
        before t1 + tgrace(to avoid end effects). t1 should correspond to the
        the end of the stimulus
        Version using a list of numpy arrays for cvisi
    """
    cvisit = np.arange(0, t1, binwidth) # build time bins
    cvisi = [[]]*len(cvisit)
    for i in range(0, len(splist)): # for all the traces
        if len(splist[i]) < 2: # need at least 2 spikes
            # printf(f"Trial {i:d} fewer than 2 spikes")
            continue
        isib = np.floor(np.array(splist[i])[0:-2]/binwidth) # begining spike times for each interval
        isii = np.diff(splist[i]) # associated intervals
        for j in range(0, len(isib)): # loop over spikes
            if (splist[i][j] < t0) or (splist[i][j] > t1+tgrace) or (splist[i][j+1] > t1+tgrace):  # first spike of the pair is outside the window
                continue
            # print('spike, t1, tgrace', splist[i][j+1], t1, tgrace)
            if splist[i][j+1] > (t1 + tgrace): # second spike is after stimulus ends
                print('not in grace')
                continue
            cvisi[int(isib[j])] = np.append(cvisi[int(isib[j])], isii[j]) # and add the isi in that bin
    cvm = np.array([]) # set up numpy arrays for mean, std and time for cv analysis
    cvs = np.array([])
    cvt = np.array([])
    for i in range(0, len(cvisi)): # for each entry (possible bin)
        c = cvisi[i]
        # print(len(c))
        if len(c) >= 3: # require 3 spikes in a bin for statistics
            cvm = np.append(cvm, np.mean(c))
            cvs = np.append(cvs, np.std(c))
            cvt = np.append(cvt, i*binwidth)
    return(cvisit, cvisi, cvt, cvm, cvs)

def firing_rate(spikes):
    """
    Rate of the spike train.
    """
    if len(spikes) < 2:
        return np.nan
    return (len(spikes) - 1) / (spikes[-1] - spikes[0])


def CV(spikes):
    """
    Coefficient of variation.
    """
    if spikes == []:
        return np.nan
    ISI = diff(spikes)  # interspike intervals
    return std(ISI) / mean(ISI)

if __name__ == "__main__":

    """
    Test the calculation
    """
    import matplotlib.pyplot as mpl

    N = 100
    reps = 50
    rng = np.random.default_rng()
    spikes = rng.exponential(scale=5.0, size=(reps, N))
    T1 = [np.cumsum(s) for s in spikes]
    print(np.mean(T1[0]))
    print(f"Firing rate (first trial): {firing_rate(T1[0]):.3f}")
    for i in range(len(T1)):
        print(f"max time for trial: {i:d} = {np.max(T1[i]):.2f} ms   nspikes = {len(T1[i]):d}")
    cvisit, cvisi, cvt, cvm, cvs = isi_cv(T1, binwidth=2, t0 = 10, t1 = 210, tgrace = 25)
    print("len of isis: ", len(cvisi), len(cvisit))
    print("mean: ", np.mean(cvs/cvm))
    mpl.plot(cvt, cvs, 'bo-')
    mpl.plot(cvt, cvm, 'rs-')
    mpl.show()