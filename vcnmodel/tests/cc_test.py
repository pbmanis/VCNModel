import matplotlib.pyplot as mpl
import numpy as np
import neo
from quantities import s, Hz, ms
from elephant.spike_train_generation import homogeneous_poisson_process
import elephant.spike_train_correlation as ESTC
from elephant.conversion import BinnedSpikeTrain

from vcnmodel import spikestatistics as SPKS

st1 = homogeneous_poisson_process(
      rate=100.0*Hz, t_start=0.0*s, t_stop=10.0*s)
st2 = homogeneous_poisson_process(
      rate=50.0*Hz, t_start=0.0*s, t_stop=10.0*s)

cc_matrix = ESTC.corrcoef(BinnedSpikeTrain([st1, st2], binsize=1*ms))
print(cc_matrix[0, 1])
binned_st1 = BinnedSpikeTrain([st1], binsize=1*ms)
binned_st2 = BinnedSpikeTrain([st2], binsize=1*ms)
# print(dir(st1))
print(st1.times)
# exit()
cc_hist =  ESTC.cross_correlation_histogram(
      binned_st1, binned_st2, window=[-30,30],
      border_correction=False,
      binary=False, kernel=None, method='memory')
print('cchist: ', cc_hist[0].times.magnitude)
# print(cc_hist[0][:, 0].magnitude)
print(cc_hist[0].sampling_period.magnitude)
mpl.bar(x=cc_hist[0].times.magnitude,
        height=[c[0] for c in cc_hist[0][:, 0].magnitude],
        width=cc_hist[0].sampling_period.magnitude, linewidth=1)
        
cc = SPKS.correlogram(
        st1.times, st2.times, width=60, bin=0.001, T=None
)
print('cc: ', cc)
print(len(cc))
t = np.linspace(-30., 30., len(cc))
mpl.plot(t, cc, 'r-')
mpl.xlabel('time (' + str(cc_hist[0].times.units) + ')')
mpl.ylabel('cross-correlation histogram')
mpl.axis('tight')
mpl.show()