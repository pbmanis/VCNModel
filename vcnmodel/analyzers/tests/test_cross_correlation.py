import numpy as np
from elephant import spike_train_correlation as ESTC
from elephant.conversion import BinnedSpikeTrain
from elephant.spike_train_generation import homogeneous_poisson_process
import quantities as pq 
from vcnmodel.analyzers import spikestatistics as SPKS  # from Brian
from vcnmodel.analyzers import reverseCorrelation as RC

"""
Test different cross-correlation methods (and compare)
"""

plot_flag = True

st1 = homogeneous_poisson_process(rate=100.0 * pq.Hz, t_start=0.0 * pq.s, t_stop=10.0 * pq.s)
st2 = st1 - 0.005*pq.s

# print(dir(st1))
st2 = homogeneous_poisson_process(rate=50.0 * pq.Hz, t_start=0.0 * pq.s, t_stop=10.0 * pq.s)
# st1 = np.linspace(0., 1, 100)*pq.s
# st2 = st1 + 0.005*pq.s
# standard analysis parameters:
bw = 1.0*1e-3*pq.s
width = 30.0*1e-3*pq.s
# use simple spike correlator for analysis
cc_simple = RC.reverse_correlation(st1, st2, binwidth=bw,
    corrwindow=[-width, width])

print('revcorr simple: ', cc_simple)
# use Elephant for analysis:
# cc_matrix = ESTC.correlation_coefficient(BinnedSpikeTrain([st1, st2], binsize=bw*pq.ms))
# print("ccmatrix: ", cc_matrix[0, 1])
# binned_st1 = BinnedSpikeTrain([st1], bin_size=bw*pq.ms)
# binned_st2 = BinnedSpikeTrain([st2], bin_size=bw*pq.ms)
# # print(dir(st1))
# # print('st1 times: ', st1.times)
# # exit()
# cc_hist = ESTC.cross_correlation_histogram(
#     binned_st1,
#     binned_st2,
#     window=[-width, width],
#     border_correction=False,
#     binary=False,
#     kernel=None,
#     method="memory",
# )

# print("cchist: ", cc_hist[0].times.magnitude)
# print("cchist mag: ", cc_hist[0].sampling_period.magnitude)

# use Brian for analysis:
cc_spks = SPKS.correlogram(
        st2.times*pq.ms, st1.times*pq.ms, width=2*width, bin_width=bw, T=None)

print('cc: ', cc_spks)
print(len(cc_spks))


if plot_flag:
    from matplotlib import pyplot as mpl
    f, ax = mpl.subplots(2, 1)
    # mpl.bar(
    #     x=cc_hist[0].times.magnitude,
    #     height=[c[0] for c in cc_hist[0][:, 0].magnitude],
    #     width=cc_hist[0].sampling_period.magnitude,
    #     linewidth=1,
    # )
    ax[0].eventplot([st1.times, st2.times], linelengths=0.75)
    t = np.linspace(-width, width, len(cc_spks))
    ax[1].plot(t, cc_spks, "r-")
    ax[1].plot(t[int(len(t)/2):], cc_simple[0], 'b-')
    # mpl.xlabel("time (" + str(cc_hist[0].times.units) + ")")
    mpl.ylabel("cross-correlation histogram")
    mpl.axis("tight")
    mpl.show()
