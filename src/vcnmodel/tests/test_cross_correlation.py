import numpy as np
from elephant import spike_train_correlation as ESTC
from elephant.conversion import BinnedSpikeTrain
from elephant.spike_train_generation import homogeneous_poisson_process
from quantities import Hz, ms, s
from src.vcnmodel.analyzers import spikestatistics as SPKS

plot_flag = False

st1 = homogeneous_poisson_process(rate=100.0 * Hz, t_start=0.0 * s, t_stop=10.0 * s)
st2 = homogeneous_poisson_process(rate=50.0 * Hz, t_start=0.0 * s, t_stop=10.0 * s)

cc_matrix = ESTC.correlation_coefficient(BinnedSpikeTrain([st1, st2], bin_size=1 * ms))
print("ccmatrix: ", cc_matrix[0, 1])
binned_st1 = BinnedSpikeTrain([st1], bin_size=1 * ms)
binned_st2 = BinnedSpikeTrain([st2], bin_size=1 * ms)
# print(dir(st1))
# print('st1 times: ', st1.times)
# exit()
cc_hist = ESTC.cross_correlation_histogram(
    binned_st1,
    binned_st2,
    window=[-30, 30],
    border_correction=False,
    binary=False,
    kernel=None,
    method="memory",
)
print("cchist: ", cc_hist[0].times.magnitude)
print("cchist mag: ", cc_hist[0].sampling_period.magnitude)

cc = SPKS.correlogram(
        st1.times*ms, st2.times*ms, width=60.0*ms, bin_width=1.0*ms, T=None)
# )
# print('cc: ', cc)
# print(len(cc))

if plot_flag:
    from matplotlib import pyplot as mpl

    mpl.bar(
        x=cc_hist[0].times.magnitude,
        height=[c[0] for c in cc_hist[0][:, 0].magnitude],
        width=cc_hist[0].sampling_period.magnitude,
        linewidth=1,
    )
    t = np.linspace(-30.0, 30.0, len(cc))
    mpl.plot(t, cc, "r-")
    mpl.xlabel("time (" + str(cc_hist[0].times.units) + ")")
    mpl.ylabel("cross-correlation histogram")
    mpl.axis("tight")
    mpl.show()
