from vcnmodel.util.user_tester import UserTester
import numpy as np
import neo
import elephant.spike_train_correlation as ESTC
import elephant.conversion as EC

# import elephant.spike_train_generation as ESTG
import quantities as pq
from vcnmodel.analyzers import reverse_correlation as RC  # our local version

"""
Test different cross-correlation methods (and compare)
We test our "local" calculation in reverse_correlation.py with a reference
package in elephant.
The Brian package (in spikestatistics.py) makes different assumptions, and is
not comparable.
"""


class regular_process(object):
    def __init__(self, start, interval, duration, offset=0):
        spkt = offset * pq.s + np.arange(start, duration, interval) * pq.s
        self.times = neo.SpikeTrain(spkt, t_stop=duration, units="s")


def test_cross_correlation():
    CrossCorrelationTester(key="CrossCorrelation")


def compute_cc_data(spikes="coherent", package="local"):
    if spikes not in ["coherent", "offset", "uncorrelated"]:
        raise NotImplementedError(
            f"The spike train format for testing cc is not known: '{str(spikes):s}'"
        )
    if package not in ["local", "elephant", "elephant_sttc"]:
        raise NotImplementedError(
            f"The package must be 'local' or 'elephant', but got: {str(package):s}"
        )
    tstart = 0
    tstop = 1.0
    st1 = regular_process(
        start=0.010, duration=tstop, interval=0.010, offset=0,
    )  # all identical
    # st1 = ESTG.homogeneous_poisson_process(rate=20.0 * pq.Hz, t_start=0.0 * pq.s, t_stop=50.0 * pq.s)
    if spikes == "coherent":
        st2 = st1
    elif spikes == "offset":
        st2 = regular_process(
            start=0.010, duration=tstop, interval=0.010, offset=0.005,
        )  # al
    elif spikes == "uncorrelated":
        # st2 = ESTG.homogeneous_poisson_process(rate=20.0 * pq.Hz, t_start=0.0 * pq.s, t_stop=50.0 * pq.s)
        st2 = regular_process(
            start=(0.1 / np.pi), duration=tstop, interval=(0.1 / 2.713), offset=0
        )
    # standard analysis parameters:
    bw = 0.10 * 1e-3 * pq.s
    width = 30.0 * 1e-3 * pq.s

    if package == "local":
        cc_result, ncc = RC.reverse_correlation(
            st1.times, st2.times, binwidth=bw, corrwindow=[-width, width],
        )
        t = np.linspace(-width, width, len(cc_result), endpoint=False)
    elif package == "elephant":
        nbins = int(width / bw)
        bst_i = EC.BinnedSpikeTrain(
            spiketrains=st1.times,
            bin_size=bw,
            n_bins=None,
            t_start=tstart * pq.s,
            t_stop=tstop * pq.s,
            tolerance=1e-8,
            sparse_format="csr",
        )
        bst_j = EC.BinnedSpikeTrain(
            spiketrains=st2.times,
            bin_size=bw,
            n_bins=None,
            t_start=tstart * pq.s,
            t_stop=tstop * pq.s,
            tolerance=1e-8,
            sparse_format="csr",
        )
        cc_result, lags = ESTC.cross_correlation_histogram(
            bst_i,
            bst_j,
            window=(-nbins, nbins),
            border_correction=False,
            binary=False,
            kernel=None,
            method="speed",
            cross_correlation_coefficient=True,
        )
        t = np.linspace(-width, width, len(cc_result), endpoint=False)
    elif package == "elephant_sttc":
        nbins = int(width / bw)
        bst_i = EC.BinnedSpikeTrain(
            spiketrains=st1.times,
            bin_size=bw,
            n_bins=None,
            t_start=tstart * pq.s,
            t_stop=tstop * pq.s,
            tolerance=1e-8,
            sparse_format="csr",
        )
        bst_j = EC.BinnedSpikeTrain(
            spiketrains=st2.times,
            bin_size=bw,
            n_bins=None,
            t_start=tstart * pq.s,
            t_stop=tstop * pq.s,
            tolerance=1e-8,
            sparse_format="csr",
        )
        sttc_wins = np.arange(0.0002, 0.006, 0.00001)*pq.s
        cc_result = np.zeros_like(sttc_wins)
        for i, dt in enumerate(sttc_wins):  # run over a range of windows
            cc_result[i] = ESTC.spike_time_tiling_coefficient(
                st1.times, st2.times, dt=dt,
            )
        t = sttc_wins
    return cc_result, t, st1, st2


class CrossCorrelationTester(UserTester):
    def __init__(self, key):
        self.packages = ["local", "elephant", "elephant_sttc"]
        self.spikerelationships = ["coherent", "offset", "uncorrelated"]
        UserTester.__init__(self, key)

    def assert_test_info(self, *args, **kwds):
        try:
            super(CrossCorrelationTester, self).assert_test_info(*args, **kwds)
        finally:
            pass

    def run_test(self):
        for j, pkg in enumerate(self.packages):
            print(f"Package: {pkg:s}")
            for i, spkt in enumerate(self.spikerelationships):
                print(f"Spike timing: {spkt:s}")
                self.ccs, self.t, self.st1, self.st2 = compute_cc_data(
                    spikes=spkt, package=pkg
                )
                self.show_result(i, j)
            # if self.audit:
            #     self.show_result()
        mpl.show()

    def show_result(self, index: int = 0, pkgindex=0):
        width = 0.1
        i = index * 2
        j = pkgindex
        if i == 0 and j == 0:
            self.f, self.ax = mpl.subplots(
                2 * len(self.spikerelationships), len(self.packages), figsize=(6, 8)
            )
        self.ax[i, j].eventplot(
            [self.st1.times.tolist(), self.st2.times.tolist()],
            linelengths=[0.75, 0.75],
            lineoffsets=[0, 1],
            colors=["b", "m"],
        )
        self.ax[i + 1, j].plot(self.t, self.ccs.squeeze(), "r-")
        self.ax[i + 1, j].set_ylabel("CCF")
        self.ax[i, j].set_ylabel("Spikes")
        if i == 0:
            self.ax[i, j].set_title(self.packages[pkgindex])
        if i == 5:
            self[i, j] / set_xlabel("T (s)")
        # mpl.axis("tight")
        # mpl.show()


if __name__ == "__main__":
    from matplotlib import pyplot as mpl

    cc = CrossCorrelationTester("cc test")
    cc.run_test()
