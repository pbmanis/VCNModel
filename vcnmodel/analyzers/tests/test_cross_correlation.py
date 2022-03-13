"""
Test different cross-correlation methods (and compare)
We test our "local" calculation in reverse_correlation.py with a reference
package in elephant.
The Brian package (in spikestatistics.py) makes different assumptions, and is
not comparable.
"""

import neo
import numpy as np
import matplotlib.pyplot as mpl
import elephant.spike_train_generation as ESTG
import quantities as pq
from elephant import conversion as EC
from elephant import spike_train_correlation as ESTC
import elephant.spike_train_generation as ESTG

from vcnmodel.analyzers import reverse_correlation as RC  # our local version
from vcnmodel.util.user_tester import UserTester
from vcnmodel.analyzers import sttc

STTC = sttc.STTC(seed=0, engine="numba")

class RegularProcess:
    """
    Gemerate a regular spike train process -determinstic - for testing
    """

    def __init__(self, start, interval, duration, offset=0):
        """
        Parameters
        ----------
        start: float
            start time for regular process
        interval: float
            intervals between events in process
        duration: float
            End time of the process (max time)
        offset: float
            offset time (added to all times in process)

        Returns
        -------
            nothing
         """
        spkt = offset * pq.s + np.arange(start, duration, interval) * pq.s
        self.times = neo.SpikeTrain(spkt, t_stop=duration, units="s")


def test_cross_correlation():
    """
    Invoke the cross-correlation tests
    """
    CrossCorrelationTester(key="CrossCorrelation")


def compute_cc_data(spikes="coherent", package="local"):
    """
    Compute the cross correlation for one set of spike trains using
    a particular package method

    Parameters
    ----------
    spikes: string
        The name of the type of spike trains to compare
    package: string
        The name of the method for testing correolation

    """
    if spikes not in ["coherent", "offset", "uncorrelated", "poisson"]:
        raise NotImplementedError(
            f"The spike train format for testing cc is not known: '{str(spikes):s}'"
        )
    if package not in ["local", "elephant", "local_sttc", "elephant_sttc"]:
        raise NotImplementedError(
            f"The package must be 'local' or 'elephant', but got: {str(package):s}"
        )
    np.random.seed(92311)  # force constant starting state
    # TODO: THis should use rng = numpy.random.default_rng(92311), but
    # then need to modify elephant poisson process to use the new
    # generator as well.

    tstart = 0
    if spikes == "poisson":
        tstop = 50.
    else:
        tstop = 1.0
        st1 = RegularProcess(
            start=0.010, duration=tstop, interval=0.010, offset=0,
        )  # all identical
    if spikes == "coherent":
        st2 = st1
    elif spikes == "offset":
        st2 = RegularProcess(
            start=0.010, duration=tstop, interval=0.010, offset=0.005,
        )  # al
    elif spikes == "uncorrelated":
        st2 = RegularProcess(
            start=(0.1 / np.pi), duration=tstop, interval=(0.1 / 2.713), offset=0
        )
    elif spikes == "poisson":
        st1 = ESTG.homogeneous_poisson_process(rate=20.0 * pq.Hz,
                    t_start=0.0 * pq.s, t_stop=tstop * pq.s)
        st2 = ESTG.homogeneous_poisson_process(rate=20.0 * pq.Hz,
                    t_start=0.0 * pq.s, t_stop=tstop * pq.s)
    st1 = neo.SpikeTrain(st1.times, t_stop=tstop, units='s')
    st2 = neo.SpikeTrain(st2.times, t_stop=tstop, units='s')
    # standard analysis parameters:
    bin_width = 0.10 * 1e-3 * pq.s
    width = 30.0 * 1e-3 * pq.s
    bst_i = EC.BinnedSpikeTrain(
        spiketrains=st1,
        bin_size=bin_width,
        n_bins=None,
        t_start=tstart * pq.s,
        t_stop=tstop * pq.s,
        tolerance=None,
        sparse_format="csr",
    )
    bst_j = EC.BinnedSpikeTrain(
        spiketrains=st2,
        bin_size=bin_width,
        n_bins=None,
        t_start=tstart * pq.s,
        t_stop=tstop * pq.s,
        tolerance=None,
        sparse_format="csr",
    )
    if package == "local":
        cc_result, ncc = RC.reverse_correlation(
            st1.times, st2.times, binwidth=bin_width, corrwindow=[-width, width],
        )
        cc_times = np.linspace(-width, width, len(cc_result), endpoint=False)
    elif package == "elephant":
        nbins = int(width / bin_width)
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
        cc_times = np.linspace(-width, width, len(cc_result), endpoint=False)
    elif package == "elephant_sttc":
        nbins = int(width / bin_width)
        sttc_wins = np.arange(0.000, 0.020, 0.001)
        cc_result = np.zeros_like(sttc_wins)
        for i, sttc_win in enumerate(sttc_wins):  # run over a range of windows
            cc_result[i] = ESTC.spike_time_tiling_coefficient(
                st1, #neo.SpikeTrain(st1.times, t_stop=tstop*pq.s),
                st2, #neo.SpikeTrain(st2.times, t_stop=tstop*pq.s),
                dt=sttc_win*pq.s,
            )
        cc_times = sttc_wins
    elif package == 'local_sttc':
        sttc_wins = np.arange(0.000, 0.020, 0.001)
        cc_result = np.zeros_like(sttc_wins)
        samplerate = 0.0001
        for i, sttc_win in enumerate(sttc_wins):  # run over a range of windows
            STTC.set_spikes(samplerate, st1.times, st2.times,
                            tilewindow=sttc_win)
            cc_result[i] = STTC.calc_sttc()
        cc_times = sttc_wins
    return cc_result, cc_times, st1, st2


class CrossCorrelationTester(UserTester):
    """
    Instantiate the unit test for cross correlations.
    """

    def __init__(self, key):
        """

        Parameters
        ----------
        key : string
            Key to identify this test.

        Returns
        -------
        None.

        """

        self.packages = ["local", "elephant", "local_sttc", "elephant_sttc"]
        self.spikerelationships = ["coherent", "offset", "uncorrelated", "poisson"]
        self.cc_times = None
        self.ccs = None
        self.st1 = None
        self.st2 = None
        self.fig = None
        self.axes = None

        UserTester.__init__(self, key)

    def assert_test_info(self, *args, **kwds):
        try:
            super(CrossCorrelationTester, self).assert_test_info(*args, **kwds)
        finally:
            pass

    def run_test(self, plotflag=False):
        for j, pkg in enumerate(self.packages):
            # print(f"Package: {pkg:s}")
            for i, spkt in enumerate(self.spikerelationships):
                # print(f"Spike timing: {spkt:s}")
                self.ccs, self.cc_times, self.st1, self.st2 = compute_cc_data(
                    spikes=spkt, package=pkg
                )
                if self.audit and plotflag:
                    self.show_result(i, j)
        if self.audit and plotflag: 
            mpl.show()
        info = {}

        return info

    def show_result(self, index: int = 0, pkgindex: int = 0):
        """
        Plot the results

        Parameters
        ----------
        index : int, optional
            Index into the spike relationship list. The default is 0.
        pkgindex : int, optional
            Index into the package list. The default is 0.

        Returns
        -------
        None.

        """
        from matplotlib import pyplot as mpl
        # width = 0.1
        i = index * 2
        j = pkgindex
        if i == 0 and j == 0:
            self.fig, self.axes = mpl.subplots(
                2 * len(self.spikerelationships), len(self.packages), figsize=(10, 10)
            )
        self.axes[i, j].eventplot(
            [self.st1.times.tolist(), self.st2.times.tolist()],
            linelengths=[0.75, 0.75],
            lineoffsets=[0, 1],
            colors=["b", "m"],
        )
        self.axes[i + 1, j].plot(self.cc_times, self.ccs.squeeze(), "r-")
        self.axes[i + 1, j].set_ylabel("CCF")
        self.axes[i, j].set_ylabel("Spikes")
        if i == 0:
            self.axes[i, j].set_title(self.packages[pkgindex])
        if i == 5:
            self.axes[i, j].set_xlabel("T (s)")
        # mpl.axis("tight")
        # mpl.show()


if __name__ == "__main__":

    cc = CrossCorrelationTester("CrossCorrelation")
    cc.run_test()
