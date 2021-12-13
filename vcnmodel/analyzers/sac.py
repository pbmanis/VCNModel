# import datetime
import functools
import time
from dataclasses import dataclass

import numpy as np
import pyqtgraph as pg
import pyximport
from numba import njit
from numba.typed import List
from pylibrary.plotting import pyqtgraph_plothelpers as pgh

pyximport.install()
from vcnmodel.analyzers import sac_cython

# Shuffled autocorrelation function
# Based on Louage et al., 2004
# X is an array of N responses x length(t) points in time. The times are in
# msec. Since the number of spikes in X may vary from trial to trial,
# X is assumed to be a cell array of spike times.
# The autocorrelation is calculated for every event in X that is not within
# is not within a single spike train.
# twin is the time window for the autocorrelation to be calculated (in
# msec)
# binw is the bin width for the histogram, in msec.
# delay is the delay to the start of measurement in X; in msec.
# dur is the duration of the measurement window in X; in msec.
# (these last two parameters allow you to pass a response and only analyze
# the portion that is relevant).
#
# y is the output autocorrelation function. It consists of a list of
# the correlation intervals, not corrected at all.
#
# yh is the histogram bin amplitudes, with bin width hx.
# This histogram is normalized.
# mr is the mean rate
#
# 10/5/04 Paul B. Manis, Ph.D.
# 10/5/04 some fixes for an spike trains.
# Note: calling the routine with letters as the first argument will
# generate the specific plots
# Replicates Figure 2 of Louage et al. J. Neurophysiol. 2004.

# Python version, adapted from Matlab version, 6/16/2015 pbm.

# Updated 12/13/2021, using numba and cython.
# as cython is faster, we stick with that.
# see sac_cython for implementations of sac and xac in cython.
# SAC tested against original python and numba - all produce
# the same result. (No numba version for xac was coded).


@dataclass
class SACPars:
    twin: float = 5.0
    binw: float = 0.05
    delay: float = 0.0
    dur: float = 1000.0
    ddur: float = 10.0
    ntestrep: int = 20
    baseper: float = 4 / 3.0
    stimdur: float = 1000


def time_func(func):
    """
    A decorator to show execution time of functions.
    Place inside (after) winprint if using
    Output is to terminal.
    """

    @functools.wraps(func)
    def wrapper_timer(self, *args, **kwargs):
        print(f"Starting : {func.__name__!r}")
        start_time = time.perf_counter()  # 1
        value = func(self, *args, **kwargs)
        end_time = time.perf_counter()  # 2
        run_time = end_time - start_time  # 3
        print(f"Finished {func.__name__!r} in {run_time:.4f} secs")
        return value

    return wrapper_timer


@time_func
@njit(parallel=False, cache=True)
def nb_SAC_Calc(X, twin, binw, delay, dur):
    N = len(X)
    maxn = 10000000  # pre-allocation for spike trains

    yc = np.nan * np.zeros(maxn)
    spcount = np.zeros(N)
    ns = 0
    for n in range(N):  # for each trial
        n_index = [np.where((X[n] >= delay) & (X[n] < (delay + dur)))][0][0]
        spcount[n] = len(n_index)
        for m in range(N):  # against each other trial, except:
            if m == n:
                continue  # skip identical trains
            m_index = [np.where((X[m] >= delay) & (X[m] < (delay + dur)))][0][0]
            for mi in m_index:  # cross correlate all spikes in Yi
                xd = X[n][n_index] - X[m][mi]  # distance in time between the two spikes
                xd = xd[
                    np.where((xd >= (-twin - binw / 2.0)) & (xd <= (twin + binw / 2.0)))
                ]  # limit window
                nd = len(xd)
                yc[ns : ns + nd] = xd
                ns += nd
    y = yc[~np.isnan(yc)]  # clean it up.
    return y, spcount


@time_func
@njit(parallel=False, cache=True)
def nb_XAC_Calc(X, Y, twin, binw, delay, dur):
    M = len(X)
    N = len(Y)
    maxn = 10000000  # pre-allocation for spike trains

    yc = np.nan * np.zeros(maxn)
    spike_count_x = np.zeros(M)
    spike_count_y = np.zeros(N)
    ns = 0
    for n in range(N):  # for each trial
        n_index = [np.where((X[i] >= delay) & (X[i] < (delay + dur)))][0][0]
        spike_count_x[n] = len(n_index)
        for m in range(M):  # against each other trial, except:
            m_index = [np.where((X[m] >= delay) & (X[m] < (delay + dur)))][0][0]
            spike_count_y[m] = len(m_index)
            for mi in m_index:  # cross correlate all spikes in the windows
                xd = X[n][n_index] - X[m][mi]
                # all differences from i'th spike in m'th trial to spikes in nth trial
                xd = xd[
                    np.where((xd >= (-twin - binw / 2.0)) & (xd <= (twin + binw / 2.0)))
                ]  # limit window
                nd = len(xd)
                yc[ns : (ns + nd)] = xd  # store
                ns += nd  # keep track of storage location
    y = yc[~np.isnan(yc)]  # clean it up.
    return y, spike_count_x, spike_count_y


class SAC(object):
    def __init__(self):
        pass

    def makeTests(self, test):
        """
        Set up for the individual tests
        Parameters
        ----------
        test : str
            Letter for test, according to the panels in Louage et al. 2004
        """

        self.SPars = SACPars(
            twin=5.0,
            binw=0.05,
            delay=0.0,
            dur=1000.0,
            ntestrep=20,
            ddur=10.0,
            baseper=4.0 / 3.0,
        )

        X = [None] * self.SPars.ntestrep

        if test == "A":
            for i in range(self.SPars.ntestrep):
                X[i] = np.arange(
                    self.SPars.baseper, self.SPars.stimdur, self.SPars.baseper
                )  # all identical

        elif test == "B":
            for i in range(self.SPars.ntestrep):
                X[i] = np.arange(
                    self.SPars.baseper, self.SPars.stimdur, self.SPars.baseper
                )  # all identical
                X[i] = X[i] + np.random.normal(0.0, 0.08, len(X[i]))

        elif test == "C":
            for i in range(self.SPars.ntestrep):
                X[i] = np.arange(
                    self.SPars.baseper, self.SPars.stimdur, self.SPars.baseper
                )  # all identical
                X[i] = X[i] + np.random.normal(0.0, 0.08, len(X[i]))
                n = len(X[i])
                sig = self.SPars.stimdur * np.sort(np.random.rand(int(n / 4)))
                X[i] = np.sort(np.append(X[i], sig))

        elif test == "D":
            for i in range(self.SPars.ntestrep):
                X[i] = np.arange(
                    self.SPars.baseper, self.SPars.stimdur, self.SPars.baseper
                )  # all identical
                X[i] = X[i] + np.random.normal(0.0, 0.170, len(X[i]))

        elif test == "E":
            for i in range(self.SPars.ntestrep):
                bx = np.random.permutation(
                    np.arange(
                        self.SPars.baseper, self.SPars.stimdur, self.SPars.baseper
                    )
                )
                bx = bx[: int(len(bx) / 2)]
                X[i] = sorted(bx)  # elimnate half by random selection
                X[i] = X[i] + np.random.normal(0.0, 0.170, len(X[i]))

        elif test == "F":
            X = [None] * self.SPars.ntestrep
            for i in range(self.SPars.ntestrep):
                X[i] = self.SPars.stimdur * np.sort(
                    np.random.rand(120)
                )  # corresponds to about 120 s/s
            self.SPars.ddur = 100.0

        elif test == "G":
            self.SPars.binw = 0.15
            sig = self.SPars.stimdur * np.sort(np.random.rand(120))
            for i in range(self.SPars.ntestrep):
                X[i] = sig  # same signal every time
            self.SPars.ddur = 100

        return X

    @time_func
    def SAC_Calc(
        self, X: np.ndarray, twin: float, binw: float, delay: float, dur: float
    ):
        N = len(X)
        maxn = 10000000  # pre-allocation for spike trains

        yc = np.nan * np.zeros(maxn)
        n = 1
        # window data
        Y = [[]] * N
        spcount = np.zeros(N)
        for i in range(N):
            Y[i] = X[i][np.where((X[i] >= delay) & (X[i] < (delay + dur)))]
            spcount[i] = len(Y[i])
        ns = 0
        for n in range(N):  # for each trial
            for m in range(N):  # against each other trial, except:
                if m == n:
                    continue  # skip identical trains
                for i in range(len(Y[m])):  # cross correlate against all spikes in Yi
                    xd = (
                        Y[n] - Y[m][i]
                    )  # all differences from i'th spike in m'th trial to spikes in nth trial
                    xd = xd[
                        np.where(
                            (xd >= (-twin - binw / 2.0)) & (xd <= (twin + binw / 2.0))
                        )
                    ]  # limit window
                    yc[ns : (ns + len(xd))] = xd  # store
                    ns += len(xd)  # keep track of storage location
        y = yc[~np.isnan(yc)]  # clean it up.
        return y, spcount

    @time_func
    def c_SAC_Calc(
        self, X: np.ndarray, twin: float, binw: float, delay: float, dur: float
    ):
        y = np.full(10000000, np.nan)  # need to size array
        spcount = np.zeros(np.array(X).shape[0], dtype=np.int)
        ns = 0
        sac_cython.sac_cython(
            np.array(X),  # data array (input)
            self.SPars.twin,  # time win single float input
            self.SPars.binw,
            self.SPars.delay,
            self.SPars.dur,
            y,  # result SAC (output)
            spcount,  # spike count, int, output
            ns,
        )
        y = y[~np.isnan(y)]
        return y, spcount

    def SAC(self, X, pars, engine="cython"):
        """
        Compute the SAC on X, given the parameter set

        Parameters
        ----------
        pars : dict or dataclass of SACPars
            Minimially must specify twin, binw, delay and dur

        """
        # now the SAC calculations.
        if isinstance(pars, dict):
            self.SPars = SACPars(
                twin=pars["twin"],
                binw=pars["binw"],
                delay=pars["delay"],
                dur=pars["dur"],
            )
        elif isinstance(pars, SACPars):
            self.SPars = pars
        else:
            raise ValueError(
                "SAC: parameters must be dict or dataclass of type SACPars"
            )

        if engine in ["python", "Python"]:
            y, spcount = self.SAC_Calc(
                X,
                twin=self.SPars.twin,
                binw=self.SPars.binw,
                delay=self.SPars.delay,
                dur=self.SPars.dur,
            )
        elif engine == "numba":
            y, spcount = nb_SAC_Calc(
                np.array(X),
                twin=self.SPars.twin,
                binw=self.SPars.binw,
                delay=self.SPars.delay,
                dur=self.SPars.dur,
            )
        elif engine == "cython":
            y, spcount = self.c_SAC_Calc(
                X,
                twin=self.SPars.twin,
                binw=self.SPars.binw,
                delay=self.SPars.delay,
                dur=self.SPars.dur,
            )

        return y, spcount

    def SAC_with_histo(self, X, pars, engine="cython"):
        y, spcount = self.SAC(X, pars=pars, engine=engine)
        yh, bins = self.SAC_make_histogram(y, spcount)
        return yh, bins

    def SAC_make_histogram(self, y, spcount):
        # now calculate the normalized histogram.
        # normalization is N*(N-1)*r^2*deltat*D
        # where N is the number of presentations (lx), r is the rate (sp/sec),
        # deltat is the bin width for the histogram (sec), and D is the duration of
        # the trace.
        # normalization goes to 1 as the time gets larger...
        #
        rate = np.sum(spcount) / (self.SPars.ntestrep * self.SPars.dur / 1000.0)
        #        print'Mean firing rate: %8.1f' % rate
        N = len(X)
        nfac = (
            N
            * (N - 1)
            * rate
            * rate
            * (self.SPars.binw / 1000.0)
            * (self.SPars.dur / 1000.0)
        )  # correction factor
        yh, bins = np.histogram(
            y,
            bins=np.linspace(
                -self.SPars.ddur,
                self.SPars.ddur,
                num=2 * int(self.SPars.ddur / self.SPars.binw),
            ),
            density=False,
        )
        yh = yh / nfac  # to convert to rate, spikes/second
        return yh, bins

    def XAC(self, X, Y, pars, engine="python"):
        """
        Cross correlation SAC
        """
        if isinstance(pars, dict):
            self.SPars = SACPars(
                twin=pars["twin"],
                binw=pars["binw"],
                delay=pars["delay"],
                dur=pars["dur"],
            )
        elif not isinstance(pars, SACPars):
            raise ValueError(
                "SAC: parameters must be dict or dataclass of type SACPars"
            )

        if engine in ["python", "Python"]:
            y, spcount_x, spcount_y = self.XAC_Calc(
                X,
                Y,
                twin=self.SPars.twin,
                binw=self.SPars.binw,
                delay=self.SPars.delay,
                dur=self.SPars.dur,
            )
        elif engine == "cython":
            y, spcount_x, spcount_y = self.c_XAC_Calc(
                X,
                Y,
                twin=self.SPars.twin,
                binw=self.SPars.binw,
                delay=self.SPars.delay,
                dur=self.SPars.dur,
            )
        return y, spcount_x, spcount_y

    @time_func
    def c_XAC_Calc(
        self,
        X: np.ndarray,
        Y: np.ndarray,
        twin: float,
        binw: float,
        delay: float,
        dur: float,
    ):
        y = np.full(10000000, np.nan)  # need to size array
        spcount_x = np.zeros(np.array(X).shape[0], dtype=np.int)
        spcount_y = np.zeros(np.array(Y).shape[0], dtype=np.int)
        ns = 0
        sac_cython.xac_cython(
            np.array(X),  # data array (input)
            np.array(Y),
            self.SPars.twin,  # time win single float input
            self.SPars.binw,
            self.SPars.delay,
            self.SPars.dur,
            y,  # result SAC (output)
            spcount_x,  # spike count, int, output
            spcount_y,
            ns,
        )
        y = y[~np.isnan(y)]
        return y, spcount_x, spcount_y

    def XAC_with_histo(self, X, Y, pars, engine="cython"):
        y, spcount_x, spcount_y = self.XAC(X, Y, pars=pars, engine=engine)
        yh, bins = self.XAC_make_histogram(spcount_x, spcount_y)
        return yh, bins

    @time_func
    def XAC_Calc(
        self,
        X: np.ndarray,
        Y: np.ndarray,
        twin: float,
        binw: float,
        delay: float,
        dur: float,
    ):

        N = len(Y)
        M = len(X)

        # pre-allocation for spike trains diffs to avoid extending array repeatedly
        maxn = 10000000

        yc = np.nan * np.zeros(maxn)
        n = 1
        # window data
        Xa = [[]] * M
        Ya = [[]] * N
        spike_count_x = np.zeros(M)
        spike_count_y = np.zeros(N)
        for i in range(M):
            Xa[i] = X[i][
                np.where(
                    (X[i] >= self.SPars.delay)
                    & (X[i] < (self.SPars.delay + self.SPars.dur))
                )
            ]
            spike_count_x[i] = len(Xa[i])
        for i in range(N):
            Ya[i] = Y[i][
                np.where(
                    (Y[i] >= self.SPars.delay)
                    & (Y[i] < (self.SPars.delay + self.SPars.dur))
                )
            ]
            spike_count_y[i] = len(Ya[i])
        ns = 0
        for n in range(N):  # for each trial
            for m in range(M):  # against each other trial, except:
                for i in range(len(Xa[m])):  # cross correlate against all spikes in Yi
                    xd = (
                        Ya[n] - Xa[m][i]
                    )  # all differences from i'th spike in m'th trial to spikes in nth trial
                    xd = xd[
                        np.where(
                            (xd >= (-self.SPars.twin - self.SPars.binw / 2.0))
                            & (xd <= (self.SPars.twin + self.SPars.binw / 2.0))
                        )
                    ]  # limit window
                    yc[ns : (ns + len(xd))] = xd  # store
                    ns += len(xd)  # keep track of storage location
        y = yc[~np.isnan(yc)]  # clean it up.
        return y, spike_count_x, spike_count_y

    def XAC_make_histogram(self, spike_count_x, spike_count_y):
        # now calculate the normalized histogram.
        # normalization is N*(N-1)*r^2*deltat*D
        # where N is the number of presentations (lx), r is the rate (sp/sec),
        # deltat is the bin width for the histogram (sec), and D is the duration of
        # the trace.
        # normalization goes to 1 as the time gets larger...
        #
        rate_x = np.sum(spike_count_x) / (self.SPars.ntestrep * self.SPars.dur / 1000.0)
        rate_y = np.sum(spike_count_y) / (self.SPars.ntestrep * self.SPars.dur / 1000.0)
        #        print'Mean firing rate: %8.1f' % rate
        nfac = (
            N
            * M
            * rate_x
            * rate_y
            * (self.SPars.binw / 1000.0)
            * (self.SPars.dur / 1000.0)
        )  # correction factor

        yh, bins = np.histogram(
            y,
            bins=np.linspace(
                -self.SPars.ddur,
                self.SPars.ddur,
                num=2 * int(self.SPars.ddur / self.SPars.binw),
            ),
            density=False,
        )
        yh = yh / nfac  # to convert to rate, spikes/second
        return yh, bins

    def diffcorr(self, X, Y, pars):
        ysac, binsac = self.SAC_with_histo(X, pars)
        yxac, binxac = self.XAC_with_histo(X, Y, pars)
        diffcorr = ysac - yxac
        return diffcorr, binsac


if __name__ == "__main__":

    """
    Test the SAC method, regenerating Figure 2 of Louage et al. J. Neurophys, 2004
    """
    sac = SAC()
    tests = ["A", "B", "C", "D", "E", "F", "G"]
    ymax = {"A": 30, "B": 6, "C": 6, "D": 6, "E": 6, "F": 6, "G": 56}
    win = pgh.figure(title="AN Inputs")
    layout = pgh.LayoutMaker(
        cols=2, rows=len(tests), win=win, labelEdges=True, ticks="talbot"
    )
    if len(tests) <= 2:
        win.resize(600, 300)
    else:
        win.resize(600, 125 * len(tests))
    for i, t in enumerate(tests):
        X = sac.makeTests(t)
        yhp, binsp = sac.SAC_with_histo(X, sac.SPars, engine="python")
        yhn, binsn = sac.SAC_with_histo(X, sac.SPars, engine="numba")
        yh, bins = sac.SAC_with_histo(X, sac.SPars, engine="cython")
        assert np.array_equal(yh, yhn)
        assert np.array_equal(yh, yhp)
        assert np.array_equal(yhp, yhn)  # make sure all results are the same
        SAC_with_histo = pg.PlotCurveItem(
            binsn, yhn, stepMode=True, fillLevel=0, brush=(255, 0, 255, 255), pen=None
        )
        layout.getPlot((i, 1)).addItem(SAC_with_histo)
        layout.getPlot((i, 1)).setXRange(-5.0, 5.0)
        layout.getPlot((i, 1)).setYRange(0, ymax[t])
        size = 4
        Y = [[]] * len(X)
        for j in range(len(X)):
            Y[j] = (j + 1) * np.ones(len(X[j]))
        raster = pg.ScatterPlotItem(
            x=np.array(X).flatten(),
            y=np.array(Y).flatten(),
            pen="w",
            brush="b",
            size=size,
            pxMode=True,
        )
        layout.getPlot((i, 0)).addItem(raster)
        if t not in ["F", "G"]:
            layout.getPlot((i, 0)).setXRange(0.0, 10.0)
        else:
            layout.getPlot((i, 0)).setXRange(0.0, 100.0)

    pgh.show()
