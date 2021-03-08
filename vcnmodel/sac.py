from dataclasses import dataclass
import numpy as np
import pyqtgraph as pg
from pylibrary.plotting import pyqtgraph_plothelpers as pgh

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
# Note: calling the routine with letters as the first argument will generate the plots
# from Figure 2 of Louage et al. J. Neurophysiol. 2004.

# Python version, adapted from Matlab version, 6/16/2015 pbm.


@dataclass
class SACPars:
    twin: float = 5.0
    binw: float = 0.05
    delay: float = 0.0
    dur: float = 1000.0
    ddur: float = 10.0
    ntestrep: int = 20
    baseper: float = 4/3.
    stimdur: float = 1000


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
            binw = 0.15
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

    def SAC(self, X, pars):
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
        elif not isinstance(pars, SACPars):
            raise ValueError(
                f"SAC: parameters must be dict or dataclass of type SACPars"
            )

        N = len(X)
        maxn = 10000000  # pre-allocation for spike trains 

        yc = np.nan * np.zeros(maxn)
        n = 1
        # window data
        Y = [[]] * N
        spcount = np.zeros(N)
        for i in range(N):
            Y[i] = X[i][
                np.where(
                    (X[i] >= self.SPars.delay)
                    & (X[i] < (self.SPars.delay + self.SPars.dur))
                )
            ]
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
                            (xd >= (-self.SPars.twin - self.SPars.binw / 2.0))
                            & (xd <= (self.SPars.twin + self.SPars.binw / 2.0))
                        )
                    ]  # limit window
                    yc[ns : (ns + len(xd))] = xd  # store
                    ns += len(xd)  # keep track of storage location
        y = yc[~np.isnan(yc)]  # clean it up.

        # now calculate the normalized histogram.
        # normalization is N*(N-1)*r^2*deltat*D
        # where N is the number of presentations (lx), r is the rate (sp/sec),
        # deltat is the bin width for the histogram (sec), and D is the duration of
        # the trace.
        # normalization goes to 1 as the time gets larger...
        #
        rate = np.sum(spcount) / (self.SPars.ntestrep * self.SPars.dur / 1000.0)
        #        print'Mean firing rate: %8.1f' % rate
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
                -self.SPars.ddur, self.SPars.ddur, num=2 * int(self.SPars.ddur / self.SPars.binw)
            ),
            density=False,
        )
        yh = yh / nfac  # to convert to rate, spikes/second
        return yh, bins

    def XAC(self, X, Y, pars):
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
                f"SAC: parameters must be dict or dataclass of type SACPars"
            )
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
                # if m == n:
                #     continue  # skip identical trains
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
                -self.SPars.ddur, self.SPars.ddur, num=2 * int(self.SPars.ddur / self.SPars.binw)
            ),
            density=False,
        )
        yh = yh / nfac  # to convert to rate, spikes/second
        return yh, bins

    def diffcorr(self, X, Y, pars):
        ysac, binsac = self.SAC(X, pars)
        yxac, binxac = self.XAC(X, Y, pars)
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
        yh, bins = sac.SAC(X, sac.SPars)
        sach = pg.PlotCurveItem(
            bins, yh, stepMode=True, fillLevel=0, brush=(255, 0, 255, 255), pen=None
        )
        layout.getPlot((i, 1)).addItem(sach)
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
