# import datetime
import functools
import time
from dataclasses import dataclass
from typing import Union
import numpy as np
import scipy
import pyqtgraph as pg
import pyximport
from numba import njit, jit
from numba.typed import List
from pylibrary.plotting import pyqtgraph_plothelpers as pgh

pyximport.install()
from vcnmodel.analyzers import sac_cython

"""
Shuffled autocorrelation function
Based on Louage et al., 2004
X is an array of N responses x length(t) points in time. The times are in
msec. Since the number of spikes in X may vary from trial to trial,
X is assumed to be a cell array of spike times.
The autocorrelation is calculated for every event in X that is not within
is not within a single spike train.
twin is the time window for the autocorrelation to be calculated (in
msec)
binw is the bin width for the histogram, in msec.
delay is the delay to the start of measurement in X; in msec.
dur is the duration of the measurement window in X; in msec.
(these last two parameters allow you to pass a response and only analyze
the portion that is relevant).

y is the output autocorrelation function. It consists of a list of
the correlation intervals, not corrected at all.

yh is the histogram bin amplitudes, with bin width hx.
This histogram is normalized.
mr is the mean rate

10/5/04 Paul B. Manis, Ph.D.
10/5/04 some fixes for an spike trains.
Note: calling the routine with letters as the first argument will
generate the specific plots
Replicates Figure 2 of Louage et al. J. Neurophysiol. 2004.

Python version, adapted from Matlab version, 6/16/2015 pbm.

Updated 12/13/2021, using numba and cython.
as cython is faster, we stick with that.
see sac_cython for implementations of sac and xac in cython.
SAC tested against original python and numba - all produce
the same result. (No numba version for xac was coded).
"""


@dataclass
class SACPars:
    twin: float = 0.005  # time window for sac calculation
    binw: float = 0.00005  # sac calculation bin width
    delay: float = 0.0  # delay into spike train to start calculation
    dur: float = 1.0  # duration of calulcation window in spike train
    displayDuration: float = 0.010  # display duration
    maxn: int = 100000000  # number of potential points in the SAC
    nrep: int = 20  # number of repetitions for test cases
    baseper: float = 4 / 3.0  # baseperiod for test cases
    stimdur: float = 1.0   # stimulus duration for test cases


@dataclass
class SACResult:
    SAC_peak : float = np.nan  # max of sac
    SAC_minimum: float = np.nan  # min of sac
    SAC_HW: float = np.nan  # half-width of SAC peak around 0 time
    n_spikes : int = 0  # number of spikes
    
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


# @time_func
@jit(nopython=False, parallel=False, cache=True)
def nb_SAC_Calc(X, event_lengths, twin, binw, maxn):
    yc = np.nan * np.zeros(maxn)
    ns = 0
    N = len(X)
    for n in range(N):  # for each trial
        for m in range(1, N):  # against each other trial, except:
            if m == n:
                continue  # skip identical trains
            for mi in range(event_lengths[m]):  # cross correlate all spikes in Yi
                xd = X[n][:event_lengths[n]] - X[m][mi]  # distance in time between the two spikes
                xd = xd[
                    np.where((xd >= (-twin - binw / 2.0)) & (xd <= (twin + binw / 2.0)))
                ]  # limit window
                nd = len(xd)
                yc[ns : ns + nd] = xd
                ns += nd
    y = yc[~np.isnan(yc)]  # clean it up.
    return y, ns


# @time_func
@njit(parallel=False, cache=True)
def nb_XAC_Calc(X, Y, event_lengths_x, event_lengths_y, twin, binw,
        maxn):
    M = len(X)
    N = len(Y)

    yc = np.nan * np.zeros(maxn)  # preallocate for spike trains
    ns = 0
    for n in range(N):  # for each trial
        for m in range(M):  # against each other trial, except:
            for mi in range(event_lengths_y):
                xd = X[n][:event_lengths_x[n]] - Y[m][mi]
                # all differences from i'th spike in m'th trial to spikes in nth trial
                xd = xd[
                    np.where((xd >= (-twin - binw / 2.0)) & (xd <= (twin + binw / 2.0)))
                ]  # limit window
                nd = len(xd)
                yc[ns : (ns + nd)] = xd  # store
                ns += nd  # keep track of storage location
    y = yc[~np.isnan(yc)]  # clean it up.
    return y, ns


# @time_func
def c_SAC_Calc(
    X: Union[list, np.ndarray], event_lengths:np.array, twin: float, binw: float,
        maxn: int,
):
    """
    Interface to cython version of SAC calculation.
    """
    y = np.full(maxn, np.nan)  # need to size array
    ns = 0
    event_lengths = np.array(event_lengths).squeeze()  # make sure is a 1-D np array

    sac_cython.sac_cython(
        X,  # data array (input)
        event_lengths,
        twin,  # time win single float input
        binw,
        y,  # result SAC (output)
        ns,
    )
    y = y[~np.isnan(y)]
    return y, ns
        

class SAC(object):
    def __init__(self):
        pass

    def makeTests(self, test, digitize=True):
        """
        Set up for the individual tests
        Parameters
        ----------
        test : str
            Letter for test, according to the panels in Louage et al. 2004
        """

        self.SPars = SACPars(
            twin=0.005,
            binw=0.00005,
            delay=0.0,
            dur=1.0,
            nrep=20,
            displayDuration=0.050,
            baseper=0.001*(4.0 / 3.0),
        )

        X = [None] * self.SPars.nrep

        if test == "A":
            for i in range(self.SPars.nrep):
                X[i] = np.arange(
                    self.SPars.baseper, self.SPars.stimdur, self.SPars.baseper
                )  # all identical

        elif test == "B":
            for i in range(self.SPars.nrep):
                X[i] = np.arange(
                    self.SPars.baseper, self.SPars.stimdur, self.SPars.baseper
                )  # all identical
                X[i] = X[i] + np.random.normal(0.0, 0.00008, len(X[i]))

        elif test == "C":
            for i in range(self.SPars.nrep):
                X[i] = np.arange(
                    self.SPars.baseper, self.SPars.stimdur, self.SPars.baseper
                )  # all identical
                X[i] = X[i] + np.random.normal(0.0, 0.00008, len(X[i]))
                n = len(X[i])
                sig = self.SPars.stimdur * np.sort(np.random.rand(int(n / 4)))
                X[i] = np.sort(np.append(X[i], sig))

        elif test == "D":
            for i in range(self.SPars.nrep):
                X[i] = np.arange(
                    self.SPars.baseper, self.SPars.stimdur, self.SPars.baseper
                )  # all identical
                X[i] = X[i] + np.random.normal(0.0, 0.000170, len(X[i]))

        elif test == "E":
            for i in range(self.SPars.nrep):
                bx = np.random.permutation(
                    np.arange(
                        self.SPars.baseper, self.SPars.stimdur, self.SPars.baseper
                    )
                )
                bx = bx[: int(len(bx) / 2)]
                X[i] = sorted(bx)  # elimnate half by random selection
                X[i] = X[i] + np.random.normal(0.0, 0.000170, len(X[i]))

        elif test == "F":
            X = [None] * self.SPars.nrep
            for i in range(self.SPars.nrep):
                X[i] = self.SPars.stimdur * np.sort(
                    np.random.rand(120)
                )  # corresponds to about 120 s/s
            self.SPars.displayDuration = 0.100

        elif test == "G":
            self.SPars.binw = 0.00015
            sig = self.SPars.stimdur * np.sort(np.random.rand(120))
            for i in range(self.SPars.nrep):
                X[i] = sig  # same signal every time
            self.SPars.displayDuration = 100

        # if digitize:  # digitize spike times to 25 usec bins
        #     t = np.arange(0., self.SPars.dur, 50e-6)
        #     X = t[np.digitize(X, t)]
        return X

    @time_func
    def py_SAC_Calc(
        self, X: np.ndarray, event_lengths: np.ndarray, twin: float, binw: float,
        maxn:int):
        """
        Python version of SAC calculation (for reference; this is slow)
        """
        N = X.shape[0]

        yc = np.nan * np.zeros(maxn)
        n = 1
        t0 = -twin - binw / 2.0
        t1 = twin + binw / 2.0
        # window data
        ns = 0
        for n in range(N):  # for each trial
            for m in range(1, N):  # against each other trial, except:
                if m == n:
                    continue  # skip identical trains
                for i in range(event_lengths[m]):  # cross correlate against all spikes in Yi
                    xd = (
                        X[n][:event_lengths[n]] - X[m][i]
                    )  # all differences from i'th spike in m'th trial to spikes in nth trial
                    xd = xd[
                        np.where(
                            (xd >= t0) & (xd < t1)
                        )
                    ]  # limit window
                    try:
                        yc[ns : (ns + len(xd))] = xd  # store
                    except:
                        print("py_SAC_Calc: extending yc array")
                        yc.append(np.nan*np.zeros(maxn))  # double array size
                        yc[ns : (ns + len(xd))] = xd  # store
                        # print('len yc ns, lenxd : ', len(yc), ns, len(xd))
                    
                    ns += len(xd)  # keep track of storage location
        y = yc[~np.isnan(yc)]  # clean it up.
        return y, ns



    def SAC(self, X, pars, engine:str="cython"):
        """
        Compute the SAC on X, given the parameter set

        Parameters
        ----------
        pars : dict or dataclass of SACPars
            Minimially must specify twin, binw, delay and dur

        """
        # now the SAC calculations.
        if isinstance(pars, dict):  # if dict, create dataclass version
            self.SPars = SACPars(
                twin=pars["twin"],
                binw=pars["binw"],
                delay=pars["delay"],
                dur=pars["dur"],
                displayDuration=pars["displayDuration"],
                maxn=pars["maxn"],
                
            )
        elif isinstance(pars, SACPars):
            self.SPars = pars
        else:
            raise ValueError(
                "SAC: parameters must be dict or dataclass of type SACPars"
            )

        # X, the list of spike times, may be a list of lists, which forms a "ragged" array
        # However, numba and cython work only with "regular" rectangular
        # arrays.
        # So, we reformat X into a rectangular array with NaNs 
        # following the end of the array of spike times for each trail
        # To simplify computations, we window the events here before passing
        # to the calculations.
        # We also make a 1-D array that holds the number of events in each trial.

        N = len(X)
        event_lengths = np.full(N, 0, dtype=int)
        event_lengths = [len(x) for x in X]
        XC = np.nan*np.zeros((N, int(np.max(event_lengths))))
        for i in range(N):
            okev = np.array(X[i][np.argwhere((X[i]>=self.SPars.delay) & (X[i] < (self.SPars.delay+self.SPars.dur)))]).flatten()
            XC[i,:len(okev)] = okev
            event_lengths[i] = len(okev)
        if engine in [ "python", "Python"]:
            y, ns = self.py_SAC_Calc(
                XC,
                event_lengths=event_lengths,
                twin=self.SPars.twin,
                binw=self.SPars.binw,
                maxn=self.SPars.maxn,
            )
        elif engine == "numba":
            y, ns = nb_SAC_Calc(
                XC,
                event_lengths=event_lengths,
                twin=self.SPars.twin,
                binw=self.SPars.binw,
                maxn=self.SPars.maxn,
            )
        elif engine == "cython":
            y, ns = c_SAC_Calc(
                XC,
                event_lengths=event_lengths,
                twin=self.SPars.twin,
                binw=self.SPars.binw,
                maxn=self.SPars.maxn,
            )
        else:
            raise ValueError(f"SAC engine specification was invalid: got {engine:s} but only 'python', 'numba' and 'cython' adre permitted")

        SACR = SACResult()
        if not np.isnan(np.sum(y)) and np.sum(y) > 0.0:
            SACR.SAC_peak = np.max(y)  # max of sac
            SACR.SAC_minimum = np.min(y) # min of sac
            SACR.SAC_HW = 0. # half-width of SAC peak around 0 time
            SACR.n_spikes = event_lengths  # number of spikes
        else:
            SACR.SAC_peak = np.nan
            SACR.SAC_minimum = np.nan
            SACR.n_spikes = 0
            
        return y, event_lengths, SACR

    def SAC_with_histo(self, X, pars, engine="cython", binsize=0.0, dither=0.0):
        if binsize > 0.0:  # digitize spike times to 25 usec bins
            t = np.arange(0., pars["dur"], binsize)
            X = t[np.digitize(X, t)]
        if dither > 0.0:
            X = [x + np.random.uniform(-dither, dither, size=len(x)) for x in X]
        y, spcount, SAC_R = self.SAC(X, pars=pars, engine=engine)
        yh, bins = self.SAC_make_histogram(y, spcount)
        return yh, bins

    def SAC_make_histogram(self, y, spcount):
        """
        Calculate the normalized histogram.
        normalization is N*(N-1)*r^2*deltat*D
        where:
            N is the number of presentations
            r is the rate (sp/sec),
            deltat is the bin width for the histogram (sec)
            and D is the duration of
        the trace.
        normalization goes to 1 as the time gets larger...
        """
        rate = np.mean(spcount) / self.SPars.dur
        N = len(spcount)  # number of presentations
        nfac = (
            N
            * (N - 1)
            * rate
            * rate
            * self.SPars.binw
            * self.SPars.dur
        )  # correction factor
        yh, bins = np.histogram(
            y,
            bins=np.linspace(
                -self.SPars.displayDuration,
                self.SPars.displayDuration,
                num=2 * int(self.SPars.displayDuration / self.SPars.binw),
            ),
            density=False,
        )
        bins = (bins[1:] + bins[:-1])/2.0
        yh = yh / nfac  # to convert to rate, spikes/second
        return yh, bins

    def SAC_measures(self, yh, bins, twin=0.003):
        """
        Measure the central peak CI and half-width
        from the SAC histogram
        
        Parameters
        ----------
        yh : the histogram (get from SAC_with_histo or SAC_make_histogramt)
        bins: the bins for the histogram (used for time window)
        twin : time window (in seconds) on either side of the central peak. 
            Pass the stimulus period for regular (SAM or phase-locked) stimuli
        
        Returns
        -------
        CI : peak value of CI
        HW : half-widths of the CI (and side lobes if present)
            This array has 4 elements for each peak:
                the width (in samples)
                the height (true value)
                left intersection point of horizontal line at half height (in samples)
                right intersection pont of horizontal ine at half height (in samples)
        FW : Full-width of CI, (and side lobes if present)
            same as HW, except for full height
        """
        if np.isnan(np.sum(yh)):
            return np.nan, [np.nan], [np.nan], [np.nan]# cannot compute, no spikes.
        win_data = np.where((-twin <= bins) & (bins <= twin))[0]
        win_data = yh[win_data]
        CI = np.max(win_data)
        peaks, _ = scipy.signal.find_peaks(win_data)
        center = [int((len(win_data)-1)/2)]  # only analyze the center peak
        half_widths = scipy.signal.peak_widths(win_data, center, rel_height=0.5)
        full_widths = scipy.signal.peak_widths(win_data, center, rel_height=1.0)
        return CI, peaks, half_widths, full_widths
        
        
        
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
        maxn: int=20000000,
    ):

        N = len(Y)
        M = len(X)

        # pre-allocation for spike trains diffs to avoid extending array repeatedly

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
        rate_x = np.sum(spike_count_x) / (self.SPars.nrep * self.SPars.dur)
        rate_y = np.sum(spike_count_y) / (self.SPars.nrep * self.SPars.dur)
        #        print'Mean firing rate: %8.1f' % rate
        nfac = (
            N
            * M
            * rate_x
            * rate_y
            * (self.SPars.binw)
            * (self.SPars.dur)
        )  # correction factor

        yh, bins = np.histogram(
            y,
            bins=np.linspace(
                -self.SPars.displayDuration,
                self.SPars.displayDuration,
                num=2 * int(self.SPars.displayDuration / self.SPars.binw),
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
    We also test all versions of the SAC calculation (different computation engines)
    with timings, and assert that they give identical results for
    all versions.
    """
    sac = SAC()
    tests = ["A", "B", "C", "D", "E", "F", "G"]
    ymax = {"A": 30, "B": 6, "C": 6, "D": 6, "E": 6, "F": 6, "G": 56}
    win = pgh.figure(title="AN Inputs")
    layout = pgh.LayoutMaker(
        cols=2, rows=len(tests), win=win, labelEdges=True, ticks="talbot"
    )
    maxn = 25000000
    if len(tests) <= 2:
        win.resize(600, 300)
    else:
        win.resize(600, 125 * len(tests))
    for i, t in enumerate(tests):
        X = sac.makeTests(t)
        
        yhp, binsp = sac.SAC_with_histo(X, sac.SPars, engine="python", dither=0.)
        yhn, binsn = sac.SAC_with_histo(X, sac.SPars, engine="numba", dither=0.)
        yhc, binsc = sac.SAC_with_histo(X, sac.SPars, engine="cython", dither=0.)
        # assert np.array_equal(yhc, yhn)
        # assert np.array_equal(yhp, yhc)
        # assert np.array_equal(yhp, yhn)  # make sure all results are the same

        SAC_with_histo = pg.PlotDataItem(
            binsp, yhp, stepMode=True, fillLevel=0, brush=(128, 128, 128, 255), pen=None
        )
        sacn = pg.PlotDataItem(binsn, yhn, stepMode=True, brush=None, pen=pg.mkPen("m")
        )
        sacc = pg.PlotDataItem(binsc, yhc, stepMode=True, brush=None, pen=pg.mkPen("r")
        )
        layout.getPlot((i, 1)).addItem(SAC_with_histo)
        layout.getPlot((i, 1)).addItem(sacn)
        layout.getPlot((i, 1)).addItem(sacc)
        layout.getPlot((i, 1)).setXRange(-0.005, 0.005)
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
            layout.getPlot((i, 0)).setXRange(0.0, 0.010)
        else:
            layout.getPlot((i, 0)).setXRange(0.0, 1.0)

    pgh.show()
