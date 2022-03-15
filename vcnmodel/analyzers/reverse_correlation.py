import datetime
import quantities as pq
from dataclasses import dataclass, field
from typing import List, Tuple, Union

import numpy as np

from pylibrary.tools import cprint as CP
import elephant.spike_train_correlation as ESTC
import elephant.conversion as EC
import neo
from vcnmodel.analyzers import spikestatistics as SPKS
from vcnmodel.analyzers import sttc as STTC

cprint = CP.cprint

"""
Compute the reverse correlation between two spike trains
For all of the spikes in st1, compute the time difference
with each spike in st2, find where the difference is within a
time window we are interested in, then get and return the indices
for all of those spikes.

We define 3 data structures first:
RevCorrPars are the parameters for the calculations, and inlcude
the stimulus and runinfo parameters from the simulations.
RevCorrData holds the results of the calculations.
SpikeData holds information about the individual postsynaptic spikes.

reverse_correlation is the python implementation

"""


def def_empty_np():
    return np.array(0)


def def_empty_list():
    return []


def def_empty_dict():
    return {}


@dataclass()
class RevCorrPars:
    ntrials: int = 1  # number of trials
    ninputs: int = 1  # number of inputs for the cell under study
    algorithm: str = "RevcorrSPKS"  # save the algorithm name
    # clip trace to avoid end effects
    min_time: float = 10.0  # msec to allow the system to settlt  # this window needs to be at least as long as minwin
    max_time: float = 250.0  # this window needs to be at least as long as maxwin
    binw: float = 0.1  # binwidth, ms
    minwin: float = -5.0  # start of revcorr display/calculation, ms relative to reference event
    maxwin: float = 2.5  # end of display/calculation relative to reference event
    amax: float = 0.0  # peak amplitude
    si: dict = field(default_factory=def_empty_dict)  # stimulus parameters
    ri: dict = field(default_factory=def_empty_dict)  # runInfon parameters


@dataclass()
class RevCorrData:
    C: list = field(default_factory=def_empty_list)  # using Brian 1.4 correlation
    CB: list = field(
        default_factory=def_empty_list
    )  # using elephant/ binned correlation
    CBT: list = field(
        default_factory=def_empty_list
    )  # using elephant/neo, binned correlation
    TC: float = 0.0  # list = field(default_factory=def_empty_list)
    st: np.array = field(default_factory=def_empty_np)  # spike times
    tx: np.array = field(default_factory=def_empty_np)  #
    ti: np.array = field(default_factory=def_empty_np)
    ti_avg: np.array = field(
        default_factory=def_empty_np
    )  # timebase for average revcorr
    sv_all: list = field(
        default_factory=def_empty_list
    )  # np.array = field(default_factory=def_empty_np)  #
    sv_avg: np.array = field(default_factory=def_empty_np)
    sv_trials: list = field(default_factory=def_empty_list)
    sites: np.array = field(default_factory=def_empty_np)
    nsp_avg: int = 0
    npost_spikes: int = 0
    npre_spikes: int = 0
    mean_pre_intervals: np.array = field(default_factory=def_empty_np)
    mean_post_intervals: float = 0.0
    max_coin_rate: float = 0
    participation: np.array = field(default_factory=def_empty_np)
    s_pair: float = 0.0
    ynspike: np.array = field(default_factory=def_empty_np)
    pre_w: list = field(default_factory=def_empty_list)
    pre_st: list = field(default_factory=def_empty_list)


@dataclass
class SpikeData:
    """
    Data class to hold information about each spike
    """

    trial: int = -1  # the trial the spike came from
    time_index: int = 0  # index into the time array for this spike
    dt: float = 0.025  # sample interval, msec
    start_win: float = -5.0
    end_win: float = 5.0
    waveform: np.array = None  # the waveform of this postspike, clipped
    prespikes: np.array = None  # time indices to pre spikes


def reverse_correlation(
    st1: Union[np.ndarray, List] = None,
    st2: Union[np.ndarray, List] = None,
    binwidth: float = 0.1,
    corrwindow: List = [
        -5.0,
        1.0,
    ],  # time window to examine correlation relative to st1
) -> Tuple[np.ndarray, int]:
    """
    Basic reverse correlation between two sets of event times
    Bins are set up so that the center bin straddles 0 time
    
    """

    if st1 is None or st2 is None:
        raise ValueError(
            "reverse_correlation: reference and comparator must be defined"
        )

    nbins = int((corrwindow[1] - corrwindow[0]) / binwidth) + 1
    bins = np.arange(
        corrwindow[0] - binwidth, corrwindow[1] + binwidth * 3, binwidth
    )  # bin centers
    bins -= binwidth / 2.0  # shift to let bins mark left edges
    xds = np.zeros_like(bins)  # (int((corrwindow[1] - corrwindow[0]) / binwidth))
    for i, sp in enumerate(st1):
        diff = st2 - sp
        # differences of spike times within the window
        v = diff[np.where((corrwindow[0] <= diff) & (diff <= corrwindow[1]))]
        if len(v) > 0:
            indxs = np.digitize(v, bins)
            # indxs = [int((vx-binwidth/2.) / binwidth) for vx in v]
            xds[indxs] = xds[indxs] + 1

    return xds, len(st1)  # return the n postsynaptic spikes

    # @trace_calls.time_func


def _count_spikes_in_window(data, trial, site, s, RCP, pre_w):
    an_i = data["Results"][trial]["inputSpikeTimes"][
        site
    ]  # input spike times for one input
    an_i = an_i[
        (an_i > RCP.min_time) & (an_i < RCP.max_time)
    ]  # restrict to those only within the response window
    an_i = an_i - s  # get relative latency from spike to it's inputs
    # print('ani: ', ani)
    pre_indx = np.asarray((an_i >= pre_w[0]) & (an_i <= pre_w[1])).nonzero()[0]
    pre_times = [an_i[k] for k in pre_indx]
    if len(pre_times) > 0:
        # print("pre times, prewindow: ", pre_times, pre_w)

        npre_i = len(pre_times)
    else:
        npre_i = 0
        pre_times = []
    return npre_i, pre_times


def revcorr(model_data:object=None, nbins: int = 0, revcorrtype: str = ""):
    maxtc = 0
    nspk_plot = 0
    data = model_data.data
    AR = model_data.AR
    RCP = model_data.RCP
    RCD = model_data.RCD
    post_intervals = []
    pre_intervals = [[] for x in range(RCP.ninputs)]
    start_time = datetime.datetime.now()
    if RCP.algorithm == "RevcorrSTTC":
        sttccorr = STTC.STTC()

    srate = (
        model_data.SI.dtIC * 1e-3
    )  # this needs to be adjusted by the date of the run, somewhere...
    # for runs prior to spring 2021, the 1e-3 is NOT needed.
    # for runs after that, the value is held in milliseconds, so needs to be
    # converted to seconds

    # sum across trials, and sort by inputs
    #
    # print("starting loop")
    RCD.sv_trials = []
    RCD.npost_spikes = 0
    for trial in range(RCP.ntrials):  # sum across trials
        # spiketimes is in msec, so leave si.dtIC in msec
        spikeindex = [int(t / (srate)) for t in data["Results"][trial]["spikeTimes"]]
        stx = AR.MC.time_base[spikeindex]
        stx = stx[  # get postsynaptic spikes and trim to analysis window
            (stx > RCP.min_time) & (stx < RCP.max_time)
        ]
        if len(stx) == 0:
            continue
        RCD.npost_spikes += len(stx)
        post_intervals.extend(np.diff(stx))
        # accumulate spikes and calculate average spike
        for n in range(len(stx)):  # for all the spikes that were detected
            reltime = np.around(RCD.ti * 1e3, 5) - np.around(stx[n], 5)
            areltime = np.argwhere(
                (RCP.minwin <= reltime) & (reltime <= RCP.maxwin)
            ).squeeze()

            if RCD.nsp_avg == 0:  # init arrays
                RCD.sv_avg = data["Results"][trial]["somaVoltage"][areltime]
                RCD.nsp_avg = 1
                RCD.ti_avg = 1e3 * RCD.ti[0 : len(areltime)] + RCP.minwin
            else:
                if len(areltime) > len(RCD.sv_avg):
                    areltime = areltime[0 : len(RCD.sv_avg)]
                if len(areltime) < len(RCD.sv_avg):
                    nextend = len(RCD.sv_avg) - len(areltime)
                    areltime = np.append(
                        areltime,
                        np.arange(areltime[-1] + 1, areltime[-1] + nextend + 1),
                    )
            if n == 0:
                RCD.sv_all = np.zeros(
                    (len(stx), RCD.ti_avg.shape[0])
                )  # initialize the array
                RCD.sv_avg += data["Results"][trial]["somaVoltage"][areltime]
                RCD.nsp_avg += 1
            RCD.sv_all[n, :] = data["Results"][trial]["somaVoltage"][areltime]

            nspk_plot += RCD.nsp_avg
        RCD.sv_trials.append(RCD.sv_all)
        max_spike_time =1e3*(data['runInfo'].pip_duration + data['runInfo'].pip_start)

        # Now get reverse  correlation for each input
        for isite in range(RCP.ninputs):  # for each ANF input
            # print(len( data["Results"][trial]["inputSpikeTimes"]), isite)
            anx = data["Results"][trial]["inputSpikeTimes"][
                isite
            ]  # get input AN spikes and trim list to window
            anx = anx[(anx > RCP.min_time) & (anx < RCP.max_time)]
            if len(anx) == 0:
                continue
            RCD.npre_spikes += len(anx)  # count up pre spikes.
            pre_intervals[isite].extend(
                np.diff(anx)
            )  # keep track of presynapit intervals by input
            if revcorrtype == "RevcorrSPKS":
                RCD.C[isite] += SPKS.correlogram(
                    stx * pq.ms,
                    anx * pq.ms,
                    width=-RCP.minwin * pq.ms,
                    bin_width=RCP.binw * pq.ms,
                    T=None,
                )
            elif revcorrtype == "RevcorrEleph":
                nbins = len(np.arange(RCP.minwin, -RCP.minwin, RCP.binw*2))
                bst_i = EC.BinnedSpikeTrain(
                    spiketrains=neo.SpikeTrain(stx*pq.ms, t_stop=max_spike_time*pq.s),
                    bin_size=RCP.binw*pq.ms,
                    n_bins=None,
                    t_start=data['runInfo'].pip_start* pq.s,
                    t_stop=(data['runInfo'].pip_start+data['runInfo'].pip_duration)* pq.s,
                    tolerance=None,
                    sparse_format="csr",
                )
                bst_j = EC.BinnedSpikeTrain(
                    spiketrains=neo.SpikeTrain(anx*pq.ms, t_stop=max_spike_time*pq.s),
                    bin_size=RCP.binw*pq.ms,
                    n_bins=None,
                    t_start=data['runInfo'].pip_start* pq.s,
                    t_stop=(data['runInfo'].pip_start+data['runInfo'].pip_duration)* pq.s,
                    tolerance=None,
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
                RCD.C[isite] += cc_result.squeeze()[:-1]
            elif revcorrtype == "RevcorrSimple":
                """
                The result returned from this version is not corrected
                for the binning or spike rate
                """
                refcorr, npost = reverse_correlation(
                    stx,
                    anx,
                    binwidth=RCP.binw,
                    # datawindow=[RCP.min_time, RCP.max_time], # stx and anx already have this done
                    corrwindow=[RCP.minwin, RCP.maxwin],
                )
                RCD.CB[isite] = RCD.CB[isite] + refcorr[:len(RCD.CB[isite])]

            elif revcorrtype == "RevcorrSTTC":
                sttccorr.set_spikes(RCP.binw, stx, anx, RCP.binw)
                refcorr = sttccorr.calc_ccf_sttc(
                    corrwindow=[RCP.minwin, RCP.maxwin], binwidth=RCP.binw
                )
                # print(len(refcorr))
                # print(len(RCD.STTC[isite]))
                # print(RCD.STTC[isite])
                RCD.STTC[isite] = RCD.STTC[isite] + refcorr

            tct = SPKS.total_correlation(anx, stx, width=-RCP.minwin, T=None)
            if ~np.isnan(tct):
                RCD.TC += tct

        if isinstance(RCD.TC, float) and RCD.TC > maxtc:
            maxtc = RCD.TC
        else:
            maxtc = 1.0

    if RCD.nsp_avg > 0:
        RCD.sv_avg /= RCD.nsp_avg
    RCD.mean_pre_intervals = [0] * RCP.ninputs
    RCD.pre_st = [[]] * RCP.ninputs
    for isite in range(RCP.ninputs):
        RCD.mean_pre_intervals[isite] = np.mean(pre_intervals[isite])
        RCD.pre_st[isite] = pre_intervals[isite]
    RCD.mean_post_intervals = np.mean(post_intervals)
    elapsed_time = datetime.datetime.now() - start_time
    print(f"    Time for calculation: {str(elapsed_time):s}")
    return model_data


def pairwise(model_data):
    allspikes = []
    data = model_data.data
    AR = model_data.AR
    RCP = model_data.RCP
    RCD = model_data.RCD
    RCD.pairwise = np.zeros((RCP.ninputs, RCP.ninputs))
    RCD.participation = np.zeros(RCP.ninputs)
    pre_spike_counts = np.zeros(
        RCP.ninputs + 1
    )  # there could be 0, or up to RCP.ninputs pre spikes
    pre_solo_spikes = np.zeros(RCP.ninputs + 1)
    srate = (
        model_data.SI.dtIC * 1e-3
    )  # this needs to be adjusted by the date of the run, somewhere...
    # for runs prior to spring 2021, the 1e-3 is NOT needed.
    # for runs after that, the value is held in milliseconds, so needs to be
    # converted to seconds
    nperspike = []
    nfilt_spikes = 0
    nfilt2_spikes = 0
    filttable = []
    filttable2 = []
    RCD.nspikes = 0
    sellist = [True] * RCP.ninputs
    # if ri.Spirou == "largestonly":
    #     for i in range(1, len(sellist)):
    #         sellist[i] = False
    # elif ri.Spirou == "twolargest":
    #     for i in range(2, len(sellist)):
    #         sellist[i] = False
    # elif ri.Spirou == "removelargest":
    #     sellist[0] = False
    print("RCP.ninputs: ", RCP.ninputs, RCD.sites)
    lack_largest = 0
    lack_two_largest = 0
    lack_three_largest = 0
    only_largest = 0
    only_two_largest = 0
    only_three_largest = 0

    for trial in range(RCP.ntrials):  # accumulate across all trials
        # spiketimes is in msec, so si.dtIC should be in msec
        spikeindex = [int(t / (srate)) for t in data["Results"][trial]["spikeTimes"]]
        spks = AR.MC.time_base[spikeindex]  # get postsynaptic spikes for the trial
        for n, s in enumerate(spks):  # for each postsynaptic spike
            if (
                s < RCP.min_time or s > RCP.max_time
            ):  # restrict post spikes to those only in a response window
                continue
            reltime = np.around(RCD.ti, 5) - np.around(s, 5)
            areltime = np.argwhere(
                (RCP.minwin <= reltime) & (reltime <= RCP.maxwin)
            ).squeeze()
            spikedata = SpikeData()  # store spike waveform using a dataclass
            spikedata.trial = trial
            spikedata.waveform = data["Results"][trial]["somaVoltage"][areltime]

            spikedata.time_index = n
            spikedata.prespikes = [[np.nan] for x in range(RCP.ninputs)]
            RCD.nspikes += 1  # number of post spikes evaluated
            n_active_inputs = (
                0  # number of active inputs associated with this post spike
            )
            spike_pattern = np.zeros(RCP.ninputs)
            solo = np.zeros(RCP.ninputs)
            # the 0'th site is the largest site.....
            for isite in range(RCP.ninputs):  # examine each input
                if not sellist[isite]:
                    continue
                npre_i, pre_times = _count_spikes_in_window(
                    data, trial, isite, s, RCP, RCD.pre_w
                )
                n_active_inputs += min(1, npre_i)  # only count once
                if npre_i > 0:

                    spikedata.prespikes[isite] = pre_times[
                        0
                    ]  # save all event times even if more than one
                    RCD.participation[
                        isite
                    ] += 1  # any spikes in the window = participation (but only count as 1)
                    spike_pattern[isite] += 1
                    # print(' spk: ', s, 'isite: ', isite)
                    for jsite in range(
                        isite + 1, RCP.ninputs
                    ):  # now do for joint combinations with other remaining inputs
                        npre_j, pre_times = _count_spikes_in_window(
                            data, trial, jsite, s, RCP, RCD.pre_w
                        )
                        # print(npre_j)
                        if npre_j > 0:  # accumulate if coincident for this pair
                            RCD.pairwise[isite, jsite] += 1

            # increment the number of times there were npre_spikes input to this post spike
            pre_spike_counts[n_active_inputs] += 1

            # check for different spike patterns (# of inputs, ordered)
            if np.sum(spike_pattern) == 1:  # only one input was active
                which_input = np.where(spike_pattern == 1)[0]
                pre_solo_spikes[which_input] += 1
            if spike_pattern[0] == 0 and np.sum(spike_pattern) >= 1:
                lack_largest += 1
            if np.sum(spike_pattern[0:2]) == 0 and np.sum(spike_pattern) >= 1:
                lack_two_largest += 1
            if np.sum(spike_pattern[0:3]) == 0 and np.sum(spike_pattern) >= 1:
                lack_three_largest += 1

            if np.sum(spike_pattern[0:1]) == 1 and np.sum(spike_pattern) == 1:
                only_largest += 1
            if np.sum(spike_pattern[0:2]) == 2 and np.sum(spike_pattern) == 2:
                only_two_largest += 1
            if np.sum(spike_pattern[0:3]) == 3 and np.sum(spike_pattern) == 3:
                only_three_largest += 1

            if sum(spike_pattern[0:5]) == 0:
                cprint(
                    "magenta",
                    f"{str(spike_pattern):s}, {int(np.sum(spike_pattern)):d}",
                )
                nfilt_spikes += 1
                filttable.append(spike_pattern)
            elif sum(spike_pattern[0:4]) == 0:
                cprint(
                    "cyan", f"{str(spike_pattern):s}, {int(np.sum(spike_pattern)):d}",
                )
                nfilt_spikes += 1
                filttable.append(spike_pattern)
            elif sum(spike_pattern[0:3]) == 0:
                cprint(
                    "blue", f"{str(spike_pattern):s}, {int(np.sum(spike_pattern)):d}",
                )
                nfilt_spikes += 1
                filttable.append(spike_pattern)
            elif sum(spike_pattern[0:2]) == 0:
                cprint(
                    "green", f"{str(spike_pattern):s}, {int(np.sum(spike_pattern)):d}",
                )
                nfilt_spikes += 1
                filttable.append(spike_pattern)
                if spike_pattern[2] == 1:
                    filttable2.append(spike_pattern)
                    nfilt2_spikes += 1
            # elif spike_pattern[0]== 0:
            #     cprint('yellow', f"{str(spike_pattern):s}, {int(np.sum(spike_pattern)):d}")
            # nfilt_spikes += 1
            # filttable.append(spike_pattern)
            else:
                pass
                # cprint("w", f"{str(spike_pattern):s}, {which_input[0]:d}")

            allspikes.append(spikedata)

    print("\nPairwise matrix: \n", RCD.pairwise)
    print(
        "Total prespikes: \n", sum(pre_spike_counts), [int(p) for p in pre_spike_counts]
    )
    print("\nTotal post spikes: ", RCD.nspikes)
    print("Windowd post spikes: ", RCD.npost_spikes)
    print("\nPre solo drive: ", pre_solo_spikes)
    print("\nFiltered Spikes: ", nfilt_spikes)
    print(
        f"Postspikes without the largest input active:      {lack_largest:5d} ({100.*lack_largest/RCD.nspikes:4.1f}%)"
    )
    print(
        f"Postspikes with only the largest input active:    {only_largest:5d} ({100.*only_largest/RCD.nspikes:4.1f}%)"
    )
    print(
        f"Postspikes without two largest inputs active:     {lack_two_largest:5d} ({100.*lack_two_largest/RCD.nspikes:4.1f}%)"
    )
    print(
        f"Postspikes with only two largest inputs active:   {only_two_largest:5d} ({100.*only_two_largest/RCD.nspikes:4.1f}%)"
    )
    print(
        f"Postspikes without three largest inputs active:   {lack_three_largest:5d} ({100.*lack_three_largest/RCD.nspikes:4.1f}%)"
    )
    print(
        f"Postspikes with only three largest inputs active: {only_three_largest:5d} ({100.*only_three_largest/RCD.nspikes:4.1f}%)"
    )
    print(
        f"Mean presynaptic rate: {np.mean([1./RCD.mean_pre_intervals[k] for k in range(len(RCD.mean_pre_intervals))]):f}"
    )
    print(f"Mean postsynaptic rate: {1./RCD.mean_post_intervals:f}")
    print(f"RCP min/max time: {RCP.min_time:8.3f} {RCP.max_time:8.3f}")

    # print('filttable: \n', filttable)
    filttable = np.array(filttable)
    filttable2 = np.array(filttable2)
    print("Counts: ", filttable.sum(axis=0))
    if filttable.shape[0] > 0:
        print("\nFilt Spike Proportions: ", filttable.sum(axis=0) / nfilt_spikes)
        fs1 = np.array(filttable)[:, 2:].sum(axis=0) / nfilt_spikes

    fsa = np.array(RCD.sites[2:]) / np.sum(RCD.sites[2:])
    print("Input Proportions: ", fsa)

    filttable2 = np.array(filttable2)
    # print('filttable 2 shape: ', filttable2.shape)
    if filttable2.shape[0] > 0:
        print(
            "\nFilt Spike Proportions on input #3: ",
            filttable2.sum(axis=0) / nfilt2_spikes,
        )
        fs2 = np.array(filttable2)[:, 2:].sum(axis=0) / nfilt2_spikes

    # if filttable.shape[0] > 0:
    #     f=mpl.figure()
    #     mpl.plot(fs1, fsa)
    #     mpl.show()
    # print('pre spike_count associated with a post spike: ', pre_spike_counts)
    # plot the position of the prespikes for every trial as determined by the
    # second trial loop above.
    # f, ax = mpl.subplots(1,1)
    # for i in range(len(self.allspikes)):
    #     y = i*np.ones(len(self.allspikes[i].prespikes))
    #     print(list(self.allspikes[i].prespikes))
    #     print(y)
    #     ax.plot(self.allspikes[i].prespikes, y)
    # mpl.show()
    # print(self.allspikes)
    npartipating = np.sum(RCD.participation)
    RCD.s_pair = np.sum(RCD.pairwise)
    if RCD.s_pair > 0.0:
        RCD.pairwise /= RCD.s_pair

    # print(np.unique(nperspike, return_counts=True))
    # nperspike = [n for n in nperspike if n != 0]
    # nperspike = scipy.stats.itemfreq(nperspike).T
    # print('nperspike counts: ', nperspike)
    # nperspike = np.array(np.unique(nperspike, return_counts=True))/nspikes
    # properly fill out output
    # xnspike = np.arange(RCP.ninputs)
    # ynspike = np.zeros(RCP.ninputs)
    # for j, i in enumerate(nperspike[0]):
    #     # print(i, j, nperspike[1,j])
    #     ynspike[i - 1] = nperspike[1, j]

    # ynspike = np.cumsum(ynspike / nspikes)
    cprint("r", f"prespikecounts: {np.sum(pre_spike_counts[1:]):f}")
    if np.sum(pre_spike_counts[1:]) == 0:
        return RCP, RCD, allspikes

    RCD.ynspike = np.cumsum(pre_spike_counts[1:]) / np.sum(pre_spike_counts[1:])
    return RCP, RCD, allspikes
