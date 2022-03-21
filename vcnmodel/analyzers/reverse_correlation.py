import datetime
from dataclasses import dataclass, field
from typing import List, Tuple, Union

import elephant.conversion as EC
import elephant.spike_train_correlation as ESTC
import neo
import numpy as np
import quantities as pq
from pylibrary.tools import cprint as CP
from vcnmodel.analyzers import spikestatistics as SPKS
from vcnmodel.analyzers import sttc as STTC

cprint = CP.cprint

"""
Compute the reverse correlation between two spike trains For all of the spikes
in st1, compute the time difference with each spike in st2, find where the
difference is within a time window we are interested in, then get and return the
indices for all of those spikes.

We define 3 data structures first: RevCorrPars are the parameters for the
calculations, and inlcude the stimulus and runinfo parameters from the
simulations. RevCorrData holds the results of the calculations. SpikeData holds
information about the individual postsynaptic spikes.

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
    minwin: float = (
        -5.0
    )  # start of revcorr display/calculation, ms relative to reference event
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


##########################################################################
# Reverse correlation analysis.
##########################################################################


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


def _remove_spikes(stx, anx, win=[-2.7, -0.5]):
    """Remove spikes from the input spike train (anx) that occur in the window
    specified relative to the postsynaptic spike train (stx)

    Parameters
    ----------
    stx : list
        postsynaptic spike times
    anx : list
        2D list of spike times, by input
    win : list, optional
        window prior to post spikes where, if the selected pre spike sources
        are present, the postspikes will be removed.
        default window = [-2.7, -0.5]
    """
    new_stx = []
    for i, x in enumerate(stx):
        ds = anx - x
        rmv = np.argwhere((ds >= win[0]) & (ds <= win[1]))
        if len(rmv) == 0:
            new_stx.append(x)
    #         if len(rmv) > 0 :
    #         anx[rmv] = np.nan
    # anx = anx[np.argwhere(~np.isnan(anx))]
    return new_stx


def revcorr(
    model_data: object = None,
    nbins: int = 0,
    revcorrtype: str = "",
    ASAs: list = [],
    revcorr_params: Union[dict, None] = None,
):
    """compute a reverse correlation, optinally with
    removal of selected postsynaptic spikes if certain inputs
    are active (e.g., > a minimum ASA)

    Parameters
    ----------
    model_data : object, optional
        the model_data returned from readmodel, default None
    nbins : int, optional
        number of bins in the revcorr by default 0
    revcorrtype : str, optional
        one of the possible types: "RevCorr SPKS" (from Brian1.4),
        "RevCorr Eleph" (from elephant),
        "RevCorr Simple" (dumb implementation), by default ""
    ASAs : list of floats, optional
        The ASAs in an ordered list, largest to smallest, for deselection, by default []
    revcorr_params : Union[dict, None], optional
        A few parameters: identifies deselection criterion by ASA threshold,
        whether deseleciton will be done, and the window used to count inputs
        and deselect them, by default None

    Returns
    -------
    _type_
        _description_
    """
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

    if revcorr_params is None:  # set some defaults
        revcorr_params = {
            "deselect": False,
            "threshold": 500.0,
            "window": [-2.7, -0.5],
        }

    if revcorr_params["deselect"]:
        deselect_sites = [
            (asa > revcorr_params["threshold"]) for i, asa in enumerate(ASAs)
        ]
    else:
        deselect_sites = [False] * len(ASAs)
    print("ASAs: ", ASAs)
    srate = (
        model_data.SI.dtIC * 1e-3
    )  # this needs to be adjusted by the date of the run, somewhere...
    # for runs prior to spring 2021, the 1e-3 is NOT needed.
    # for runs after that, the value is held in milliseconds, so needs to be
    # converted to seconds

    # sum across trials, and sort by inputs
    filtered_stx_by_trial = [
        []
    ] * RCP.ntrials  # postsynaptic spikes after filtering against inputs

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
        max_spike_time = 1e3 * (
            data["runInfo"].pip_duration + data["runInfo"].pip_start
        )
        an_spikes = [[]] * RCP.ninputs
        for isite in range(RCP.ninputs):  # for each ANF input
            anx = np.array(
                data["Results"][trial]["inputSpikeTimes"][isite]
            )  # get input AN spikes for this trial and input, and trim list to window
            anx = anx[np.argwhere((anx > RCP.min_time) & (anx < RCP.max_time))]
            an_spikes[isite] = anx
            # check if we want to examine postspikes without specific inputs
            if revcorr_params["deselect"] and deselect_sites[isite]:
                stx = _remove_spikes(stx, anx, revcorr_params["window"])
        filtered_stx_by_trial[trial] = stx  # save for later

        # Now get reverse  correlation for each input
        for isite in range(RCP.ninputs):  # for each ANF input
            # print(len( data["Results"][trial]["inputSpikeTimes"]), isite)
            anx = an_spikes[isite]
            # anx = anx[(anx > RCP.min_time) & (anx < RCP.max_time)]
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
                nbins = len(np.arange(RCP.minwin, -RCP.minwin, RCP.binw * 2))
                bst_i = EC.BinnedSpikeTrain(
                    spiketrains=neo.SpikeTrain(
                        stx * pq.ms, t_stop=max_spike_time * pq.s
                    ),
                    bin_size=RCP.binw * pq.ms,
                    n_bins=None,
                    t_start=data["runInfo"].pip_start * pq.s,
                    t_stop=(data["runInfo"].pip_start + data["runInfo"].pip_duration)
                    * pq.s,
                    tolerance=None,
                    sparse_format="csr",
                )
                bst_j = EC.BinnedSpikeTrain(
                    spiketrains=neo.SpikeTrain(
                        anx * pq.ms, t_stop=max_spike_time * pq.s
                    ),
                    bin_size=RCP.binw * pq.ms,
                    n_bins=None,
                    t_start=data["runInfo"].pip_start * pq.s,
                    t_stop=(data["runInfo"].pip_start + data["runInfo"].pip_duration)
                    * pq.s,
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
                RCD.CB[isite] = RCD.CB[isite] + refcorr[: len(RCD.CB[isite])]

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
    return model_data, filtered_stx_by_trial


##########################################################################
# Input pattern analysis
#
# In the following section, the analysis focusses on both
# pairwise interactions between inputs, and on the patterns
# of inputs. The questions revolve around how often the large,
# suprathreshold inputs drive a spike versus conincidence of
# smaller subthreshold inputs, under different conditions.
# We take advantage of the fact that we *know* which inputs
# are active prior to a spike.
#
##########################################################################


def _count_spikes_in_window(
    spike: float,
    data: dict,
    trial: int = 0,
    isite: int = 0,
    pre_w: list = [-2.7, -0.5],
    RCP: Union[object, None] = None,
):
    """Count presynaptic spikes in a window before a postsynaptic site

    Parameters
    ----------
    spike : float
        Time of the reference (postsynaptic) spike
    data : dict
        model data dictionary to access input spikes
    trial : int
        trial to look for input data
    isite : int
        input site number
    pre_w : list (len = 2)
        window before spike to look for input spikes
    RCP : object
        parameters for analysis - use to get time into the trace

    Returns
    -------
    npre_i : int
        number of presynaptic spikes that met criteria
    pre_times: a list of the presynaptic spikes that met criteria.
    """
    an_i = data["Results"][trial]["inputSpikeTimes"][
        isite
    ]  # get input spike times for one input in one trial
    an_i = an_i[
        (an_i > RCP.min_time) & (an_i < RCP.max_time)
    ]  # restrict to those only within the response window
    an_i = an_i - spike  # get relative latency from spike to it's inputs
    pre_indx = np.asarray((an_i >= pre_w[0]) & (an_i <= pre_w[1])).nonzero()[0]
    pre_times = [an_i[k] for k in pre_indx]
    if len(pre_times) > 0:
        npre_i = len(pre_times)
    else:
        npre_i = 0
        pre_times = []
    return npre_i, pre_times


def _get_spike_pattern(
    spike: float,
    data: dict,
    trial: int = 0,
    selected_inputs: list = [],
    RCP: object = None,
    RCD: object = None,
):
    """Compute the pattern of input spikes associated with an individual
     postsynaptic spike.

     Return the pattern (list, by input)

     Parameters
     ----------
     data : _type_
         _description_
     trial : _type_
         _description_
     selected_inputs : _type_
         _description_
     spikedata : _type_
         _description_
     RCP : _type_
         _description_
     RCD : _type_
         _description_

     Returns
     -------
     spike_pattern: A list by trial of lists by input of inputs
         that were associated with a postsynaptic spike
    pair_wise: A list by trial of 2D matrices indicating the
         paired occurance of inputs for a given postsynaptic spike
     total_pre: int. Total count of presynaptic spikes tested.
    """
    n_active_inputs = 0
    total_pre = 0
    pair_wise = np.zeros((RCP.ninputs, RCP.ninputs))
    spike_pattern = np.zeros(RCP.ninputs)
    # the 0'th site is the largest site.....
    for isite, site_flag in enumerate(selected_inputs):  # examine each input
        npre_i, pre_times = _count_spikes_in_window(
            spike,
            data=data,
            trial=trial,
            isite=isite,
            pre_w=RCD.pre_w,
            RCP=RCP,
        )
        total_pre += npre_i
        n_active_inputs += min(1, npre_i)  # only count once
        if npre_i > 0:
            spike_pattern[isite] += 1
            for jsite in range(
                isite + 1, RCP.ninputs
            ):  # now do for joint combinations with other remaining inputs
                npre_j, pre_times = _count_spikes_in_window(
                    spike,
                    data=data,
                    trial=trial,
                    isite=jsite,
                    pre_w=RCD.pre_w,
                    RCP=RCP,
                )
                if npre_j > 0:  # accumulate if coincident for this pair
                    pair_wise[isite, jsite] += 1
    return spike_pattern, pair_wise, total_pre


def _spike_filter(RCP, spks):
    # filter to analysis window.
    spks = spks[(spks >= RCP.min_time) & (spks < RCP.max_time)]
    return spks


@dataclass
class PatternData:
    """This dataclass holds
    information about the spike patterns analyzed
    The mask is a bit mask for the patterns
    """

    name: str = ""
    mask: int = 0
    logic: str = 'exact' # set to 'exact' to require mask match
                         # 'except' for selected inputs to be *inactive*
                         # 'atleast' for for selected inputs to be active (but others are ok too)
    sumtrue: int = 0  # given the logic, this counts those that match
    sumfalse: int = 0  # and those that do not match


def _print_pattern(name, pat):
    mask = f"{pat.mask:016b}"[::-1]
    print(
        f"    Pattern name: {name:18s} mask: {mask:16s}  Logic: {str(pat.logic):5s}  True: {pat.sumtrue:5d}  False: {pat.sumfalse:5d}"
    )


def _spike_pattern_analysis(
    spikes: Union[list, np.ndarray], spike_pattern: list, ninputs: int
):
    """_summary_

    Parameters
    ----------
    spikes : Union[list, np.ndarray]
        the list of all spike times (do we need actually this?)
    spike_pattern : list
        The patterns of presynaptic spikes for all of the spikes
    ninputs : int
        The number of inputs in the spike_pattern

    Returns
    -------
    patterns : dict
        The patterns dictionary that was used, with the values
        filled in in the PatternData dataclass
    solo spikes: list
        A summary of which input gave rise on it's own to a postsynaptic spike

    """
    pre_spike_counts = np.zeros(ninputs + 1)
    # there could be 0, or up to RCP.ninputs pre spikes
    pre_solo_spikes = np.zeros(ninputs + 1)
    # patterns are [bitmask(largest input in lowest/LSB position), logic, counts lack input, count have_input]
    patterns = {  # some patterns to match conditions
        "1_largest": PatternData(mask=0x01, logic='exact'),
        "2_largest": PatternData(mask=0x03, logic='exact'),
        "3_largest": PatternData(mask=0x07, logic='exact'),
        "2nd_largest": PatternData(mask=0x02, logic='exact'),
        "3rd_largest": PatternData(mask=0x04, logic='exact'),
        "4th_largest": PatternData(mask=0x08, logic='exact'),
        "1+2+5": PatternData(mask=0b00010011, logic='exact'),
        "4+5": PatternData(mask=  0b00011000, logic='exact'),
        "4+5+6": PatternData(mask=0b00111000, logic='exact'),
        "6+7+8": PatternData(mask=0b11100000, logic='exact'),
        "None": PatternData(mask=0x00, logic='exact'),
        "not_1_largest": PatternData(mask=0x01, logic='except'), # not the largest
        "not_2_largest": PatternData(mask=0x03, logic='except'), # none of the 2 largest
        "not_3_largest": PatternData(mask=0x07, logic='except'), # none of the 3 largest
        "at_least": PatternData(mask=0x55, logic='atleast'),  # at least some in a pattern
        "at_least2": PatternData(mask=0x07, logic='atleast'),  # at least the 3 largest
    }
    for p in patterns.keys():
        patterns[p].name = p
    b_spk = 0
    for i, spk in enumerate(spike_pattern):
        if spk == 1:
            b_spk |= 2**i
    for n, spike in enumerate(spikes):  # for each postsynaptic spike
        n_active_inputs = 0  # number of active inputs associated with this post spike
        # increment the number of times there were npre_spikes for this post spike
        pre_spike_counts[n_active_inputs] += 1
        #
        # # check for different spike patterns (# of inputs, ordered)
        # Test if one and only one input was active (bspk is a power of 2)
        if (b_spk > 0) & (b_spk & (b_spk - 1) == 0):
            which_input = int(np.log2(b_spk))
            pre_solo_spikes[which_input] += 1
        # now check the specific input patterns.
        for patname, pat in patterns.items():
            if pat.logic == 'exact':  # check case in which the selected inputs are active
                if b_spk == pat.mask:  # has exact pattern of input(s) only
                    pat.sumtrue += 1
                else:
                    pat.sumfalse += 1  # other inputs perhaps... 
            elif pat.logic == 'except':  # evaluate the case in which the zeroed inputs are NOT active
                if b_spk & pat.mask:  # any active input is ok, except the 0'd ones.
                    pat.sumfalse += 1
                else:
                    pat.sumtrue += 1
            elif pat.logic == 'atleast':
                combo = b_spk & pat.mask
                if combo >= pat.mask:  # any combination will 
                    pat.sumtrue += 1
                else:
                    pat.sumfalse += 1
            else:
                raise ValueError(f"Logic must be one of: 'exact', 'except' or 'atleast', got: {pat.logic:s}")
    
    return patterns, pre_solo_spikes

def _assert_patterns(sp, spike_patterns, r):
    n_fails = 0
    n_successes = 0
    for pattern in spike_patterns:
        if pattern == sp:
            if r[pattern].logic == 'exact':
                try:
                    assert r[pattern].sumtrue == 1 and r[pattern].sumfalse == 0
                    n_successes += 1
                except AssertionError:
                    cprint('r', "    Assert failed for selected pattern with exact")
                    cprint('r', f"       sp: {sp:s}   vs {pattern:s}")
                    _print_pattern(pattern, r[pattern])
                    print("      ", r[pattern])
                    n_fails += 1
            elif r[pattern].logic == 'except':
                try:
                    assert r[pattern].sumfalse == 1 and r[pattern].sumtrue == 0
                    n_successes += 1
                except AssertionError:
                    cprint('c', "    Assert failed for non-selected pattern with except")
                    cprint('c', f"       sp: {sp:s}   vs {pattern:s}  with logic: {r[pattern].logic:s}")
                    _print_pattern(pattern, r[pattern])
                    n_fails += 1
            elif r[pattern].logic == 'atleast':
                try:
                    assert r[pattern].sumtrue == 1 and r[pattern].sumfalse == 0
                    n_successes += 1
                except AssertionError:
                    cprint('m', "    Assert failed for non-selected pattern with except")
                    cprint('m', f"       sp: {sp:s}   vs {pattern:s}  with logic: {r[pattern].logic:s}")
                    _print_pattern(pattern, r[pattern])
                    n_fails += 1
    if n_fails > 0:
    #     cprint("g", f"     Asserts passed for {sp:s} with {n_successes:d}")
    # else:
        cprint("r", f"     {n_fails:d} Asserts failed for {sp:s}")
    return n_fails


def test_spike_patterns():
    spikes = [64] * 1
    spike_patterns = {  # some test spike patterns
        "1_largest": [1, 0, 0, 0, 0, 0, 0, 0],  # largest only
        "2_largest": [1, 1, 0, 0, 0, 0, 0, 0],  # 2 largest
        "3_largest": [1, 1, 1, 0, 0, 0, 0, 0],  # 3largest
        "2nd_largest": [0, 1, 0, 0, 0, 0, 0, 0],  # second only
        "3rd_largest": [0, 0, 1, 0, 0, 0, 0, 0],  # 3rd only
        "4th_largest": [0, 0, 0, 1, 0, 0, 0, 0],  # 4th only
        "1+2+5": [1, 1, 0, 0, 1, 0, 0, 0],  # 2 largest and a smaller one
        "4+5": [0, 0, 0, 1, 1, 0, 0, 0],  # 4th and 5th  Should pass 0 largest
        "4+5+6": [0, 0, 0, 1, 1, 1, 0, 0],
        "6+7+8": [0, 0, 0, 0, 0, 1, 1, 1],  # bunch of little ones - pass 0 largest
        "None": [0, 0, 0, 0, 0, 0, 0, 0],
        "not_1_largest": [1, 0, 0, 0, 0, 0, 0, 0],  # to test the ones we are against... 
        "not_2_largest": [1, 1, 0, 0, 0, 0, 0, 0],
        "not_3_largest": [1, 1, 1, 0, 0, 0, 0, 0],
        "at_least": [1, 0, 1, 0, 1, 0, 1, 0],
        "at_least2": [1, 1, 1, 0, 1, 0, 1, 0],
    }
    ninputs = 8
    for sp in spike_patterns:
        r, solo = _spike_pattern_analysis(spikes, spike_patterns[sp], ninputs)
        _assert_patterns(sp, spike_patterns, r)

        print("solo: ", solo)

    return

    # if spike_pattern[0] == 0 and np.sum(spike_pattern) >= 1:
    #     lack_largest += 1
    # if np.sum(spike_pattern[0:2]) == 0 and np.sum(spike_pattern) >= 1:
    #     lack_two_largest += 1
    # if np.sum(spike_pattern[0:3]) == 0 and np.sum(spike_pattern) >= 1:
    #     lack_three_largest += 1

    # if np.sum(spike_pattern[0:1]) == 1 and np.sum(spike_pattern) == 1:
    #     only_largest += 1
    # if np.sum(spike_pattern[0:2]) == 2 and np.sum(spike_pattern) == 2:
    #     only_two_largest += 1
    # if np.sum(spike_pattern[0:3]) == 3 and np.sum(spike_pattern) == 3:
    #     only_three_largest += 1
    # print_patterns = False
    # if sum(spike_pattern[0:5]) == 0:
    #     if print_patterns:
    #         cprint(
    #         "magenta",
    #         f"{str(spike_pattern):s}, {int(np.sum(spike_pattern)):d}",
    #     )
    #     filter_table1.append(spike_pattern)
    # elif sum(spike_pattern[0:4]) == 0:
    #     if print_patterns:
    #         cprint(
    #         "cyan", f"{str(spike_pattern):s}, {int(np.sum(spike_pattern)):d}",
    #         )

    #     filter_table2.append(spike_pattern)
    # elif sum(spike_pattern[0:3]) == 0:
    #     if print_patterns:
    #         cprint(
    #         "blue", f"{str(spike_pattern):s}, {int(np.sum(spike_pattern)):d}",
    #     )
    #     filter_table1.append(spike_pattern)
    # elif sum(spike_pattern[0:2]) == 0:
    #     if print_patterns:
    #         cprint(
    #         "green", f"{str(spike_pattern):s}, {int(np.sum(spike_pattern)):d}",
    #     )
    #     filter_table2.append(spike_pattern)
    # if spike_pattern[2] == 1:
    #         filter_table2.append(spike_pattern)
    # # elif spike_pattern[0]== 0:
    # #     cprint('yellow', f"{str(spike_pattern):s}, {int(np.sum(spike_pattern)):d}")
    # # nfilt_spikes += 1
    # # filttable.append(spike_pattern)
    # else:
    #     pass
    #     # cprint("w", f"{str(spike_pattern):s}, {which_input[0]:d}")

    # return patterns  # , filter_table1, filter_table2

    # def print_pairwise(RCD, RCP, PAT):
    #     print("\nPairwise matrix: \n", RCD.pairwise)
    #     print(
    #         "Total prespikes: \n", sum(pre_spike_counts), [int(p) for p in pre_spike_counts]
    #     )
    #     print("\nTotal post spikes: ", RCD.nspikes)
    #     print("Windowd post spikes: ", RCD.npost_spikes)
    #     print("\nPre solo drive: ", PAT.pre_solo_spikes)
    #     print("\nFiltered Spikes: ", nfilt_spikes)
    #     print(
    #         f"Postspikes without the largest input active:      {lack_largest:5d} ({100.*lack_largest/RCD.nspikes:4.1f}%)"
    #     )
    #     print(
    #         f"Postspikes with only the largest input active:    {only_largest:5d} ({100.*only_largest/RCD.nspikes:4.1f}%)"
    #     )
    #     print(
    #         f"Postspikes without two largest inputs active:     {lack_two_largest:5d} ({100.*lack_two_largest/RCD.nspikes:4.1f}%)"
    #     )
    #     print(
    #         f"Postspikes with only two largest inputs active:   {only_two_largest:5d} ({100.*only_two_largest/RCD.nspikes:4.1f}%)"
    #     )
    #     print(
    #         f"Postspikes without three largest inputs active:   {lack_three_largest:5d} ({100.*lack_three_largest/RCD.nspikes:4.1f}%)"
    #     )
    #     print(
    #         f"Postspikes with only three largest inputs active: {only_three_largest:5d} ({100.*only_three_largest/RCD.nspikes:4.1f}%)"
    #     )
    #     print(
    #         f"Mean presynaptic rate: {np.mean([1./RCD.mean_pre_intervals[k] for k in range(len(RCD.mean_pre_intervals))]):f}"
    #     )
    #     print(f"Mean postsynaptic rate: {1./RCD.mean_post_intervals:f}")
    #     print(f"RCP min/max time: {RCP.min_time:8.3f} {RCP.max_time:8.3f}")

    # print('filttable: \n', filttable)
    # filttable = np.array(filttable)
    # filttable2 = np.array(filttable2)
    # print("Counts: ", filttable.sum(axis=0))
    # if filttable.shape[0] > 0:
    #     print("\nFilt Spike Proportions: ", filttable.sum(axis=0) / nfilt_spikes)
    #     fs1 = np.array(filttable)[:, 2:].sum(axis=0) / nfilt_spikes

    # fsa = np.array(RCD.sites[2:]) / np.sum(RCD.sites[2:])
    # print("Input Proportions: ", fsa)

    # filttable2 = np.array(filttable2)
    # # print('filttable 2 shape: ', filttable2.shape)
    # if filttable2.shape[0] > 0:
    #     print(
    #         "\nFilt Spike Proportions on input #3: ",
    #         filttable2.sum(axis=0) / nfilt2_spikes,
    #     )
    #     fs2 = np.array(filttable2)[:, 2:].sum(axis=0) / nfilt2_spikes

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


@dataclass()
class Patterns:
    """
    Hold results for different spike patterns
    """

    selected_inputs: field(default_factory=def_empty_list)
    event_spike_pattern: field(default_factory=def_empty_list)
    event_pair_wise: field(default_factory=def_empty_list)
    filter_table: field(default_factory=def_empty_list)
    filter_table2: field(default_factory=def_empty_list)
    all_spikes: field(default_factory=def_empty_list)
    n_post: int = 0
    n_pre: int = 0
    lack_largest: int = 0
    lack_two_largest: int = 0
    lack_three_largest: int = 0
    only_largest: int = 0
    only_two_largest: int = 0
    only_three_largest: int = 0


def spike_pattern_analysis(model_data):
    """Analyze input spike patterns prior to postsynaptic spike

    Parameters
    ----------
    model_data : object
        The model_data class returned from readmodel

    Returns
    -------
    RCP: reverse correlation parameter data structure
    RCD: reverse correlation result data structure
    PAT.all_spikes: pattern data structure for all spikes
    """
    data = model_data.data
    AR = model_data.AR
    RCP = model_data.RCP
    RCD = model_data.RCD
    RCD.pairwise = np.zeros((RCP.ninputs, RCP.ninputs))
    RCD.participation = np.zeros(RCP.ninputs)

    srate = (
        model_data.SI.dtIC * 1e-3
    )  # this needs to be adjusted by the date of the run, somewhere...
    # for runs prior to spring 2021, the 1e-3 is NOT needed.
    # for runs after that, the value is held in milliseconds, so needs to be
    # converted to seconds
    PAT = Patterns
    PAT.selected_inputs = [True] * RCP.ninputs  # all inputs
    PAT.event_spike_pattern = []
    PAT.event_pair_wise = []
    PAT.n_post = np.zeros(RCP.ntrials)
    PAT.n_pre = np.zeros(RCP.ntrials)
    PAT.filttable = []
    PAT.filttable2 = []
    PAT.all_spikes = []
    for trial in range(RCP.ntrials):  # accumulate across all trials

        spikeindex = [int(t / (srate)) for t in data["Results"][trial]["spikeTimes"]]
        spikes = AR.MC.time_base[spikeindex]  # get postsynaptic spikes for the trial
        spikes = _spike_filter(RCP, spikes)
        PAT.n_post[trial] = len(spikes)
        for _, spike in enumerate(spikes):  # for each spike, get the input pattern
            spike_pattern, pair_wise, total_pre = _get_spike_pattern(
                spike=spike,
                data=data,
                trial=trial,
                selected_inputs=PAT.selected_inputs,
                RCP=RCP,
                RCD=RCD,
            )
            PAT.n_pre[trial] = total_pre
            PAT.event_spike_pattern.append(spike_pattern)
            PAT.event_pair_wise.append(pair_wise)
            PAT.all_spikes.append(spike)  # match up the spike with the pattern
            patterns = _spike_pattern_analysis(spikes, spike_pattern, RCP.ninputs)
            PAT.filttable.append(patterns)
    print(
        "sum of spikes: ",
        np.sum(PAT.event_spike_pattern, axis=0),
        np.sum(PAT.event_spike_pattern),
    )
    RCD.participation = np.sum(PAT.event_spike_pattern, axis=0) / np.sum(
        PAT.event_spike_pattern
    )
    print(f"Participation: {str(RCD.participation):s}")
    nparticipating = np.sum(RCD.participation)
    print("n participating: ", nparticipating)
    print("PAT.n_post and PAT.n_pre: ", PAT.n_post, PAT.n_pre)
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
    cprint("r", f"prespikecounts: {np.sum(PAT.n_pre[1:]):f}")
    if np.sum(PAT.n_pre[1:]) == 0:
        return RCP, RCD, PAT.allspikes

    RCD.ynspike = np.cumsum(np.sum(PAT.n_pre[1:]) / np.sum(PAT.n_pre[1:]))
    return RCP, RCD, PAT.all_spikes


# def pairwise_old(model_data):
#     allspikes = []
#     data = model_data.data
#     AR = model_data.AR
#     RCP = model_data.RCP
#     RCD = model_data.RCD
#     RCD.pairwise = np.zeros((RCP.ninputs, RCP.ninputs))
#     RCD.participation = np.zeros(RCP.ninputs)
#     pre_spike_counts = np.zeros(
#         RCP.ninputs + 1
#     )  # there could be 0, or up to RCP.ninputs pre spikes
#     pre_solo_spikes = np.zeros(RCP.ninputs + 1)
#     srate = (
#         model_data.SI.dtIC * 1e-3
#     )  # this needs to be adjusted by the date of the run, somewhere...
#     # for runs prior to spring 2021, the 1e-3 is NOT needed.
#     # for runs after that, the value is held in milliseconds, so needs to be
#     # converted to seconds
#     nperspike = []
#     nfilt_spikes = 0
#     nfilt2_spikes = 0
#     filttable = []
#     filttable2 = []
#     RCD.nspikes = 0
#     sellist = [True] * RCP.ninputs

#     print(f"RCP.ninputs {RCP.ninputs:d}", "sites: ", RCD.sites)
#     for trial in range(RCP.ntrials):  # accumulate across all trials
#         # spiketimes is in msec, so si.dtIC should be in msec
#         spikeindex = [int(t / (srate)) for t in data["Results"][trial]["spikeTimes"]]
#         spks = AR.MC.time_base[spikeindex]  # get postsynaptic spikes for the trial
#         allspikes = _spike_filter(data, spks)

#     print("\nPairwise matrix: \n", RCD.pairwise)
#     print(
#         "Total prespikes: \n", sum(pre_spike_counts), [int(p) for p in pre_spike_counts]
#     )
#     print("\nTotal post spikes: ", RCD.nspikes)
#     print("Windowd post spikes: ", RCD.npost_spikes)
#     print("\nPre solo drive: ", pre_solo_spikes)
#     print("\nFiltered Spikes: ", nfilt_spikes)
#     print(
#         f"Postspikes without the largest input active:      {lack_largest:5d} ({100.*lack_largest/RCD.nspikes:4.1f}%)"
#     )
#     print(
#         f"Postspikes with only the largest input active:    {only_largest:5d} ({100.*only_largest/RCD.nspikes:4.1f}%)"
#     )
#     print(
#         f"Postspikes without two largest inputs active:     {lack_two_largest:5d} ({100.*lack_two_largest/RCD.nspikes:4.1f}%)"
#     )
#     print(
#         f"Postspikes with only two largest inputs active:   {only_two_largest:5d} ({100.*only_two_largest/RCD.nspikes:4.1f}%)"
#     )
#     print(
#         f"Postspikes without three largest inputs active:   {lack_three_largest:5d} ({100.*lack_three_largest/RCD.nspikes:4.1f}%)"
#     )
#     print(
#         f"Postspikes with only three largest inputs active: {only_three_largest:5d} ({100.*only_three_largest/RCD.nspikes:4.1f}%)"
#     )
#     print(
#         f"Mean presynaptic rate: {np.mean([1./RCD.mean_pre_intervals[k] for k in range(len(RCD.mean_pre_intervals))]):f}"
#     )
#     print(f"Mean postsynaptic rate: {1./RCD.mean_post_intervals:f}")
#     print(f"RCP min/max time: {RCP.min_time:8.3f} {RCP.max_time:8.3f}")

#     # print('filttable: \n', filttable)
#     filttable = np.array(filttable)
#     filttable2 = np.array(filttable2)
#     print("Counts: ", filttable.sum(axis=0))
#     if filttable.shape[0] > 0:
#         print("\nFilt Spike Proportions: ", filttable.sum(axis=0) / nfilt_spikes)
#         fs1 = np.array(filttable)[:, 2:].sum(axis=0) / nfilt_spikes

#     fsa = np.array(RCD.sites[2:]) / np.sum(RCD.sites[2:])
#     print("Input Proportions: ", fsa)

#     filttable2 = np.array(filttable2)
#     # print('filttable 2 shape: ', filttable2.shape)
#     if filttable2.shape[0] > 0:
#         print(
#             "\nFilt Spike Proportions on input #3: ",
#             filttable2.sum(axis=0) / nfilt2_spikes,
#         )
#         fs2 = np.array(filttable2)[:, 2:].sum(axis=0) / nfilt2_spikes

#     # if filttable.shape[0] > 0:
#     #     f=mpl.figure()
#     #     mpl.plot(fs1, fsa)
#     #     mpl.show()
#     # print('pre spike_count associated with a post spike: ', pre_spike_counts)
#     # plot the position of the prespikes for every trial as determined by the
#     # second trial loop above.
#     # f, ax = mpl.subplots(1,1)
#     # for i in range(len(self.allspikes)):
#     #     y = i*np.ones(len(self.allspikes[i].prespikes))
#     #     print(list(self.allspikes[i].prespikes))
#     #     print(y)
#     #     ax.plot(self.allspikes[i].prespikes, y)
#     # mpl.show()
#     # print(self.allspikes)
#     npartipating = np.sum(RCD.participation)
#     RCD.s_pair = np.sum(RCD.pairwise)
#     if RCD.s_pair > 0.0:
#         RCD.pairwise /= RCD.s_pair

#     # print(np.unique(nperspike, return_counts=True))
#     # nperspike = [n for n in nperspike if n != 0]
#     # nperspike = scipy.stats.itemfreq(nperspike).T
#     # print('nperspike counts: ', nperspike)
#     # nperspike = np.array(np.unique(nperspike, return_counts=True))/nspikes
#     # properly fill out output
#     # xnspike = np.arange(RCP.ninputs)
#     # ynspike = np.zeros(RCP.ninputs)
#     # for j, i in enumerate(nperspike[0]):
#     #     # print(i, j, nperspike[1,j])
#     #     ynspike[i - 1] = nperspike[1, j]

#     # ynspike = np.cumsum(ynspike / nspikes)
#     cprint("r", f"prespikecounts: {np.sum(pre_spike_counts[1:]):f}")
#     if np.sum(pre_spike_counts[1:]) == 0:
#         return RCP, RCD, allspikes

#     RCD.ynspike = np.cumsum(pre_spike_counts[1:]) / np.sum(pre_spike_counts[1:])
#     return RCP, RCD, allspikes


if __name__ == "__main__":
    test_spike_patterns()
