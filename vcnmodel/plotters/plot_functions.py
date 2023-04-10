"""
"""
from typing import Union, Tuple, List

import matplotlib.pyplot as mpl
import numpy as np
import seaborn as sns
from vcnmodel.util import trace_calls
from vcnmodel.analyzers import analysis as SPKANA
import vcnmodel.analyzers.flatten_spike_array as VAFlatten
from vcnmodel.util.basic_units import radians
from vcnmodel.util import make_sound_waveforms

TRC = trace_calls.TraceCalls


def twinax(fig: object, ax1: object, pos: float = 0.0) -> object:
    """
    Create a 'twin' axis on the right side of a plot
    Note: pyqtgraph.plotting.styles can also do an inset
    which may be used instead
    """
    ax2 = fig.add_axes(ax1.get_position(True), sharex=ax1)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.set_offset_position("right")
    ax2.tick_params(direction="in", length=5.0, width=1.0, labelsize=6)
    ax2.spines["right"].set_position(("data", pos))
    ax2.spines["left"].set_color("none")
    ax2.spines["top"].set_color("none")
    # ax2.set_autoscalex_on(ax1.get_autoscalex_on())
    # ax1.yaxis.tick_left()
    ax2.xaxis.set_visible(False)
    ax2.patch.set_visible(False)
    #    PH.adjust_spines(ax2, distance=0.)
    return ax2


def add_stimbar(stimbar, zero_time, max_time, ax):
    if stimbar['sound_params'] is None and stimbar['waveform'] is None:
        return
    if stimbar['sound_params'] is not None:
        if stimbar['sound_params'].make_sound():
            wv = stimbar['sound_params'].stimWaveform
            wtb = stimbar['sound_params'].stimTimebase
    if stimbar['waveform'] is not None:
        wv = stimbar['waveform'][1]
        wtb = stimbar['waveform'][0]
    twin = np.argwhere((wtb > zero_time) & (wtb < max_time))
    wv = wv/np.max(wv)
    # now build a second axis below the first

    ax.plot(wtb[twin]-zero_time, wv[twin]-2, 'b-', clip_on=False)  # plot outside the axes
    ax.xaxis.set_tick_params(which='major', pad=50)



@TRC()
def plot_psth(
    spike_times: Union[list, np.ndarray],
    run_info: object,
    zero_time: float = 0.0,
    max_time: float = 1.0,
    bin_width: float = 1e-3,
    ax: Union[object, None] = None,
    scale: float = 1.0,
    bin_fill: bool = True,
    edge_color: str = "k",
    alpha: float = 1.0,
    xunits: str = "time",
    stimbar: Union[dict, None] = None, # {'sound_params': None, 'waveform': None},
):
    """Correctly plot PSTH with spike rate in spikes/second
    All values are in seconds (times, binwidths)

    Parameters
    ----------
    spike_times : Union[list, np.ndarray]
        A list of the spike times (in seconds)
    run_info : object
        The model_params run_info dataclass associated with this data
    zero_time : float, optional
        The onset time for the stimulus for reference, by default 0.0
    max_time : float, optional
        Maximum time for the FSL/SSL plot, by default 1.0
    bin_width : float, optional
        Width of bins for FSL/SSL plot, by default 1e-3
    ax : Union[object, None], optional
        matplotlib axis instance for this plotm by default None
        note: if None, no plot is generated
    scale : float, optional
        Weights for the histogram optional. Used to scale (for example,
        by the number of cells or fibers)
    bin_fill : bool, optional
        Fill the bins, or leave open (transparent), by default True
    edge_color : str, optional
        color of the edges of the plot bars, by default "k"
    alpha : float, optional
        Transparency of the plot, by default 1.0
    xunits: str, optional
        The text for the x axis units
    stimbar: dict{sound_params: None, waveform: None}
        A dictionary for placing a plot of the stimulus under the plot.
        Keys: 'sound_params', 'waveform'
        if 'sound_params' is not None, then it should be a sound_params dataclass
        from make_sound_waveforms, and will be used to generate the waveform.
        if 'waveform' is not None, then it is assumed that the passed
        argument is a list [time, waveform].

    Returns
    -------
    Nothing
    """
    num_trials = run_info.nReps
    spike_times_flat = VAFlatten.flatten_spike_array(spike_times,
                        time_window=(zero_time, max_time), isi_flag=False)
    print("Spike times flat len: ", len(spike_times_flat), np.min(spike_times_flat), np.max(spike_times_flat))
    spike_times_flat -= zero_time
    bins = np.arange(0.0, max_time - zero_time, bin_width)
    if bin_fill:
        face_color = edge_color
    else:
        face_color = None  # "None"
    if xunits == "radians":
        xu = radians
        bins = bins * radians
    else:
        xu = None
    
    if (
        (len(spike_times_flat) > 0)
        # and not np.isnan(np.sum(spike_times_flat)))
        and (ax is not None)
    ):
        print("facecoclor, edgecolor: ", face_color, edge_color)
        ax.hist(
            x=spike_times_flat,
            bins=bins,
            density=False,
            # weights=(1./(bin_width*num_trials)) * np.ones_like(spike_times_flat),
            histtype="stepfilled",
            facecolor=face_color,
            edgecolor=edge_color,
            linewidth=None,
            rwidth=1.0,
            alpha=alpha, # alpha,
        )
    else:
        if ax is not None:
            ax.text(
                0.5,
                0.5,
                "No Spikes",
                fontsize=14,
                color="r",
                transform=ax.transAxes,
                horizontalalignment="center",
            )
    if xunits == "radians":
        ticks = np.linspace(0, 2 * np.pi, 5)
        tick_labels = ["0", r"$\pi /2$", r"$\pi$", r"$3\pi /2$", r"$2\pi$"]
        ax.set_xticks(ticks, tick_labels)
    
    # ax.set_ylim(0, ax.get_ylim()[1])
    if stimbar is not None:
        add_stimbar(stimbar, zero_time, max_time, ax)

    return

def plot_isi(
    spike_times: Union[list, np.ndarray],
    run_info: object,
    zero_time: float = 0.0,
    max_time: float = 1.0,
    bin_width: float = 1e-3,
    ax: Union[object, None] = None,
    scale: float = 1.0,
    bin_fill: bool = True,
    edge_color: str = "k",
    alpha: float = 1.0,
    xunits: str = "time",
    stimbar: Union[dict, None] = None, # {'sound_params': None, 'waveform': None},
):
    """Correctly plot ISI distriution
    All values are in seconds (times, binwidths)

    Parameters
    ----------
    spike_times : Union[list, np.ndarray]
        A list of the spike times (in seconds)
    run_info : object
        The model_params run_info dataclass associated with this data
    zero_time : float, optional
        The onset time for the stimulus for reference, by default 0.0
    max_time : float, optional
        Maximum time for the FSL/SSL plot, by default 1.0
    bin_width : float, optional
        Width of bins for FSL/SSL plot, by default 1e-3
    ax : Union[object, None], optional
        matplotlib axis instance for this plotm by default None
        note: if None, no plot is generated
    scale : float, optional
        Weights for the histogram optional. Used to scale (for example,
        by the number of cells or fibers)
    bin_fill : bool, optional
        Fill the bins, or leave open (transparent), by default True
    edge_color : str, optional
        color of the edges of the plot bars, by default "k"
    alpha : float, optional
        Transparency of the plot, by default 1.0
    xunits: str, optional
        The text for the x axis units
    stimbar: dict{sound_params: None, waveform: None}
        A dictionary for placing a plot of the stimulus under the plot.
        Keys: 'sound_params', 'waveform'
        if 'sound_params' is not None, then it should be a sound_params dataclass
        from make_sound_waveforms, and will be used to generate the waveform.
        if 'waveform' is not None, then it is assumed that the passed
        argument is a list [time, waveform].

    Returns
    -------
    Nothing
    """
    num_trials = run_info.nReps
    spike_times_flat_isi = VAFlatten.flatten_spike_array(spike_times, 
                            time_window=(zero_time, max_time), isi_flag=True)

    bins = np.arange(0.0, max_time - zero_time, bin_width)
    if bin_fill:
        face_color = edge_color
    else:
        face_color = None  # "None"
    if xunits == "radians":
        xu = radians
        bins = bins * radians
    else:
        xu = None
    if (
        (not np.isnan(np.sum(spike_times_flat_isi)))
        and (len(spike_times_flat_isi) > 0)
        and (ax is not None)
    ):
        ax.hist(
            x=spike_times_flat_isi,
            bins=bins,
            density=False,
            weights=(1./(bin_width*num_trials)) * np.ones_like(spike_times_flat_isi),
            histtype="stepfilled",
            facecolor=face_color,
            edgecolor=edge_color,
            linewidth=0,
            alpha=alpha,
        )
    else:
        if ax is not None:
            ax.text(
                0.5,
                0.5,
                "No Spikes",
                fontsize=14,
                color="r",
                transform=ax.transAxes,
                horizontalalignment="center",
            )

    return

@TRC()
def print_AN_rates(
    spike_times: Union[list, np.ndarray],
    run_info: object,
):
    SPKANA.ANfspike(
        spike_times,
        stime = run_info.pip_start,
        nReps=run_info.nReps,
        stimdur = run_info.pip_duration
    )


@TRC()
def plot_fsl_ssl(
    spike_times: Union[list, np.ndarray],
    run_info: object,
    max_time: float,
    min_time: float = 0.0,
    fsl_win: Union[None, tuple] = None,
    bin_width: float = 1e-3,
    ax: Union[object, None] = None,
    zero_time: float = 0.0,
    offset: float = 0.0,
    cellID: Union[int, None, str] = None,
    show_values: bool = True,
) -> Tuple[object, object]:
    """Plot the first and second spike latencies as histogram into the selected axis
    All times should be specified in *seconds*

    Parameters
    ----------
    spike_times : Union[list, np.ndarray]
        A list of the spike times (in seconds)
    run_info : object
        The model_params run_info dataclass associated with this data
    max_time : float
        Maximum time for the FSL/SSL plot
    min_time : float, optional
        Minimum time for the fsl/ssl plot, by default 0.0
    fsl_win : Union[None, tuple], optional
        Time window for identifying the FSL, by default None
    bin_width : float, optional
        Width of bins for FSL/SSL plot, by default 1e-3
    ax : Union[object, None], optional
        matplotlib axis instance for this plot, by default None
    zero_time : float, optional
        The onset time for the stimulus for reference, by default 0.0
    offset : float, optional
        _description_, by default 0.0
    cellID : Union[int, None, str], optional
        The ID for this cell (Not used!), by default None
    show_values : bool, optional
        Flag to place the values in text fields on the plot, by default True

    Returns
    -------
    Tuple[fsl object, ssl object]
        First and second spike latency arrays
    """

    num_trials = run_info.nReps  # self._ntrials(spike_times)

    fsl, ssl = SPKANA.CNfspike(
        spike_times,
        run_info.pip_start - zero_time,
        nReps=run_info.nReps,
        fsl_win=fsl_win,
    )

    fsl *= 1e3  # better to show FSL/SSL in milliseconds
    ssl *= 1e3
    # if cellID is None:
    #     index_row = self.parent.selected_index_rows[0]
    #     selected = self.parent.table_manager.get_table_data(index_row)
    #     print(
    #         f"Cell: '{int(self.parent.cellID):d}','{selected.synapseExperiment:s}',",
    #         end="",
    #     )
    # else:
    #     print(f"Cell: '{int(cellID):d}','{selected.synapseExperiment:s}',", end="")

    if len(fsl) > 0:
        print(f"{np.nanmean(fsl):.3f},{np.nanstd(fsl):.3f},", end="")
    else:
        print("No FSL spike data. ", end="")
    if len(ssl) > 0:
        print(f"{np.nanmean(ssl):.3f},{np.nanstd(ssl):.3f}")
    else:
        print("No SSL spike data")
    if ax is not None:
        latency_hbins = np.arange(min_time, max_time, bin_width)  # use msec here
        if (len(fsl)) > 0.0:
            ax.hist(
                fsl,
                bins=latency_hbins,
                facecolor="b",
                edgecolor="b",
                alpha=0.6,
                bottom=offset,
            )
        if (len(ssl)) > 0.0:
            ax.hist(
                ssl,
                bins=latency_hbins,
                facecolor="r",
                edgecolor="r",
                alpha=0.6,
                bottom=offset,
            )
        ax.set_xlim(min_time, max_time)
        if show_values:
            if len(fsl) > 0:
                fsl_text = f"FSL: {np.nanmean(fsl):.3f} (SD {np.nanstd(fsl):.3f})"
            else:
                fsl_text = "No spikes for FSL measure. "
            if len(ssl) > 0:
                fsl_text += f"\nSSL: {np.nanmean(ssl):.3f} (SD {np.nanstd(ssl):.3f})"
            else:
                fsl_text = "No spikes for SSL measure. "
            ax.text(
                0.30,
                0.95,
                fsl_text,
                fontsize=7,
                color="k",
                # fontfamily="monospace",
                transform=ax.transAxes,
                horizontalalignment="left",
                verticalalignment="top",
            )

    return (fsl, ssl)


def plot_spike_raster(
    P: object,
    mode: str='postsynaptic',
    n_inputs: int = 0,
    i_trial: int = 0,
    n_trials: int = 10, 
    spike_times: List = [],
    data_window: List = [],
    panel: Union[str, None] = None,
):
    """Plot spike rasters

    Parameters
    ----------
    P : object
        pylibrary plothelpers Plot object
    mode : str
        either "presynaptic" or "postsynaptic"
        changes the way the spike raster is plotted.
        for postysynaptic, just the spikes in each trial
        for presynaptic, each input is plotted and trials are grouped.
    n_inputs : int, optional
        number of inputs for ell, by default 0
    i_trial : int, optional
        trial number, by default 0. Only used if mode is "presynaptic"
    n_trials : int, optional (default 10)
        number of trials, by default 10. 
    spike_times : List
        A list of spike times or a list of lists of spike times, by default []
    data_window : List, optional
        _description_, by default []
    panel : str, optional
        Panel label where the data will be plotted, by default ""
    """
    assert mode in ["presynaptic", "postsynaptic"]
    assert panel is not None
    # Raster of spikes for panel
    plot_dur = np.fabs(np.diff(data_window))

    if mode == "postsynaptic":
        for j in range(n_trials):
            ispt = [
                i
                for i in range(len(spike_times[j]))
                if spike_times[j][i] >= data_window[0] and spike_times[j][i] < data_window[1]
            ]
            P.axdict[panel].plot(
                np.array(spike_times[j][ispt]) - data_window[0],
                j * np.ones(len(ispt)),
                "|",
                markersize=1.5,
                color="b",
            )
        P.axdict[panel].set_xlim(0, plot_dur)
    elif mode == "presynaptic":
        for k in range(n_inputs):  # raster of input spikes
            if np.max(spike_times[k]) > 2.0:  # probably in miec...
                tk = np.array(spike_times[k]) * 1e-3
            else:
                tk = np.array(spike_times[k])

            y = (i_trial + 0.1 + k * 0.05) * np.ones(len(tk))
            in_spt = [
                i
                for i in range(tk.shape[0])
                if (tk[i] >= data_window[0]) and (tk[i] < data_window[1])
            ]
            y = (i_trial + 0.1 + k * 0.05) * np.ones(len(in_spt))
            if i_trial % 10 == 0:
                P.axdict[panel].plot(
                    tk[in_spt] - data_window[0],
                    y,
                    "|",
                    markersize=2.5,
                    color="k",
                    linewidth=0.5,
                )
        P.axdict[panel].set_xlim(0, plot_dur)

def plot_spiketrain_raster(
        spike_times, max_raster: int = 20, ax=None, plot_win: tuple = (0, 0.25)
    ):
    # print(len(spike_times))
    # print(plot_win)
    n_trials = len(spike_times)
    if n_trials > max_raster:
        n_trials = max_raster
    for i in range(n_trials):
        ispt = [
            j
            for j in range(len(spike_times[i]))
            if spike_times[i][j] >= plot_win[0] and spike_times[i][j] < plot_win[1]
        ]
        ax.plot(
            np.array(spike_times[i])[ispt] - plot_win[0],
            i * np.ones(len(ispt)),
            "|",
            markersize=1.5,
            color="b",
        )

def plot_stacked_spiketrain_rasters(
    spike_times_by_input,
    ax,
    si,
    syninfo,
    plot_win: tuple = (0, 0.25),
    max_trials=2,
    use_colors: bool = True,
    colormap="Set3",
    cbar_vmax: float = 300.0,
    linewidth=1.0,
):
    """
    Spike trains are plotted as a raster for all inputs in the AN data

    Parameters
    ----------
    spike_times_by_input : list
        list by input of spike times by trial

    ax : matplotlib axis to place the plat

    si : Params

    syninfo : synaptic info from SC, syninfo = plot_sims.get_synaptic_info(cell_n)
    plot_win : tuple
        time window of data to show

    max_trials : int
        number of trials to show (the first max_trials are plotted)

    use_colors : bool (default True)
        Plot the raster ticks with colors scaled by input surface area

    colormap : str
        color map to use for plotting when use_colors is True

    cbar_vmax : float
        maximum for the colorbar scale

    """
    n_inputs = len(spike_times_by_input)
    n_trials = len(spike_times_by_input[0])
    trial_spc = 1.0
    input_spc = 0.6 / n_inputs

    if use_colors:
        cmx = sns.color_palette(colormap, as_cmap=True)
        # cell_n = int(si.cellID[-2:])
        
        syn_ASA = np.array([syninfo[1][isite][0] for isite in range(n_inputs)])
        # max_ASA = np.max(syn_ASA)
    print("syn asa: ", syn_ASA)
    for i in range(n_trials):  # and by trial
        if i > max_trials - 1:
            continue
        for k in range(n_inputs):  # raster of input spikes by input
            if np.max(spike_times_by_input[k][i]) > 2.0:  # probably in msec...
                tk = np.array(spike_times_by_input[k][i]) * 1e-3
            else:
                tk = np.array(spike_times_by_input[k][i])

            in_spt = [
                i
                for i in range(tk.shape[0])
                if (tk[i] >= plot_win[0]) and (tk[i] < plot_win[1])
            ]
            y = (
                ((i + 1) * trial_spc)
                - ((n_inputs / 2.0) * input_spc)
                + (k * input_spc)
            ) * np.ones(len(in_spt))
            if use_colors:
                color = cmx.colors[int(cmx.N * syn_ASA[k] / cbar_vmax) - 1]
            else:
                color = "k"
            ax.plot(
                tk[in_spt] - plot_win[0],
                y,
                "|",
                markersize=2.5,
                color=color,
                linewidth=linewidth,
            )
            ax.yaxis.set_ticks([1, 2], ["1", "2"])
    ax.set_xlim(plot_win)