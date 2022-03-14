"""
"""
from typing import Union, Tuple

import matplotlib.pyplot as mpl
import numpy as np
from vcnmodel.util import trace_calls
from vcnmodel.analyzers import analysis as SPKANA

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
        Weights for the bins, optional and not used, by default 1.0
    bin_fill : bool, optional
        Fill the bins, or leave open (transparent), by default True
    edge_color : str, optional
        color of the edges of the plot bars, by default "k"
    alpha : float, optional
        Transparency of the plot, by default 1.0

    Returns
    -------
    Nothing
    """
    num_trials = run_info.nReps
    bins = np.arange(0.0, max_time - zero_time, bin_width)
    spf = []
    for x in spike_times:
        if isinstance(x, list):
            spf.extend(x)
        else:
            spf.extend([x])
    spike_times_flat = np.array(spf, dtype=object).ravel() - zero_time
    # h, b = np.histogram(spike_times_flat, bins=bins)
    bins = np.arange(0.0, max_time - zero_time, bin_width)
    if bin_fill:
        face_color = "k"
    else:
        face_color = "None"
    if (
        (not np.isnan(np.sum(spike_times_flat)))
        and (len(spike_times_flat) > 0)
        and (ax is not None)
    ):
        ax.hist(
            x=spike_times_flat,
            bins=bins,
            density=False,
            # weights=scale * np.ones_like(spike_times_flat),
            histtype="stepfilled",
            facecolor=face_color,
            edgecolor=edge_color,
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
    return  # h, b


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
    All times should be specified in seconds

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
