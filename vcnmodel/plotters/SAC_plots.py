
from pathlib import Path
from typing import List, Union

import matplotlib.pyplot as mpl
import numpy as np
import pandas as pd
import seaborn as sns
from lmfit.models import GaussianModel
from matplotlib import collections as mc
from pylibrary.plotting import plothelpers as PH
from pylibrary.tools import cprint as CP

import vcnmodel.group_defs as GRPDEF
import vcnmodel.plotters.figure_data as FD
import vcnmodel.plotters.plot_sims as plot_sims
import vcnmodel.util.get_data_paths as get_data_paths
import vcnmodel.util.readmodel as readmodel
from vcnmodel.analyzers import sac as SAC
from vcnmodel.util.get_data_paths import get_data_paths
from vcnmodel.util.set_figure_path import set_figure_path

ReadModel = readmodel.ReadModel()
cprint = CP.cprint

config = get_data_paths()
# print(config)
# os.chdir(config['codeDirectory'])
PS = plot_sims.PlotSims(parent=None)
sacs = FD.figure_SAC
print(sacs)

# seaborn default palette, first 3 colors
colors = [
    (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
    (1.0, 0.4980392156862745, 0.054901960784313725),
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
]

expts = ["all", "largestonly", "removelargest"]


def one_sac(cell_id, protocol, pname):
    """Compute One SAC for the specified cell and the experiment (pname)

    Parameters
    ----------
    cell_id : int
        cell identification by number
    protocol : str
        protocol name
    pname : str
        experiment name (one of expts above)

    Returns
    -------
    dict or None
        A dictionary of the results, or None if there was no data.
    """

    plot_win = [0.1, 1.0]
    plot_dur = np.fabs(np.diff(plot_win))
    time_scale = 1.0
    gbc = f"VCN_c{int(cell_id):02d}"
    simfile = None
    if pname in list(sacs[cell_id].keys()):
        simfile = sacs[cell_id][pname]
    else:
        return None

    basefn = f"{config['cellDataDirectory']:s}/{gbc:s}/Simulations/AN/"
    fndir = Path(basefn, simfile)
    files = list(fndir.glob("*.p"))
    filename = Path(fndir, files[0])
    PD = plot_sims.PData(GRPDEF.gradeACells)
    print(f"Getting data for gbc: {gbc:s}")
    SC, syninfo = PS.get_synaptic_info(gbc)
    # mtime = Path(fn).stat().st_mtime
    # timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime(
    #     "%Y-%m-%d-%H:%M"
    # )

    # changetimestamp = plot_sims.get_changetimestamp()

    # res = PS.get_data(fn, PD, changetimestamp, protocol)
    model_data = ReadModel.get_data(fn=filename, PD=PD, protocol=protocol)
    if model_data is None:
        print("no file?")
        return
    data = model_data.data
    si = model_data.SI
    ri = model_data.RI
    AR = model_data.AR
    SP = model_data.SP
    RM = model_data.RM
    RCP = model_data.RCP
    RCD = model_data.RCD
    PD.thiscell = gbc

    RCP.si = si
    RCP.ri = ri
    (
        totaldur,
        soundtype,
        pip_start,
        pip_duration,
        F0,
        dB,
        fmod,
        dmod,
    ) = PS.get_stim_info(si, ri)
    all_bu_st = []
    all_bu_st_trials = []
    ntr = len(AR.MC.traces)  # number of trials
    for i in range(ntr):  # for alls trials in the measure.
        time_base = AR.MC.time_base / 1000.0  # convert to seconds
        trd = data["Results"][i]
        dt = si.dtIC / 1000.0  # convert from msec to seconds
        idx = (int(plot_win[0] / dt), int(plot_win[1] / dt))
        if i == 0:
            waveform = trd["stimWaveform"]
        stb = trd["stimTimebase"]  # convert to seconds
        stim_dt = stb[1] - stb[0]
        if i == 0:
            n_inputs = len(trd["inputSpikeTimes"])
        if (
            len(trd["spikeTimes"]) > 0 and np.nanmax(trd["spikeTimes"]) > 2.0
        ):  # probably in msec
            time_scale = 1e-3  # so scale to seconds
        sptimes = np.array(trd["spikeTimes"]) * time_scale  # convert to seconds
        if not isinstance(trd["spikeTimes"], list) and not isinstance(
            trd["spikeTimes"], np.ndarray
        ):
            cprint("r", "spiketimes is not list")
            cprint("r", f"    {type(trd['spikeTimes'])=}")
            return
        all_bu_st.extend(sptimes)
        all_bu_st_trials.append(sptimes)
        ispt = [  # plot spike times in the SAC analysis window
            i
            for i in range(len(sptimes))
            if sptimes[i] >= pip_start and sptimes[i] < pip_duration - pip_start
        ]
        # P.axdict["B"].plot(
        #     np.array(sptimes[ispt]),
        #     i * np.ones(len(ispt)),
        #     "|",
        #     markersize=1.5,
        #     color="b",
        # )
        w_tb = np.linspace(0.0, stim_dt * len(waveform), num=len(waveform))

        i_wpt = np.where((w_tb > pip_start) & (w_tb <= pip_duration))[0]
        # P.axdict["C"].plot(w_tb[i_wpt], waveform[i_wpt], linewidth=0.33)
    if ri.soundtype.endswith("Clicks"):
        print("Clickpars")
        pars = {
            "twin": 0.002,
            "binw": 3 * dt,
            "delay": ri.clickStart + 0.2 * ri.clickTrainDuration,
            "dur": 0.8 * ri.clickTrainDuration,
            "displayDuration": 0.002,
            "nrep": len(all_bu_st_trials),
            "baseper": 1e-3 * 1.3333333333333333,
            "stimdur": pip_duration * 0.8,
            "maxn": 100000000,
        }

    else:
        print("Sam Pars")
        pars = {
            "twin": 0.020,
            "binw": 3 * dt,
            "delay": pip_start + 0.2 * pip_duration,
            "dur": 0.8 * pip_duration,
            "displayDuration": 0.050,
            "nrep": len(all_bu_st_trials),
            "baseper": 1e-3 * 1.3333333333333333,
            "stimdur": pip_duration * 0.8,
            "maxn": 100000000,
        }
    print(ri.soundtype)
    sac = SAC.SAC()
    yh, bins = sac.SAC_with_histo(
        all_bu_st_trials, pars=pars, engine="cython", dither=dt / 2.0
    )
    sac_bu_CI, peaks, HW, FW = sac.SAC_measures(yh, bins)
    if not np.isnan(sac_bu_CI):
        sac_bu_hw = HW[0][0] * pars["binw"]
        print("BU SAC Report: \n  ")
        print(
            f"    HW:    {sac_bu_hw:.6f}  CI: {sac_bu_CI:.2f}  Left IPS: {HW[2][0]:.2f}  Right IPS: {HW[3][0]:.2f}"
        )
        if ri.soundtype.endswith("Clicks"):
            sac_label = f"Expt: {ri.Spirou:14s} {ri.dB:3.0f} dBSPL  HW={1e3*sac_bu_hw:.3f} ms  CI={sac_bu_CI:6.2f}"
        else:
            sac_label = f"Expt: {ri.Spirou:14s} {ri.dB:3.0f} dBSPL Fmod={ri.fmod:5.1}fHz Dmod={ri.dmod:5.1f}%"
    else:
        print("No spikes ")
        sac_bu_hw = np.nan
    res = {
        "CI": sac_bu_CI,
        "HW": sac_bu_hw,
        "yh": yh,
        "bins": bins,
        "simfile": simfile,
        "gbc": gbc,
    }
    return res


def run_sac_analysis():
    protocol = "runANPSTH"
    rowcount = 0
    max_CI = 0.0
    df = pd.DataFrame(
        data=None, columns=["Cell", "Expt", "CI", "Halfwidth", "SAC", "bins"]
    )

    for cell_number in list(sacs.keys()):
        for ik, expt in enumerate(expts):
            res = one_sac(cell_number, protocol, expt)
            igbc = int(res["gbc"][-2:])
            if res["CI"] == 0:
                res["CI"] = np.nan
            if res is not None:
                if res["CI"] == 0:
                    res["CI"] = np.nan

                print("gbc: ", res["gbc"], " simfile: ", res["simfile"])
                # P.axarr[rowassign[igbc][0], rowassign[igbc][1]].plot(res["bins"], res["yh"], color=colors[ik])
                if res["CI"] > max_CI:
                    max_CI = res["CI"]
                df = df.append(
                    {
                        "Cell": cell_number,
                        "Expt": expt,
                        "CI": res["CI"],
                        "Halfwidth": res["HW"],
                        "SAC": res["yh"],
                        "bins": res["bins"],
                    },
                    ignore_index=True,
                )
            else:
                df = df.append(
                    {
                        "Cell": cell_number,
                        "Expt": expt,
                        "CI": np.nan,
                        "Halfwidth": np.nan,
                        "SAC": res["yh"],
                        "bins": res["bins"],
                    },
                    ignore_index=True,
                )

    return df


def do_plot(df, sacs, max_CI: float = 60.0, fits:bool=False)->object:
    """Make a supplemental plot of the SAC results for all cells

    Parameters
    ----------
    df : Pandas dataframe

    sacs : SAC data
        _description_
    max_CI : float, optional
        largest scale factor, by default 5.0
    """
    import string

    from matplotlib.gridspec import GridSpec

    labels = [s for s in string.ascii_uppercase]
    rowassign = {
        2: [0, 0],
        5: [1, 0],
        6: [2, 0],
        9: [3, 0],
        10: [4, 0],
        11: [0, 1],
        13: [1, 1],
        17: [2, 1],
        18: [3, 1],
        30: [4, 1],
    }
    # fig = mpl.figure(constrained_layout = True, figsize=(10, 8))
    # fig.set_tight_layout({"pad": 0.15})
    # subfigs = fig.subfigures(1, 2, wspace=0.15, hspace=0.08)
    # subfigs[1].figure.set_tight_layout({"pad":0.08})
    # axes_left = subfigs[0].subplots(5, 2, sharey=True)
    # axes_right = subfigs[1].subplots(2, 1)
    # print(dir(subfigs[1]))

    P = PH.regular_grid(
        5,
        2,
        order="columnsfirst",
        figsize=(10.0, 8.0),
        showgrid=False,
        verticalspacing=0.05,
        horizontalspacing=0.1,
        margins={
            "bottommargin": 0.1,
            "leftmargin": 0.07,
            "rightmargin": 0.35,
            "topmargin": 0.15,
        },
        labelposition=(-0.15, 1.02),
        parent_figure=None,
        panel_labels=labels[:10],
    )
    ax = P.axarr.ravel()
    P2 = PH.regular_grid(
        2,
        1,
        showgrid=False,
        verticalspacing=0.1,
        horizontalspacing=0.1,
        labelposition=(-0.15, 1.0),
        margins={
            "bottommargin": 0.1,
            "leftmargin": 0.72,
            "rightmargin": 0.05,
            "topmargin": 0.15,
        },
        parent_figure=P,
        panel_labels=labels[10:12],
    )
    ax2 = P2.axarr.ravel()
    for ij, j in enumerate(sorted(list(rowassign.keys()))):
        igbc = int(j)
        print("gbc: ", igbc, ij)
        this_ax = ax[ij]  # P.axarr[rowassign[igbc][0], rowassign[igbc][1]]
        for ik, k in enumerate(expts):
            if igbc == 30:
                label = k
            else:
                label = None
            bins = df[(df["Cell"] == j) & (df["Expt"] == k)]["bins"]
            yh = df[(df["Cell"] == j) & (df["Expt"] == k)]["SAC"]
            this_ax.plot(
                bins.values[0][:-1] * 1e3, yh.values[0], color=colors[ik], label=label
            )
            if fits:
                if len(df[(df["Cell"] == j) & (df["Expt"] == k)]["gfitx"]) > 0:
                    fitbins = df[(df["Cell"] == j) & (df["Expt"] == k)]["gfitx"].values[0]
                    fitgauss = df[(df["Cell"] == j) & (df["Expt"] == k)]["gfity"].values[0]
                    for rep in ["[", "]", "\n"]:
                        fitbins = fitbins.replace(rep, "")
                        fitgauss = fitgauss.replace(rep, "")
                    fitbins = np.fromstring(fitbins, sep=" ")
                    fitgauss = np.fromstring(fitgauss, sep=" ")
                    this_ax.plot(
                        fitbins * 1e3,
                        fitgauss,
                        color="r",  # colors[ik],
                        linestyle="--",
                        linewidth=0.5,
                    )
            this_ax.set_ylim(0, max_CI)
            this_ax.set_ylabel("CI")

            # this_ax.text(-0.05, 1.05, s=labels[ij], transform=this_ax.transAxes,
            #     va="bottom", ha="right", fontsize=12)
            this_ax.text(
                0.5,
                1.0,
                s=f"BC{igbc:02d}",
                transform=this_ax.transAxes,
                va="bottom",
                ha="center",
                fontsize=10,
            )
        if ij in [4, 9]:
            this_ax.set_xlabel("T (ms)")
        if ij < 5:
            this_ax.set_ylabel("CI")
        PH.nice_plot(this_ax, direction="outward", ticklength=4)
        if j == 0:
            this_ax.legend(fontsize=7)

    df["HW (ms)"] = df["Halfwidth"].values * 1000.0
    df["HW (ms)"] = df["HW (ms)"].replace(0, np.nan)
    hue_order = expts
    spalette = {}
    for ic, h in enumerate(hue_order):
        spalette[h] = colors[ic]

    sns.stripplot(
        data=df,
        x="Cell",
        y="HW (ms)",
        hue="Expt",
        ax=ax2[0],
        jitter=0.5,
        dodge=True,
        size=5,
        palette=spalette,
    )
    ax2[0].set_ylim(0, 1.0)
    # legend = ax2[0].legend(
    #         bbox_to_anchor=(0.45, 0.25),
    #         loc="upper left",
    #         borderaxespad=0,
    #         fontsize=9,
    #         markerscale=0.5,
    #         frameon=True,
    #         edgecolor="black",
    #     )

    sns.stripplot(
        data=df,
        x="Cell",
        y="CI(g)",
        hue="Expt",
        ax=ax2[1],
        jitter=0.5,
        dodge=True,
        size=5,
        palette=spalette,
    )
    ax2[1].set_ylim(0, 75)
    ax2[1].set_ylabel("CI")
    ax2[1].get_legend().remove()
    # reorder the legends
    handles, labels = ax2[0].get_legend_handles_labels()
    legdict = {"all": "All Inputs", "largestonly": "Largest input only", 
    "removelargest": "Largest input removed", "removetwolargest": "Two largest inputs removed"}
    labels = [legdict[label] for label in labels]
    order = [0, 2, 1]
    ax2[0].legend([handles[idx] for idx in order],[labels[idx] for idx in order],
            bbox_to_anchor=(0.5, 0.01),
            loc="lower center",
            borderaxespad=0,
            fontsize=9,
            markerscale=0.5,
            frameon=True,
            edgecolor="black",
            )
    # legdict = {"all": "All Inputs", "largestonly": "Largest input only", 
    # "removelargest": "Largest input removed", "removetwolargest": "Two largest inputs removed"}
    # ltexts = ax2[0].legend().get_texts()
    # for ik, this_l in enumerate(ltexts):
    #     labeltext = this_l.get_text()
    #     if labeltext in list(legdict.keys()):
    #         this_l.set_text(legdict[labeltext])

    for ax in ax2:
        PH.nice_plot(ax, direction="outward", ticklength=4)
    return P
    mpl.savefig(
        Path(
            config["baseDataDirectory"],
            "Figures",
            "Figure6",
            "Figure6_supp",
            "Figure6_Supplemental3_SAC_Clicks_SynapseConfigs_V2.pdf",
        ),
        metadata={
            "Creator": "Paul Manis",
            "Author": "Paul Manis",
            "Title": "Figure 6 Supplemental Figure 3 V2",
        },
    )
    mpl.show()
    return P


# sns.barplot(data=df, x="Expt", y="CI", hue="Cell")
def sac2plot(df):
    df["HWms"] = df["Halfwidth"].values * 1000.0
    df["HWms"] = df["HWms"].replace(0, np.nan)
    PB = PH.regular_grid(
        1,
        2,
        order="columnsfirst",
        figsize=(6.0, 4.0),
        showgrid=False,
        verticalspacing=0.05,
        horizontalspacing=0.1,
        margins={
            "bottommargin": 0.1,
            "leftmargin": 0.07,
            "rightmargin": 0.05,
            "topmargin": 0.15,
        },
        labelposition=(0.0, 0.0),
        parent_figure=None,
        panel_labels=None,
    )

    hue_order = expts

    sns.stripplot(
        data=df,
        x="Cell",
        y="HWms",
        hue="Expt",
        ax=PB.axdict["A"],
        jitter=0.5,
        dodge=True,
        size=5,
    )
    # sns.lineplot(data=df, x="Cell", y="HWms", hue="Expt", hue_order=hue_order, ax=PB.axdict["A"])

    # lines = ([[x, n] for n in group] for x, (_, group) in enumerate(df.groupby(['Cell'], sort = False)['HWms']))
    # lc = mc.LineCollection(lines, colors='red', linewidths=1)
    # PB.axdict["A"].add_collection(lc)

    sns.stripplot(
        data=df,
        x="Cell",
        y="CI(g)",
        hue="Expt",
        ax=PB.axdict["B"],
        jitter=0.5,
        dodge=True,
        size=5,
    )
    mpl.show()


def remeasure(df, twin: float = 2):
    rowassign = {
        2: [0, 0],
        5: [1, 0],
        6: [2, 0],
        9: [3, 0],
        10: [4, 0],
        11: [0, 1],
        13: [1, 1],
        17: [2, 1],
        18: [3, 1],
        30: [4, 1],
    }
    sac = SAC.SAC()

    for ij, j in enumerate(sorted(list(rowassign.keys()))):
        igbc = int(j)
        for ik, k in enumerate(["all", "removelargest", "largestonly"]):
            data_pointer = df[(df["Cell"] == j) & (df["Expt"] == k)]
            bins = data_pointer["bins"]
            bins = np.array(bins.values[0])
            yh = data_pointer["SAC"]
            yh = np.array(yh.values[0])
            if np.max(yh) > 0.0:  # must have a peak to work with
                win_indx = np.argwhere((-twin <= bins) & (bins < twin))[:-1].T
                win_indx = [s for s in win_indx][0]
                win_data = yh[win_indx]
                win_bins = bins[win_indx]
                model = GaussianModel()
                pars = model.guess(win_data, x=win_bins)
                out = model.fit(win_data, pars, x=win_bins)

                height = (
                    0.3989423
                    * out.best_values["amplitude"]
                    / max(1e-15, out.best_values["sigma"])
                )
                fwhm = 2.3548200 * out.best_values["sigma"]
                df.loc[(df["Cell"] == j) & (df["Expt"] == k), "Halfwidth"] = fwhm
                df.loc[(df["Cell"] == j) & (df["Expt"] == k), "CI(g)"] = height
                df.loc[(df["Cell"] == j) & (df["Expt"] == k), "gfity"] = str(
                    out.best_fit
                )
                df.loc[(df["Cell"] == j) & (df["Expt"] == k), "gfitx"] = str(win_bins)

            else:
                df.loc[(df["Cell"] == j) & (df["Expt"] == k), "Halfwidth"] = np.nan
                df.loc[(df["Cell"] == j) & (df["Expt"] == k), "CI(g)"] = np.nan
                df.loc[(df["Cell"] == j) & (df["Expt"] == k), "gfity"] = ""
                df.loc[(df["Cell"] == j) & (df["Expt"] == k), "gfitx"] = ""

    return df

def plot_sacs(save_fig:bool=True, figinfo: Union[object, None] = None, fits:bool=False):
    pd_file = Path(
        config["baseDataDirectory"],
        # config["figureDirectory"],
        config["figureIntermediateDirectory"],
        "SAC_Results_Clicks.pkl",
    )
    # df = run_sac_analysis()
    # df.to_pickle(pd_file)
    df = pd.read_pickle(pd_file)
    df["gfitx"] = "nan"  # df['gfitx'].astype(object)
    df["gfity"] = "nan"  # df['gfity'].astype(object)
    df2 = remeasure(df)

    P = do_plot(df2, sacs, fits=fits)
    # print(sacs)
    # sac2plot(df)

    if save_fig:
        figinfo.P = P
        figinfo.show_name = False
        figinfo.filename = set_figure_path(
            fignum=6, filedescriptor="SAC_Clicks_SynapseConfigs", suppnum=3
        )
        # "Figure7_Supplemental3_SAC_Clicks_SynapseConfigs.pdf",
        figinfo.title[
            "title"
        ] = "SBEM Project Figure 6 Modeling: Supplemental 3: SAC_Clicks_SynapseConfigs"
        title2 = {"title": f"", "x": 0.99, "y": 0.01}
        return figinfo
    else:
        mpl.show()
        return None

def main():
    plot_sacs(save_fig=False)
    pass  # keep sphinx from running this

if __name__ == "__main__":
    main()
