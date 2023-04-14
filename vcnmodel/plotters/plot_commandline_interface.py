#######################################################
# plotting from the command line - limited... 
#######################################################
import argparse
import argparse

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Tuple, Union
from vcnmodel.util.get_data_paths import get_data_paths
import numpy as np
import matplotlib.pyplot as mpl
from plot_sims import PlotSims
import vcnmodel.group_defs as GRPDEF
#  from plot_sims import grAList
from plot_sims import get_changetimestamp
from pylibrary.plotting import plothelpers as PH

@dataclass
class PData:
    """
    data class for some parameters that control what we read
    """

    gradeA: list = field(default_factory=GRPDEF.grAList())
    default_modelName: str = "XM13_nacncoop"
    soma_inflate: bool = True
    dend_inflate: bool = True
    basepath: str = ""  # config["baseDataDirectory"]
    renderpath: str = ""  # " str(Path(self.config["codeDirectory"], "Renderings"))
    thiscell: str = ""


modeltypes = ["mGBC", "XM13", "RM03", "XM13_nacncoop", "XM13A_nacncoop"]
runtypes = ["AN", "an", "IO", "IV", "iv", "gifnoise"]
experimenttypes = [
    None,
    "delays",
    "largestonly",
    "removelargest",
    "mean",
    "allmean",
    "twolargest",
]
modetypes = ["find", "singles", "IO", "multi"]
analysistypes = ["traces", "PSTH", "revcorr", "SAC", "tuning", "singles"]
dendriteChoices = [
    "normal",
    "passive",
    "active",
]

def build_parser():
    """
    create the command line parser and process the line
    """

    parser = argparse.ArgumentParser(description="Plot GBC results")
    parser.add_argument(
        dest="cell",
        action="store",
        nargs="+",
        type=str,
        default=None,
        help="Select the cell(s) or 'A' for all(no default)",
    )
    parser.add_argument(
        "-p",
        "--protocol",
        dest="protocol",
        action="store",
        default="IV",
        choices=runtypes,
        help=("Select the protocol (default: IV) from: %s" % runtypes),
    )
    parser.add_argument(
        "-a",
        "--analysis",
        dest="analysis",
        action="store",
        default=analysistypes[0],
        choices=analysistypes,
        help=(
            f"Select the analysis type (default: {analysistypes[0]:s}) from: {str(analysistypes):s}"
        ),
    )
    parser.add_argument(
        "-M",
        "--modeltype",
        dest="modeltype",
        action="store",
        default="XM13_nacncoop",
        help=("Select the model type (default XM13_nacncoop) from: %s " % modeltypes),
    )

    parser.add_argument(
        "-s",
        "--scaled",
        dest="scaled",
        action="store_true",
        help=("use scaled data or not"),
    )

    parser.add_argument(
        "-e",
        "--experiment",
        dest="experiment",
        action="store",
        default="delays",
        choices=experimenttypes,
        help=("Select the experiment type from: %s " % experimenttypes),
    )

    parser.add_argument(
        "--dendritemode",
        dest="dendriteMode",
        default="normal",
        choices=dendriteChoices,
        help="Choose dendrite table (normal, active, passive)",
    )

    parser.add_argument(
        "-d",
        "--dB",
        dest="dbspl",
        type=float,
        action="store",
        default=None,
        help=("Select the models at specific intensity"),
    )

    parser.add_argument(
        "-r",
        "--nreps",
        dest="nreps",
        type=int,
        action="store",
        default=None,
        help=("Select the models with # reps"),
    )

    parser.add_argument(
        "-c",
        "--check",
        dest="check",
        action="store_true",
        help=("Just check selection criteria and return"),
    )

    return parser


def cmdline_display(args, PD):
    """
    Display analysis and traces
    when called from the commandline
    """

    PS = PlotSims(parent=None)
    config = get_data_paths()
    args.protocol = args.protocol.upper()
    changetimestamp = get_changetimestamp()
    # PD.gradeA = [cn for cn in args.cell]
    print("args.cell: ", args.cell)
    if args.cell[0] in ["A", "a"]:
        pass  # (all grade a cells defined in pd dataclass)
    else:
        PD.gradeA = [int(c) for c in args.cell]
    rows, cols = PH.getLayoutDimensions(len(PD.gradeA))
    plabels = [f"VCN_c{g:02d}" for g in PD.gradeA]
    for i in range(rows * cols - len(PD.gradeA)):
        plabels.append(f"empty{i:d}")

    sizex = cols * 3
    sizey = rows * 2.5
    if rows * cols < 4:
        sizex *= 2
        sizey *= 2

    chan = "_"  # for certainity in selection, include trailing underscore in this designation

    modelName = args.modeltype
    print("Scaling: ", args.scaled)
    stitle = "unknown scaling"
    if args.scaled:
        PD.soma_inflate = True
        PD.dend_inflate = True
        stitle = "scaled"
    else:
        PD.soma_inflate = False
        PD.dend_inflate = False
        stitle = "notScaled"

    if args.analysis == "tuning":

        fng = []
        for ig, gbc in enumerate(PD.gradeA):
            basefn = f"{config['cellDataDirectory']:s}/VCN_c{gbc:02d}/Simulations/{args.protocol:s}/"
            pgbc = f"VCN_c{gbc:02d}"
            allf = list(Path(basefn).glob("*"))
            allfiles = sorted([f for f in allf if f.is_dir()])
            # print(allfiles)
            fnlast3 = allfiles[-3:]
        for f in fnlast3:
            fn = list(f.glob("*"))
            fng.append(fn[0])
        P = PS.plot_tuning(args, filename=None, filenames=fng)

    else:
        for ig, gbc in enumerate(PD.gradeA):
            basefn = f"{config['cellDataDirectory']:s}/VCN_c{gbc:02d}/Simulations/{args.protocol:s}/"
            pgbc = f"VCN_c{gbc:02d}"
            if args.protocol == "IV":
                name_start = f"IV_Result_VCN_c{gbc:02d}_inp=self_{modelName:s}*.p"
                args.experiment = None
            elif args.protocol == "AN":
                name_start = f"AN_Result_VCN_c{gbc:02d}_*.p"

            # print(f"Searching for:  {str(Path(basefn, name_start)):s}")
            fng = list(Path(basefn).glob(name_start))
            # print(f"Found: {len(fng):d} files.")
            """ cull list by experiment and db """

            if args.check:
                return
            print("\nConditions: soma= ", PD.soma_inflate, "  dend=", PD.dend_inflate)

            if args.analysis in ["traces", "revcorr"]:
                fng = PS.select_filenames(fng, args)
                times = np.zeros(len(fng))
                for i, f in enumerate(fng):
                    times[i] = f.stat().st_mtime
                # pick most recent file = this should be better managed (master file information)

                ix = np.argmax(times)
                fng = [fng[ix]]

                if ig == 0:
                    P = PH.regular_grid(
                        rows,
                        cols,
                        order="rowsfirst",
                        figsize=(sizex, sizey),
                        panel_labels=plabels,
                        labelposition=(0.05, 0.95),
                        margins={
                            "leftmargin": 0.1,
                            "rightmargin": 0.01,
                            "topmargin": 0.15,
                            "bottommargin": 0.15,
                        },
                    )

                ax = P.axdict[pgbc]
                for k, fn in enumerate(fng):
                    if args.analysis == "traces":
                        PS.plot_traces(ax, fn, PD, args.protocol, nax=k)
                    elif args.analysis == "revcorr":
                        res = PS.compute_revcorr(ax, pgbc, fn, PD, args.protocol)
                        if res is None:
                            return
            elif args.analysis == "singles":
                fna = PS.select_filenames(fng, args)

            elif args.analysis == "tuning":
                P = PS.plot_tuning(args, filename=None, filenames=fng)
    if chan == "_":
        chan = "nav11"
    else:
        chan = chan[:-1]
    P.figure_handle.suptitle(
        f"Model: {modelName:s}  Na Ch: {chan:s} Scaling: {stitle:s} ", fontsize=7
    )
    mpl.show()


def getCommands():
    parser = build_parser()
    args = parser.parse_args()

    # if args.configfile is not None:
    #     config = None
    #     if args.configfile is not None:
    #         if ".json" in args.configfile:
    #             # The escaping of "\t" in the config file is necesarry as
    #             # otherwise Python will try to treat is as the string escape
    #             # sequence for ASCII Horizontal Tab when it encounters it
    #             # during json.load
    #             config = json.load(open(args.configfile))
    #             print(f"Reading JSON configuration file: {args.configfile:s}")
    #         elif ".toml" in args.configfile:
    #             print(f"Reading TOML configuration file: {args.configfile:s}")
    #             with open("wheres_my_data.toml", "r") as fh:
    #                 self.config = toml.load(fh)
    #
    #     vargs = vars(args)  # reach into the dict to change values in namespace
    #
    #     for c in config:
    #         if c in args:
    #             # print("Getting parser variable: ", c)
    #             vargs[c] = config[c]
    #         else:
    #             raise ValueError(
    #                 f"config variable {c:s} does not match with comand parser variables"
    #             )
    #
    #     print("   ... All configuration file variables read OK")
    # now copy into the dataclass
    PD = PData(gradeA=GRPDEF.gradeACells)
    # for key, value in vargs.items():
    #      if key in parnames:
    #          # print('key: ', key)
    #          # print(str(value))
    #          exec(f"PD.{key:s} = {value!r}")
    #      # elif key in runnames:
    #      #     exec(f"runinfo.{key:s} = {value!r}")
    cmdline_display(args, PD)
    return (args, PD)


def main():
    args, PD = getCommands()


if __name__ == "__main__":
    main()