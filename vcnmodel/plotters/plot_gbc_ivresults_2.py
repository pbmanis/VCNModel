"""
Version 2: read the new format filenames
"""
import sys
import argparse
import pickle
from pathlib import Path
import datetime
import matplotlib

from ephys.ephysanalysis import MakeClamps
from ephys.ephysanalysis import RmTauAnalysis
from ephys.ephysanalysis import SpikeAnalysis

AR = MakeClamps.MakeClamps()
SP = SpikeAnalysis.SpikeAnalysis()
RM = RmTauAnalysis.RmTauAnalysis()

from matplotlib import rc

rc("text", usetex=False)
import matplotlib.pyplot as mpl
import pylibrary.plotting.plothelpers as PH

modeltypes = ["mGBC", "XM13", "RM03"]
runtypes = ["AN", "an", "IO", "IV", "iv", "gifnoise"]
# modetypes = ['find', 'singles', 'IO', 'multi']


def main():

    parser = argparse.ArgumentParser(description="Plot GBC results")
    parser.add_argument(
        dest="cell",
        action="store",
        nargs="+",
        type=str,
        default=None,
        help="Select the cell(s) (no default)",
    )
    parser.add_argument(
        "-p",
        "--protocol",
        dest="runtype",
        action="store",
        default="IV",
        help=("Select the run type (default: IV) from: %s" % runtypes),
    )
    # parser.add_argument('-m', '--modetype', dest='modetype', action='store',
    #                 default='multi', help=('Select the mode type  (default: multi) from: %s' % modetypes))
    parser.add_argument(
        "-M",
        "--modeltype",
        dest="modeltype",
        action="store",
        default="XM13",
        help=("Select the model type (default XM13) from: %s " % modeltypes),
    )

    parser.add_argument(
        "-s",
        "--scaled",
        dest="scaled",
        action="store_true",
        help=("use scaled data or not"),
    )
    args = parser.parse_args()

    gbc_names = [cn for cn in args.cell]
    if gbc_names[0] == 'A':
        # list the grade A cells
        gbc_names = ['02', '05', '06', '09', '10', '11', '13', '17', '18', '24', '30', '31']

    rows, cols = PH.getLayoutDimensions(len(args.cell))
    P = PH.regular_grid(
        rows,
        cols,
        figsize=(6, 8),
        panel_labels=gbc_names,
        labelposition=(0.05, 0.95),
        margins={
            "leftmargin": 0.07,
            "rightmargin": 0.05,
            "topmargin": 0.12,
            "bottommargin": 0.1,
        },
    )
    chan = "_"  # for certainity in selection, include trailing underscore in this designation

    default_modelName = "XM13_nacncoop"
    # default_modelName = 'XM13'
    if len(sys.argv) > 1:
        modelName = sys.argv[1]
    else:
        modelName = default_modelName
    if len(modelName) > 4:
        chan = modelName[4:] + "_"
        modelName = modelName[:4]
    if args.scaled:
        soma = True
        dend = True
    else:
        soma = False
        dend = False

    stitle = "unknown scaling"

    for gbc in gbc_names:
        basefn = f"/Users/pbmanis/Desktop/Python/VCN-SBEM-Data/VCN_Cells/VCN_c{gbc:2s}/Simulations/IV/"
        # ivdatafile = Path(f"VCN_Cells/VCN_c{gbc:2s}/Simulations/IV/VCN_c{gbc:2s}_pulse_XM13{chan:s}_II_soma=*_monitor.p")
        # print(list(Path(basefn).glob('*.p')))
        # print(basefn)
        print(f"search for:  VCN_c{gbc:2s}_pulse_*.p")
        fng = list(Path(basefn).glob(f"VCN_c{gbc:2s}_pulse_*.p"))
        print("fng: ", fng)
        ivdatafile = None
        print("\nConditions: soma= ", soma, "  dend=", dend)
        for fn in fng:
            fnp = Path(fn)
            print(fnp.is_file())
            if not fnp.is_file():
                print(f"file: {str(fnp):s} NOT FOUND")
                continue
            print("checking file: ", fn)
            with (open(fnp, "rb")) as fh:
                d = pickle.load(fh)
            if "Params" not in list(d.keys()):
                print("File missing Params; need to re-run")
                continue
            # print(d.keys())
            par = d["Params"]
            # print(par)
            if soma and dend:
                if par["soma_inflation"] > 0 and par["dendrite_inflation"] > 0.0:
                    ivdatafile = Path(fn)
                    stitle = "Soma and Dend scaled"
                    print(stitle)
                    break
            elif soma and not dend:
                if par["soma_inflation"] > 0 and par["dendrite_inflation"] < 0.0:
                    ivdatafile = Path(fn)
                    stitle = "Soma only scaled"
                    print(stitle)
                    break
            elif dend and not soma:
                if par["soma_inflation"] < 0 and par["dendrite_inflation"] > 0.0:
                    ivdatafile = Path(fn)
                    stitle = "Dend only scaled"
                    print(stitle)
                    break
            elif not par["soma_autoinflate"] and not par["dendrite_autoinflate"]:
                print("\nConditions x: soma= ", soma, "  dend=", dend)
                ivdatafile = Path(fn)
                stitle = "No scaling (S, D)"
                print(stitle)
                break
            else:
                print("continue")
                continue
        if ivdatafile is None:
            print("No simulation found that matched conditions")
            print(fng)
            continue
        print("datafile to read: ", str(ivdatafile))
        if not ivdatafile.is_file():
            print("no file? : ", str(ivdatafile))
            continue
        # ftime = ivdatafile.stat().st_mtime
        # print('  File time: ', datetime.datetime.fromtimestamp(ftime).strftime('%H:%M:%S %d-%b-%Y'))
        AR.read_pfile(ivdatafile)
        # for i in range(AR.traces.shape[0]):
        #     mpl.plot(AR.time_base, AR.traces[i])
        # mpl.show()
        bridge_offset = 0.0
        threshold = -32.0
        tgap = 0.0  # gap before fittoign taum
        # print(self.AR.tstart, self.AR.tend)
        RM.setup(AR, SP, bridge_offset=bridge_offset)
        SP.setup(
            clamps=AR,
            threshold=threshold,
            refractory=0.0001,
            peakwidth=0.001,
            interpolate=True,
            verify=False,
            mode="peak",
        )
        SP.analyzeSpikes()
        SP.analyzeSpikeShape()
        # SP.analyzeSpikes_brief(mode='baseline')
        # SP.analyzeSpikes_brief(mode='poststimulus')
        SP.fitOne(function="fitOneOriginal")
        RM.analyze(
            rmpregion=[0.0, AR.tstart - 0.001],
            tauregion=[AR.tstart, AR.tstart + (AR.tend - AR.tstart) / 5.0],
            to_peak=True,
            tgap=tgap,
        )

        RMA = RM.analysis_summary
        print(RMA)
        fh = open(ivdatafile, "rb")
        df = pickle.load(fh)
        r = df["Results"][0]

        for trial in range(len(df["Results"])):
            ds = df["Results"][trial]
            k0 = list(df["Results"][trial].keys())[0]
            dx = ds[k0]["monitor"]
            P.axdict[gbc].plot(dx["time"], dx["postsynapticV"], linewidth=1.0)
            P.axdict[gbc].set_xlim(0.0, 150.0)
            P.axdict[gbc].set_ylim(-200.0, 50.0)
        PH.calbar(
            P.axdict[gbc],
            calbar=[120.0, -95.0, 25.0, 20.0],
            axesoff=True,
            orient="left",
            unitNames={"x": "ms", "y": "mV"},
            font="Arial",
            fontsize=8,
        )
        ftname = str(ivdatafile.name)
        ip = ftname.find("_II_") + 4
        ftname = ftname[:ip] + "...\n" + ftname[ip:]
        toptitle = f"{ftname:s}\nRin={RMA['Rin']:.1f} Mohm  Taum={RMA['taum']:.2f} ms"
        P.axdict[gbc].set_title(toptitle, fontsize=5)

    if chan == "_":
        chan = "nav11"
    else:
        chan = chan[:-1]
    P.figure_handle.suptitle(
        f"Model: {modelName:s}  Na Ch: {chan:s} Scaling: {stitle:s} ", fontsize=7
    )
    mpl.show()


if __name__ == "__main__":
    main()
