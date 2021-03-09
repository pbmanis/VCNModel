#!/usr/bin/python
"""

"""
from __future__ import print_function
from subprocess import call
import os
import sys
import datetime
from pathlib import Path
import matplotlib

matplotlib.use("Qt4Agg")
import matplotlib.pyplot as mpl
import pylibrary.plotting.plothelpers as PH
from matplotlib import rc

rc("text", usetex=False)
import pickle
from collections import OrderedDict
import all_gbc_AN as oneAN

default_modelName = "XM13_nacncoop"
if len(sys.argv) > 1:
    modelName = sys.argv[1]
else:
    modelName = default_modelName
plotflag = True
import src.vcnmodel.model_run as mrun

M = (
    mrun.ModelRun()
)  # create an instance  (but do NOT run simulations - just to get setup information)

runs = True
gradeA = [2, 5, 6, 9, 10, 11, 13, 17, 24, 29, 30]


def setuppars(
    cell, modelType, modelName, protocol, nreps, inflateflag, inputpattern=None
):
    M.Params["cell"] = cell
    M.Params["Parallel"] = False
    M.Params["modelType"] = modelType
    M.Params["modelName"] = modelName
    M.Params['hocfile'] = cell+"_Full.hoc"
    M.Params["runProtocol"] = protocol
    M.Params["soma_autoinflate"] = inflateflag
    M.Params["dendrite_autoinflate"] = inflateflag
    M.Params["inputPattern"] = inputpattern

    M.Params["SGCmodelType"] = "cochlea"
    M.Params["soundtype"] = "tonepip"
    M.Params["SR"] = SR
    M.Params["SRType"] = SR
    M.Params["ANSynapseType"] = "multisite"
    M.Params["run_duration"] = 0.3
    M.Params["plotFlag"] = False
    M.Params["nReps"] = nrep

    # M.Params["hocfile"] = M.Params["cell"] + ".hoc"
    print(f"{' ':14s}setup pars HOC file: {str(M.Params['hocfile']):s}")
    M.setup_model(par_map=M.Params)
    return M


print("all_gbc_names_AN")
print("Model Name: {:s}".format(modelName))
modelType = "II"
protocol = "runANPSTH"
require_initialize = True
inflateflag = True
testing = False
forcerun = False  # set true to force a re-run of the simulation
# PH.show_figure_grid(P.figure_handle)
baseDirectory = "VCN_Cells"
simDirectory = "Simulations"
nrep = 50
SR = "MS"
testflag = False
inputPattern = None

sim_reportfile = Path("lastsim.txt")

cellnames = [f"VCN_c{k:02d}" for k in gradeA]


(r, c) = PH.getLayoutDimensions(len(cellnames), pref="height")
if r * c > len(cellnames):
    for i in range(len(cellnames), r * c):
        cellnames.append("ne")
print("cellnames: ", cellnames)
if plotflag:
    P = PH.regular_grid(
        rows=r,
        cols=c,
        order="columnsfirst",
        figsize=(8, 10),
        verticalspacing=0.04,
        panel_labels=cellnames,
        margins={
            "leftmargin": 0.07,
            "rightmargin": 0.05,
            "topmargin": 0.2,
            "bottommargin": 0.1,
        },
    )
# PH.show_figure_grid(P.figure_handle)

print("Models: ", P.axdict.keys())


for gbc in cellnames:
    if gbc == "ne" or "":
        continue
    M = mrun.ModelRun()  # create an instance
    if gbc == "ne":
        continue
    print("=" * 32)
    cell = gbc

    # an_result_file = 'AN_Result_VCN_c{0:s}_delays_N{1:03d}_040dB_4000.0_{2:2s}.p'.format(gbc, nrep, SR)
    # andatafile = Path(baseDirectory, cell, simDirectory, 'AN', an_result_file)
    # print('  an result file: {0:s}'.format(str(andatafile)))
    # print (andatafile.is_file() )
    # if not andatafile.is_file() or forcerun: # only run if no evidence we have run this already

    # ANR = oneAN.OneANRun(gbc, modelName, modelType, protocol, SR, nrep, forcerun=forcerun, testing=testing, inflateflag=inflateflag)
    runpars = {
        "gbc": gbc,
        "modelName": modelName,
        "modelType": modelType,
        "hocfile": cell + "_Full.hoc",
        "protocol": protocol,
        "SR": SR,
        "nrep": nrep,
        "initialize": require_initialize,
        "forcerun": forcerun,
        "testing": testing,
        "inflateflag": inflateflag,
        "inputPattern": inputPattern,
    }
    with open("oneanrun.p", "wb") as fh:
        pickle.dump(runpars, fh)

    print('run preped:::')
    if runs:
        call(["python", "vcnmodel/simulators/all_gbc_AN.py"])  # check initialization
        runpars["initialize"] = False
        runpars["forcerun"] = True
        with open("oneanrun.p", "wb") as fh:
            pickle.dump(runpars, fh)
        print("runpars gbc: ", runpars["gbc"])
        call(["python", "vcnmodel/simulators/all_gbc_AN.py"])
        print("call done")
    print('Run Done?')
    pars = setuppars(
        cell,
        modelType,
        modelName,
        protocol,
        nrep,
        inflateflag,
        inputpattern=inputPattern,
    )
    sfile = pars.Params["simulationFilename"]  # name from that simulation
    lastsim = None
    print(f"*** looking for simulation file: {str(sfile):s}")
    if sfile.is_file():
        print(f"*** found simulation file: {str(sfile):s}")
        lastsim = sfile
    else:
        if lastsim is None and sim_reportfile.is_file():
            lastsim = sim_reportfile.read_text()
        # with open(f"lastsim.txt", 'r') as fh:
        #     print('fh: ', fh)
        #     lastsim = fh.read()
        else:
            raise FileNotFoundError(
                "  Failed to find result file {0:s}".format(str(sim_reportfile))
            )
    print("lastsim: ", lastsim)
    with open(lastsim, "rb") as fh:
        dx = pickle.load(fh)
    for trial in list(dx["trials"].keys()):
        # print(dx['trials'][trial]['somaVoltage'])
        P.axdict[cell].plot(
            dx["trials"][trial]["time"],
            dx["trials"][trial]["somaVoltage"],
            linewidth=1.0,
        )
    P.axdict[cell].set_xlim(0.0, 250.0)
sinflate = "soma_rawHOC"
dinflate = "dend_rawHOC"
if M.Params["soma_autoinflate"]:
    sinflate = "soma_scaled"
if M.Params["dendrite_autoinflate"]:
    dinflate = "dend_scaled"
d = datetime.datetime.now()
outfile = f"GBC_IVs_{modelName:s}_{modelType:s}_{sinflate:s}_{dinflate:s}.pdf"
P.figure_handle.suptitle(
    f"{outfile:s}  {d.strftime('%H:%M:%S %d-%b-%Y'):s}", fontsize=8
)


mpl.savefig(outfile)
mpl.show()
