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
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as mpl
import pylibrary.plotting.plothelpers as PH
from matplotlib import rc
rc('text', usetex=False)
import pickle
from collections import OrderedDict
import all_gbc_AN as oneAN

default_modelName = 'XM13nacncoop'
if len(sys.argv) > 1:
    modelName = sys.argv[1]
else:
    
    modelName = default_modelName

import model_run as mrun
M = mrun.ModelRun() # create an instance  (but do NOT run simulations - just to get setup information)

runs = True

def setuppars(cell, modelType, modelName, protocol, nreps, inflateflag, inputpattern=None):
    M.Params['cell'] = cell
    M.Params['Parallel'] = False
    M.Params['modelType'] = modelType
    M.Params['modelName'] = modelName
    M.Params['runProtocol'] = protocol
    M.Params['soma_autoinflate'] = inflateflag
    M.Params['dendrite_autoinflate'] = inflateflag
    M.Params['inputPattern'] = inputpattern
    
    M.Params['SGCmodelType'] = 'cochlea'
    M.Params['soundtype'] = 'tonepip'
    M.Params['SR'] = SR
    M.Params['SRType'] = SR
    M.Params['ANSynapseType'] = 'multisite'
    M.Params['run_duration'] = 0.3
    M.Params['plotFlag'] = False
    M.Params['nReps'] = nrep
    
    M.Params['hocfile'] = M.Params['cell'] + '.hoc'
    print(f"{' ':14s}setup pars HOC file: {str(M.Params['hocfile']):s}")
    M.setup_model(par_map=M.Params)
    return(M)

print("all_gbc_names_AN")
print('Model Name: {:s}'.format(modelName))
modelType = 'II'
protocol = 'runANPSTH'
require_initialize = True
inflateflag = True
testing = False
forcerun = False  # set true to force a re-run of the simulation
#PH.show_figure_grid(P.figure_handle)
baseDirectory = 'VCN_Cells'
simDirectory = 'Simulations'
nrep = 50
SR = 'MS'
testflag = False

sim_reportfile = Path('lastsim.txt')

# gbc_names = ['08', '09', '09nd', '09h', '17', '18', '19', '20', '21', '22']
l1 = 0.08
l2 = 0.52
wid = 0.4
ht = 0.15
yp = [0.85, 0.65, 0.45, 0.25, 0.05, 0.0]
sizer = OrderedDict([('VCN_c08', {'pos': [l1, wid, yp[0], ht]}),
                     ('VCN_c09', {'pos': [l1, wid, yp[1], ht]}),
                     # ('VCN_c09h', [{'pos': l1, wid, yp[1], ht]}),
                     # ('VCN_c09nd', [l1, wid, yp[2], ht]}),
                     ('VCN_c11', {'pos': [l1, wid, yp[2], ht]}),
                     ('VCN_c14', {'pos': [l1, wid, yp[3], ht]}),
                     ('VCN_c16', {'pos': [l1, wid, yp[4], ht]}),
                     ('VCN_c17', {'pos': [l2, wid, yp[5], ht]}),
                     ('VCN_c18', {'pos': [l2, wid, yp[0], ht]}),
                     ('VCN_c19', {'pos': [l2, wid, yp[1], ht]}),
                     ('VCN_c20', {'pos': [l2, wid, yp[3], ht]}), 
                     ('VCN_c21', {'pos': [l1, wid, yp[4], ht]}),
                     ('VCN_c22', {'pos': [l2, wid, yp[5], ht]}),
                     ('None', {}),
])  # dict elements are [left, width, bottom, height] for the axes in the plot.
# gr = [(a, a+1, 0, 1) for a in range(0, 8)]   # just generate subplots - shape does not matter

gbc_names = [s[-2:] for s in sizer.keys()]
gbc_names = ['09']

P = PH.regular_grid(6, 2, figsize=(6,8), panel_labels=list(sizer.keys()), labelposition=(0.05, 0.95))
# P = PH.Plotter(rcshape=sizer, label=False, figsize=(6, 8), labeloffset=[0.6, 0.])

inputPattern = "VCN_c10" # None # 'VCN_c10'

for gbc in gbc_names:
    if gbc == 'ne':
        continue
    print ('='*32)
    cell = 'VCN_c{0:s}'.format(gbc)
    
    # an_result_file = 'AN_Result_VCN_c{0:s}_delays_N{1:03d}_040dB_4000.0_{2:2s}.p'.format(gbc, nrep, SR)
    # andatafile = Path(baseDirectory, cell, simDirectory, 'AN', an_result_file)
    # print('  an result file: {0:s}'.format(str(andatafile)))
    # print (andatafile.is_file() )
    # if not andatafile.is_file() or forcerun: # only run if no evidence we have run this already
    
    # ANR = oneAN.OneANRun(gbc, modelName, modelType, protocol, SR, nrep, forcerun=forcerun, testing=testing, inflateflag=inflateflag)
    runpars = {'gbc': gbc, 'modelName': modelName, 'modelType': modelType, 'protocol': protocol,
    'SR': SR, 'nrep': nrep, 'initialize': require_initialize, 'forcerun': forcerun, 'testing': testing, 
    'inflateflag': inflateflag, 'inputPattern': inputPattern}
    with open('oneanrun.p', 'wb') as fh:
        pickle.dump(runpars, fh)

    if runs:
        call(['python', 'all_gbc_AN.py'])  # check initialization
        runpars['initialize'] = False
        runpars['forcerun'] = True
        with open('oneanrun.p', 'wb') as fh:
            pickle.dump(runpars, fh)
        print('runpars gbc: ', runpars['gbc'])
        call(['python', 'all_gbc_AN.py'])
        print('call done')
    
    pars = setuppars(cell, modelType, modelName, protocol, nrep, inflateflag, inputpattern=inputPattern)
    sfile = pars.Params['simulationFilename']  # name from that simulation
    lastsim = None
    print(f"*** looking for simulation file: {str(sfile):s}")
    if sfile.is_file():
        print(f'*** found simulation file: {str(sfile):s}')
        lastsim = sfile
    else:
        if lastsim is None and sim_reportfile.is_file():
            lastsim = sim_reportfile.read_text()
        # with open(f"lastsim.txt", 'r') as fh:
        #     print('fh: ', fh)
        #     lastsim = fh.read()
        else:
            raise FileNotFoundError('  Failed to find result file {0:s}'.format(str(sim_reportfile)))
    print('lastsim: ', lastsim)
    with open(lastsim, 'rb') as fh:
        dx = pickle.load(fh)
    for trial in list(dx['trials'].keys()):
        # print(dx['trials'][trial]['somaVoltage'])
        P.axdict[cell].plot(dx['trials'][trial]['time'], dx['trials'][trial]['somaVoltage'], linewidth=1.0)
    P.axdict[cell].set_xlim(0., 250.)
sinflate='soma_rawHOC'
dinflate='dend_rawHOC'
if M.Params['soma_autoinflate']:
    sinflate = 'soma_scaled'
if M.Params['dendrite_autoinflate']:
    dinflate='dend_scaled'
d = datetime.datetime.now()
outfile = f'GBC_IVs_{modelName:s}_{modelType:s}_{sinflate:s}_{dinflate:s}.pdf'
P.figure_handle.suptitle(f"{outfile:s}  {d.strftime('%H:%M:%S %d-%b-%Y'):s}", fontsize=8)


mpl.savefig(outfile)
mpl.show()