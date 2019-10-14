#!/usr/bin/python
"""
Run all gbc IV's as a batch

"""
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib import rc
rc('text', usetex=True)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import sys
from pathlib import Path
from collections import OrderedDict
import datetime
import numpy as np

import pickle
import matplotlib.pyplot as mpl
import pylibrary.PlotHelpers as PH
import model_run as mrun
from ephysanalysis import MakeClamps
from ephysanalysis import RmTauAnalysis
from ephysanalysis import SpikeAnalysis

AR = MakeClamps.MakeClamps()
SP = SpikeAnalysis.SpikeAnalysis()
RM = RmTauAnalysis.RmTauAnalysis()

plotflag = True

default_modelName = 'XM13_nacncoop'
#default_modelName = 'XM13'
if len(sys.argv) > 1:
    modelName = sys.argv[1]
else:
    modelName = default_modelName

print('Model Name: {:s}'.format(modelName))
modelType = 'II'


testing = False
forcerun = False

# set up plots
l1 = 0.08
l2 = 0.55
wid = 0.38
ymargin = 0.05
numrows = 6
ht = (1. - 2*(ymargin))/numrows
yp = np.arange(0.05, numrows*ht, ht)
yp = np.flipud(yp)
# print ('yp:', yp)
lpos = [0.5, 0.95]
sizer = OrderedDict([('VCNc08', {'pos': [l1, wid, yp[0], ht], 'labelpos': lpos, 'ylabel': 'mV'}),
                    ('VCNc09',  {'pos': [l1, wid, yp[1], ht], 'labelpos': lpos, 'ylabel': 'mV'}),
                    ('VCNc11',  {'pos': [l1, wid, yp[2], ht], 'labelpos': lpos, 'ylabel': 'mV'}),
                    ('VCNc14',  {'pos': [l1, wid, yp[3], ht], 'labelpos': lpos, 'ylabel': 'mV'}),
                    ('VCNc16',  {'pos': [l1, wid, yp[4], ht], 'labelpos': lpos, 'ylabel': 'mV'}),
                    ('VCNc17',  {'pos': [l1, wid, yp[5], ht], 'labelpos': lpos, 'ylabel': 'mV'}),
                    ('VCNc18',  {'pos': [l2, wid, yp[0], ht], 'labelpos': lpos, 'ylabel': 'mV'}), 
                    ('VCNc19',  {'pos': [l2, wid, yp[1], ht], 'labelpos': lpos, 'ylabel': 'mV'}),
                    ('VCNc20',  {'pos': [l2, wid, yp[2], ht], 'labelpos': lpos, 'ylabel': 'mV'}), 
                    ('VCNc21',  {'pos': [l2, wid, yp[3], ht], 'labelpos': lpos, 'ylabel': 'mV'}),
                    ('VCNc22',  {'pos': [l2, wid, yp[4], ht], 'labelpos': lpos, 'ylabel': 'mV'}),
                    ('None', {}),
])  # dict elements are [left, width, bottom, height] for the axes in the plot.
gbc_names = [s[-2:] for s in sizer.keys()]

gr = [(a, a+1, 0, 1) for a in range(0, 8)]   # just generate subplots - shape does not matter
axmap = OrderedDict(zip(sizer.keys(), gr))
if plotflag:
    P = PH.regular_grid(rows=6, cols=2, figsize=(6, 8), verticalspacing=0.04, panel_labels=list(sizer.keys()),
        margins={'leftmargin': 0.07, 'rightmargin': 0.05, 'topmargin': 0.2, 'bottommargin': 0.1})
#PH.show_figure_grid(P.figure_handle)

print('Models: ', P.axdict.keys())


gbc_names = [s[-2:] for s in sizer.keys()]
# gbc_names = ['16']
gbc_names = ['09', '11', '17', '18']

for n in gbc_names:
    if n == 'ne':
        continue
    M = mrun.ModelRun() # create an instance
    # create paths to the simulation runs to check for existing IV initialization
    ivinitfile = Path(M.initDirectory, f"IVneuronState_{modelName:s}_{modelType:s}.dat")
    print(f'Initfile name: {str(ivinitfile):s}')
    cell = f'VCN_c{n:s}'
    cell_ax = f'VCNc{n:s}'
    M.Params['cell'] = cell
    M.Params['Parallel'] = True
    M.Params['modelType'] = modelType
    M.Params['modelName'] = modelName
    M.Params['soma_autoinflate'] = True
    M.Params['dendrite_autoinflate'] = True
    M.Params['sequence'] = '[-1., 2.01, 0.2]'
    
    if M.Params['hocfile'] == None: # just use the matching hoc file
        M.Params['hocfile'] = M.Params['cell'] + '.hoc'
    M.setup_model(par_map=M.Params)
    sinflate='soma_rawHOC'
    dinflate='dend_rawHOC'
    if M.Params['soma_autoinflate']:
        sinflate = 'soma_scaled'
    if M.Params['dendrite_autoinflate']:
        dinflate='dend_scaled'
    
    d = datetime.datetime.now()
    outfile = f'GBC_IVs_{modelName:s}_{modelType:s}_{sinflate:s}_{dinflate:s}.pdf'
    if plotflag:
        P.figure_handle.suptitle(f"{outfile:s}  {d.strftime('%H:%M:%S %d-%b-%Y'):s}", fontsize=8)
    print(M.Params['hocfile'])
    
    initf = M.Params['initIVStateFile']
    if initf is None:
        continue
    if not initf.is_file() or forcerun:
        print(f'creating new init file for cell, did not find {str(initf):s}')
        M.Params['runProtocol'] = 'initIV'
        if not testing:
            M.run_model(par_map = M.Params)
    else:
        print(f'    Initialization file for {cell:s} exists')
    if not initf.is_file():
        raise ValueError(f'Failed to create the IV init state file {str(initf):s}')

    ivdatafile = M.Params['simulationFilename']
    # ivdatafile = Path(M.baseDirectory, cell, M.simDirectory, 'IV', f"{cell:s}_pulse_{modelName:s}_{modelType:s}_monitor.p")
    print(' file exists: ', ivdatafile.is_file())
    if not ivdatafile.is_file() or forcerun is True:
        print (f'Creating ivdatafile: {str(ivdatafile):s}\n')
        M.Params['runProtocol'] = 'runIV'
        if not testing:
            M.run_model(par_map=M.Params)
    else:
        print('    IV file for {:s} exists'.format(cell))
    print(' IV data file: ', ivdatafile)

# analyse the trace

    # fn = Path('VCN_Cells/VCN_c22/Simulations/IV/VCN_c22_pulse_XM13nacn_II_monitor.p')

    AR.read_pfile(ivdatafile)
    bridge_offset = 0.0
    threshold = -40.
    tgap = 0.  # gap before fittoign taum
        # print(self.AR.tstart, self.AR.tend)
    RM.setup(AR, SP, bridge_offset=bridge_offset)
    SP.setup(clamps=AR, threshold=threshold, 
            refractory=0.0001, peakwidth=0.001, interpolate=True, verify=False, mode='peak')
    SP.analyzeSpikes()
    SP.analyzeSpikeShape()
    # SP.analyzeSpikes_brief(mode='baseline')
    # SP.analyzeSpikes_brief(mode='poststimulus')
    SP.fitOne(function='fitOneOriginal')
    RM.analyze(rmpregion=[0., AR.tstart-0.001],
                tauregion=[AR.tstart, AR.tstart + (AR.tend-AR.tstart)/5.],
                to_peak=True, tgap=tgap)

    RMA = RM.analysis_summary
    print('axis: ', cell_ax)
    if plotflag:
        fh = open(ivdatafile, 'rb')
        df = pickle.load(fh)
        r = df['Results'][0]

        for trial in range(len(df['Results'])):
            ds = df['Results'][trial]
            k0 = list(df['Results'][trial].keys())[0]
            dx = ds[k0]['monitor']
            P.axdict[cell_ax].plot(dx['time'], dx['postsynapticV'], linewidth=1.0)
            P.axdict[cell_ax].set_xlim(0., 150.)
            P.axdict[cell_ax].set_ylim(-200., 50.)
        PH.calbar(P.axdict[cell_ax], calbar=[120., -95., 25., 20.], axesoff=True, orient='left', 
                unitNames={'x': 'ms', 'y': 'mV'}, font='Arial', fontsize=8)
    # P.axdict[cell_ax].title(toptitle, fontsize=7)
if plotflag:
    titletext = f"Model: {modelType:s}  Na Ch: {modelName:s} Scaling: {sinflate:s}_{dinflate:s}" # .replace('_', '\_')
    P.figure_handle.suptitle(titletext, fontsize=9)
    print(str(outfile))
    mpl.savefig(outfile)
    mpl.show()

