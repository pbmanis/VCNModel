#!/usr/bin/python
"""
Run all gbc IV's as a batch


"""
from __future__ import print_function
import sys
import os
from collections import OrderedDict
import numpy as np
from matplotlib import rc
rc('text', usetex=True)
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import pickle
import matplotlib.pyplot as mpl
import pylibrary.PlotHelpers as PH
import model_run as mrun

default_modelName = 'mGBC'
if len(sys.argv) > 1:
    modelName = sys.argv[1]
else:
    modelName = default_modelName

print('Model Name: {:s}'.format(modelName))
modelType = 'II'

testing = False
forcerun = True

# set up plots
l1 = 0.08
l2 = 0.55
wid = 0.38
ymargin = 0.05
numrows = 6
ht = (1. - 2*(ymargin))/numrows
yp = np.arange(0.05, numrows*ht, ht)
yp = np.flipud(yp)
print ('yp:', yp)
lpos = [0.5, 0.95]
sizer = OrderedDict([('VCNc09nd', {'pos': [l1, wid, yp[0], ht], 'labelpos': lpos, 'ylabel': 'mV'}),
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
])  # dict elements are [left, width, bottom, height] for the axes in the plot.
gr = [(a, a+1, 0, 1) for a in range(0, 8)]   # just generate subplots - shape does not matter
axmap = OrderedDict(zip(sizer.keys(), gr))
P = PH.Plotter(rcshape=sizer, label=False, figsize=(6, 8), labeloffset=[0.6, 0.])
#PH.show_figure_grid(P.figure_handle)

gbc_names = ['09nd', '09', '11', '14', '16', '17', '18', '19', '20', '21', '22']


for n in gbc_names:
    M = mrun.ModelRun() # create an instance
    # create paths to the simulation runs to check for existing IV initialization
    ivinitfile = os.path.join(M.initDirectory, 'IVneuronState_{0:s}_{1:s}.dat'.format(modelName, modelType))
    cell = 'VCN_c{0:s}'.format(n)
    cell_ax = 'VCNc{0:s}'.format(n)
    initf = os.path.join(M.baseDirectory, cell, ivinitfile)
    print ('\nRetrieving data for cell {0:s}'.format(cell))
    M.Params['cell'] = cell
    M.Params['modelType'] = modelType
    M.Params['hocfile'] = M.Params['cell'] + '.hoc'
    if M.Params['hocfile'] == None: # just use the matching hoc file
        M.Params['hocfile'] = M.Params['cell'] + '.hoc'
    if not os.path.isfile(initf):
        print('creating new init file for cell, did not find {:s}'.format(initf))
        M.Params['runProtocol'] = 'initIV'
        print(M.Params)
        if not testing:
            M.run_model(par_map = M.Params)
    else:
        print('    Initialization file for {:s} exists'.format(cell))
    if not os.path.isfile(initf):
        raise ValueError('Failed to create the IV init state file %s' % initf)
    ivdatafile = os.path.join(M.baseDirectory, cell, M.simDirectory, 'IV', '_'.join((cell, 'pulse', modelName, modelType, 'monitor.p')))
    if not os.path.isfile(ivdatafile) or forcerun is True:
        print ('Creating ivdatafile: {:s}\n'.format(ivdatafile))
        M.Params['runProtocol'] = 'runIV'
        if not testing:
            M.run_model(par_map = M.Params)
    else:
        print('    IV file for {:s} exists'.format(cell))
    print(' IV data file: ', ivdatafile)
    fh = open(ivdatafile, 'rb')
    df = pickle.load(fh)
    r = df['Results'][0]

    for trial in range(len(df['Results'])):
        ds = df['Results'][trial]
        k0 = list(ds.keys())[0]
        dx = ds[k0]['monitor']
        P.axdict[cell_ax].plot(dx['time'], dx['postsynapticV'], linewidth=1.0)
        P.axdict[cell_ax].set_xlim(0., 150.)
        P.axdict[cell_ax].set_ylim(-200., 50.)
    PH.calbar(P.axdict[cell_ax], calbar=[120., -95., 25., 20.], axesoff=True, orient='left', 
            unitNames={'x': 'ms', 'y': 'mV'}, font='Arial', fontsize=8)

mpl.savefig('GBC_IVs_{0:s}_{1:s}.pdf'.format(modelName, modelType))
mpl.show()

