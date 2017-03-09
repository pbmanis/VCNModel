#!/usr/bin/python
"""
Run all gbc IV's as a batch


"""
from __future__ import print_function
import os
from collections import OrderedDict
import matplotlib.pyplot as mpl
import pylibrary.PlotHelpers as PH
from matplotlib import rc
rc('text', usetex=False)
import pickle

import model_run as mrun

l1 = 0.08
l2 = 0.52
wid = 0.4
ht = 0.18
yp = [0.8, 0.55, 0.3, 0.05]
sizer = OrderedDict([('VCN_c08', [l1, wid, yp[0], ht]), ('VCN_c09', [l1, wid, yp[1], ht]),
        ('VCN_c17', [l1, wid, yp[2], ht]), ('VCN_c18', [l1, wid, yp[3], ht]), 
        ('VCN_c19', [l2, wid, yp[0], ht]), ('VCN_c20', [l2, wid, yp[1], ht]), 
        ('VCN_c21', [l2, wid, yp[2], ht]), ('VCN_c22', [l2, wid, yp[3], ht]),
])  # dict elements are [left, width, bottom, height] for the axes in the plot.
gr = [(a, a+1, 0, 1) for a in range(0, 8)]   # just generate subplots - shape does not matter
axmap = OrderedDict(zip(sizer.keys(), gr))
P = PH.Plotter(rcshape=sizer, label=False, figsize=(6, 8), labeloffset=[0.6, 0.])
#PH.show_figure_grid(P.figure_handle)

M = mrun.ModelRun() # create an instance

gbc_nos = [8, 9, 17, 18, 19, 20, 21, 22]

# create paths to the simulation runs to check for existing IV initialization
ivinitfile = os.path.join(M.initDirectory, 'IVneuronState.dat')

for n in gbc_nos:
    cell = 'VCN_c{0:02d}'.format(n)
    initf = os.path.join(M.baseDirectory, cell, ivinitfile)
    print ('\nRunning Cell {0:s}: '.format(cell))
    M.Params['cell'] = cell
    M.Params['infile'] = M.Params['cell'] + '.hoc'    
    if not os.path.isfile(initf):
        print('creating new init file for cell, did not find {:s}'.format(initf))
        M.Params['runProtocol'] = 'initIV'
#        try:
        M.run_model(par_map = M.Params)
#        except:
#            print('Failed to run model for {:s}'.format(initf))
    else:
        print('    Initialization file for {:s} exists'.format(cell))
    if not os.path.isfile(initf):
        raise ValueError('Failed to create the IV init state file')
    ivdatafile = os.path.join(M.baseDirectory, cell, M.simDirectory, 'IV', cell+'.p')
    if not os.path.isfile(ivdatafile):
        print ('Creating ivdatafile: {:s}\n'.format(ivdatafile))
        M.Params['runProtocol'] = 'runIV'
        M.run_model(par_map = M.Params)
    else:
        print('    IV file for {:s} exists'.format(cell))
    
    fh = open(ivdatafile)
    df = pickle.load(fh)
    r = df['Results'][0]
    for trial in range(len(df['Results'])):
        ds = df['Results'][trial]
        k0 = df['Results'][trial].keys()[0]
        dx = ds[k0]['monitor']
        P.axdict[cell].plot(dx['time'], dx['postsynapticV'], linewidth=1.0)
        P.axdict[cell].set_xlim(0., 150.)

mpl.show()

