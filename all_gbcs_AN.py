#!/usr/bin/python
"""

"""
from __future__ import print_function
from subprocess import call
import os

import matplotlib.pyplot as mpl
import pylibrary.PlotHelpers as PH
from matplotlib import rc
rc('text', usetex=False)
import pickle
from collections import OrderedDict

forcerun = False  # set true to force a re-run of the simulation

gbc_names = ['08', '09', '09nd', '09h', '17', '18', '19', '20', '21', '22']
l1 = 0.08
l2 = 0.52
wid = 0.4
ht = 0.15
yp = [0.85, 0.65, 0.45, 0.25, 0.05]
sizer = OrderedDict([('VCN_c09', [l1, wid, yp[0], ht]), ('VCN_c17', [l2, wid, yp[0], ht]),
                     ('VCN_c09h', [l1, wid, yp[1], ht]), ('VCN_c18', [l2, wid, yp[1], ht]),
                     ('VCN_c09nd', [l1, wid, yp[2], ht]), ('VCN_c19', [l2, wid, yp[2], ht]), 
                     ('VCN_c08', [l1, wid, yp[3], ht]), ('VCN_c20', [l2, wid, yp[3], ht]), 
                     ('VCN_c21', [l1, wid, yp[4], ht]), ('VCN_c22', [l2, wid, yp[4], ht]),
])  # dict elements are [left, width, bottom, height] for the axes in the plot.
gr = [(a, a+1, 0, 1) for a in range(0, 8)]   # just generate subplots - shape does not matter
axmap = OrderedDict(zip(sizer.keys(), gr))
P = PH.Plotter(rcshape=sizer, label=False, figsize=(6, 8), labeloffset=[0.6, 0.])
#PH.show_figure_grid(P.figure_handle)
baseDirectory = 'VCN_Cells'
simDirectory = 'Simulations'
nrep = 50
SR = 'MS'
testflag = False
print("all_gbc_names_AN")
for gbc in gbc_names:
    print()
    print ('='*32)
    cell = 'VCN_c{0:s}'.format(gbc)
    an_result_file = 'AN_Result_VCN_c{0:s}_delays_N{1:03d}_040dB_4000.0_{2:2s}.p'.format(gbc, nrep, SR)
    andatafile = os.path.join(baseDirectory, cell, simDirectory, 'AN', an_result_file)
    print('  an result file: {0:s}'.format(andatafile))
    print (os.path.isfile(andatafile) )
    if not os.path.isfile(andatafile) or forcerun: # only run if no evidence we have run this already 
        cmdlist = ["python", "all_gbc_AN.py", "%s" % gbc, "%s"%SR, '%d'%nrep]
        if testflag:
            cmdlist.append('test')
        if forcerun:
            cmdlist.append('forcerun')
        print(cmdlist)
        call(cmdlist)
    try:
        with open(andatafile, 'r') as fh:
            dx = pickle.load(fh)
            for k in dx['somaVoltage'].keys():
                P.axdict[cell].plot(dx['time'], dx['somaVoltage'][k], linewidth=1.0)
            P.axdict[cell].set_xlim(0., 250.)
    except IOError:
            print('  Failed to find result file {0:s}'.format(andatafile))
mpl.show()