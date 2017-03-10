#!/usr/bin/python
"""

"""
from subprocess import call
import os

import matplotlib.pyplot as mpl
import pylibrary.PlotHelpers as PH
from matplotlib import rc
rc('text', usetex=False)
import pickle
from collections import OrderedDict

forcerun = True  # set tru to force a re-run of the simulation

gbcs = [8, 9, 17, 18, 19, 20, 21, 22]
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
baseDirectory = 'VCN_Cells'
simDirectory = 'Simulations'
nrep = 50
SR = 'MS'
for gbc in gbcs:
    cell = 'VCN_c{0:02d}'.format(gbc)
    an_result_file = 'AN_Result_VCN_c{0:02d}_delays_N{1:03d}_040dB_4000.0_{2:2s}.p'.format(gbc, nrep, SR)
    andatafile = os.path.join(baseDirectory, cell, simDirectory, 'AN', an_result_file)
    print('an result: ', andatafile)
    print (os.path.isfile(andatafile) )
    if not os.path.isfile(andatafile) or forcerun: # only run if no evidence we have run this already 
        if not forcerun:
            call(["python", "all_gbc_AN.py", "%d"%gbc, "%s"%SR, '%d'%nrep])
        else:
            call(["python", "all_gbc_AN.py", "%d"%gbc, "%s"%SR, '%d'%nrep, "forcerun"])
    fh = open(andatafile)
    dx = pickle.load(fh)
    
    for k in dx['somaVoltage'].keys():
        P.axdict[cell].plot(dx['time'], dx['somaVoltage'][k], linewidth=1.0)
    P.axdict[cell].set_xlim(0., 250.)

mpl.show()