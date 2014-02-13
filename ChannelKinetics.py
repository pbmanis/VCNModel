"""
ChannelKinetics.py
Simple function to display the responses of .mod files to voltage steps in voltage
clamp in a model cell.
This code is primarily for visual verification of model function.
2/3/2014 PB Manis
"""

import neuron as h
from neuron import *

import numpy as np
import scipy as sp
import nrnlibrary
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import pylibrary

class ChannelKinetics():
    def __init__(self, args):
        self.app = pg.mkQApp()
        self.win = pg.GraphicsWindow(title="VC Plots")
        self.win.show()
        self.win.resize(800,600)

        self.p1 = self.win.addPlot(title="I (VC)")
        self.win.nextCol()
        self.p2 = self.win.addPlot(title="I_ss")
        self.win.nextRow()
        self.p3 = self.win.addPlot(title="Vcmd")
        self.win.nextCol()
        self.p4 = self.win.addPlot(title="I_min")
        if isinstance(args, list):
            modfile = args[0] # must be string, not list...
        else:
            modfile = args # 'CaPCalyx'
        self.run(modfile=modfile)
        self.win.setWindowTitle('VC Plots: ' + modfile)
        self.show()

    def show(self):
        QtGui.QApplication.instance().exec_()


    def run(self, modfile='CaPCalyx'):
        if isinstance(modfile, list):
            modfile = modfile[0]
        Channel = nrnlibrary.nrnutils.Mechanism(modfile)
        leak = nrnlibrary.nrnutils.Mechanism('leak')
        #leak.set_paramters({'gbar': 1e-6})

        soma = nrnlibrary.nrnutils.Section(L=20, diam=20, mechanisms=[Channel, leak])
#        Channel.insert_into(soma)
#        leak.insert_into(soma)

        h.celsius = 22 # set the temperature.
        ca_init = 70e-6

        vec={}
        for var in ['time', 'V', 'ICa', 'Vcmd']:
            vec[var] = h.Vector()

        h.dt = 0.025
        v_init = -65.
        tstop = 255.

        clampV = v_init
        vcPost = h.SEClamp(0.5, sec=soma)
        vcPost.dur1 = 5
        vcPost.amp1 = clampV
        vcPost.dur2 = 200.
        vcPost.amp2 = clampV-0.0 # just a tiny step to keep the system honest
        vcPost.dur3 = 50.0
        vcPost.amp3 = clampV
        vcPost.rs = 1e-9
        print "soma: ", soma
        print 'vcpost sec: ', vcPost.Section()

        stimamp = np.linspace(-100, 60, num=35, endpoint=True)
        ivss = np.zeros((2, stimamp.shape[0]))
        ivmin = np.zeros((2, stimamp.shape[0]))

        for i, V in enumerate(stimamp):
            stim={}
            stim['NP'] = 1
            stim['Sfreq'] = 1 # stimulus frequency
            stim['delay'] = 5
            stim['dur'] = 100
            stim['amp'] = V
            stim['PT'] = 0.0
            vcPost.amp2=V
            vec['ICa'].record(vcPost._ref_i, sec=soma)
            vec['V'].record(soma()._ref_v, sec=soma)
            vec['time'].record(h._ref_t)
            print 'V = ', V
            h.tstop = tstop
            h.finitialize(v_init)
            h.run()
            t = np.array(vec['time'])
            ica = np.array(vec['ICa'])
            self.p1.plot(t, ica)
            self.p3.plot(t, np.array(vec['V']))
            (ivss[1,i], r2) = pylibrary.Utility.measure('mean', t, ica, 190., 200.)
            (ivmin[1,i], r2) = pylibrary.Utility.measure('min', t, ica, 5.2, 8.)
            ivss[0,i] = V
            ivmin[0,i] = V

        self.p2.plot(ivss[0,:], ivss[1,:], symbol='o', symbolsize=2.0, pen= pg.mkPen('r'))
        self.p4.plot(ivmin[0,:], ivmin[1,:], symbol='s', symbolsize=2.0, pen=pg.mkPen('b'))

if __name__ == "__main__":

    ChannelKinetics(sys.argv[1:])
