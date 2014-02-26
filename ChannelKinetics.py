"""
ChannelKinetics.py
Simple function to display the responses of .mod files to voltage steps in voltage
clamp in a model cell.
This code is primarily for visual verification of model function.
Call: ChannelKinetics modfile
Note: only modfiles that implement voltage-dependent ion channel models make sense to run
with his routine
2/3/2014 PB Manis
"""

import neuron as h
from neuron import *
import gc

import numpy as np
import scipy as sp
import nrnlibrary
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import pylibrary

class ChannelKinetics():
    def __init__(self, args):
        if isinstance(args, list):
            modfile = args[0] # must be string, not list...
        else:
            modfile = args # 'CaPCalyx'
        doKinetics = False
        if len(args) > 1:
            doKinetics = True
        print len(args)
        print doKinetics
        self.app = pg.mkQApp()
        self.win = pg.GraphicsWindow(title="VC Plots")
       # self.win.show()
        self.win.resize(800,600)

        self.p1 = self.win.addPlot(title="I (VC)")
        self.win.nextCol()
        self.p2 = self.win.addPlot(title="I_ss")
        self.win.nextRow()
        self.p3 = self.win.addPlot(title="Vcmd")
        self.win.nextCol()
        self.p4 = self.win.addPlot(title="I_min")
        self.tdur = {'CaPCalyx': [20., 10.], 'hcno': [1000., 200.], 'ih': [1000., 200.], 'ihvcn': [1000., 200.],
                     'nav11': [10., 5.], 'jsrna': [10., 5.], 'kht':[200., 20.], 'klt': [200., 20.], 'nacn': [10., 5.]}

        self.run(modfile=modfile)
        self.win.setWindowTitle('VC Plots: ' + modfile)
        gc.collect()

        if doKinetics:
            self.win2 = pg.GraphicsWindow(title='KineticPlots')
            self.win2.resize(800, 600)
            self.kp1 = self.win.adPlot(title="htau")
            self.computeKinetics('nav11')

        self.show()

    def show(self):
        self.win.show()
        QtGui.QApplication.instance().exec_()


    def run(self, modfile='CaPCalyx'):
        if isinstance(modfile, list):
            modfile = modfile[0]
        if modfile in self.tdur:
            tstep = self.tdur[modfile]
        else:
            tstep = [200., 50.]
        tdelay = 5.0
        Channel = nrnlibrary.nrnutils.Mechanism(modfile)
        leak = nrnlibrary.nrnutils.Mechanism('leak')
        Channel.set_parameters({'gbar': 1})
        leak.set_parameters({'gbar': 1e-12})

        self.soma = nrnlibrary.nrnutils.Section(L=10, diam=10, mechanisms=[Channel, leak])
#        Channel.insert_into(soma)
#        leak.insert_into(soma)
        print dir(self.soma)
        h.celsius = 22 # set the temperature.
        ca_init = 70e-6

        self.vec={}
        for var in ['time', 'V', 'IChan', 'Vcmd']:
            self.vec[var] = h.Vector()

        h.dt = 0.025
        v_init = -65.
        h.celsius=22.

        clampV = v_init
        self.vcPost = h.SEClamp(0.5, sec=self.soma)
        self.vcPost.dur1 = tdelay
        self.vcPost.amp1 = clampV
        self.vcPost.dur2 = tstep[0]
        self.vcPost.amp2 = clampV-0.0 # just a tiny step to keep the system honest
        self.vcPost.dur3 = tstep[1]
        self.vcPost.amp3 = clampV
        self.vcPost.rs = 1e-6
        print "soma: ", self.soma
        print 'vcpost sec: ', self.vcPost.Section()

        stimamp = np.linspace(-100, 60, num=35, endpoint=True)
        self.ivss = np.zeros((2, stimamp.shape[0]))
        self.ivmin = np.zeros((2, stimamp.shape[0]))

        for i, V in enumerate(stimamp):
            stim={}
            stim['NP'] = 1
            stim['Sfreq'] = 1 # stimulus frequency
            stim['delay'] = 5
            stim['dur'] = 100
            stim['amp'] = V
            stim['PT'] = 0.0
            self.vcPost.amp2=V
            self.vec['IChan'].record(self.vcPost._ref_i, sec=self.soma)
            self.vec['V'].record(self.soma()._ref_v, sec=self.soma)
            self.vec['time'].record(h._ref_t)
            print 'V = ', V
            h.tstop = self.vcPost.dur1+self.vcPost.dur2+self.vcPost.dur3
            h.finitialize(v_init)
            h.run()
            self.t = np.array(self.vec['time'])
            self.ichan = np.array(self.vec['IChan'])
            self.v = np.array(self.vec['V'])
            self.p1.plot(self.t, self.ichan)
            self.p3.plot(self.t, self.v)
            (self.ivss[1,i], r2) = pylibrary.Utility.measure('mean', self.t, self.ichan, tdelay+tstep[0]-10., tdelay+tstep[0])
            (self.ivmin[1,i], r2) = pylibrary.Utility.measure('minormax', self.t, self.ichan, tdelay+0.1, tdelay+tstep[0]/5.0)
            self.ivss[0,i] = V
            self.ivmin[0,i] = V

        self.p2.plot(self.ivss[0,:], self.ivss[1,:], symbol='o', symbolsize=2.0, pen= pg.mkPen('r'))
        self.p4.plot(self.ivmin[0,:], self.ivmin[1,:], symbol='s', symbolsize=2.0, pen=pg.mkPen('b'))


    def computeKinetics(self, ch):
        pass

if __name__ == "__main__":

    ChannelKinetics(sys.argv[1:])
