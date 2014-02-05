import neuron as h
from neuron import *

import numpy as np
import scipy as sp
import nrnlibrary
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui

class kinetics():
    def __init__(self, args):
        self.app = pg.mkQApp()
        self.win = pg.GraphicsWindow(title="VC Plots")
        self.win.show()
        self.win.resize(800,600)

        self.p1 = self.win.addPlot(title="")
        self.win.nextCol()
        self.p2 = self.win.addPlot(title="")
        self.win.nextRow()
        self.p3 = self.win.addPlot(title="")
        self.win.nextCol()
        self.p4 = self.win.addPlot(title="")

        self.run()
        self.show()

    def show(self):
        QtGui.QApplication.instance().exec_()


    def run(self):
        CaP = nrnlibrary.nrnutils.Mechanism('CaPCalyx')
        leak = nrnlibrary.nrnutils.Mechanism('leak')
        soma = nrnlibrary.nrnutils.Section(L=20, diam=20, mechanisms=[CaP, leak])

        h.celsius = 37 # set the temperature.
        ca_init = 70e-6

        vec={}
        for var in ['time', 'V', 'ICa', 'Vcmd']:
            vec[var] = h.Vector()

        h.dt = 0.025
        v_init = -65.
        tstop = 70.

        clampV = v_init
        vcPost = h.SEClamp(0.5, sec=soma)
        vcPost.dur1 = 5
        vcPost.amp1 = clampV
        vcPost.dur2 = 50
        vcPost.amp2 = clampV-0.0 # just a tiny step to keep the system honest
        vcPost.dur3 = 1.0
        vcPost.amp3 = clampV
        vcPost.rs = 1e-6


        stimamp = np.linspace(-100, 60, num=17, endpoint=True)

        for V in stimamp:
            stim={}
            stim['NP'] = 1
            stim['Sfreq'] = 1 # stimulus frequency
            stim['delay'] = 5
            stim['dur'] = 100
            stim['amp'] = V
            stim['PT'] = 0.0
            vcPost.amp2=V
        #    (secmd, maxt, tstims) = nrnlibrary.makestim.makestim(stim, pulsetype='square', dt = h.dt)
        #    vec['Vcmd'] = h.Vector(secmd)
        #    vec['Vcmd'].play(soma()._ref_v, h.dt , 0, sec=soma)
            vec['ICa'].record(vcPost._ref_i, sec=soma)
            vec['V'].record(soma()._ref_v, sec=soma)
            vec['time'].record(h._ref_t)
            print 'V = ', V
            h.tstop = tstop
            h.finitialize(v_init)
            h.run()
            self.p1.plot(np.array(vec['time']), np.array(vec['ICa']))


if __name__ == "__main__":

    kinetics(sys.argv[1:])
