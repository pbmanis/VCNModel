__author__ = 'pbmanis'

import pprint
import os, sys
from collections import OrderedDict
import numpy as np
import pickle
import matplotlib
matplotlib.use('Qt5Agg')
rcParams = matplotlib.rcParams
rcParams['svg.fonttype'] = 'none' # No text as paths. Assume font installed.
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as mpl

from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
import pylibrary.plotting.pyqtgraph_plothelpers as pgPH
import pylibrary.plotting.plothelpers as PH  # matplotlib
import vcnmodel.analyzers.analyze_run as AR

pg.setConfigOption('background', 'w')  # set background to white
pg.setConfigOption('foreground', 'k')


class IVPlots():
    def __init__(self, title=None, mode='pg'):
        if mode in ['pg', 'mpl']:
            self.plotmode = mode
        else:
            print('IVPlots: plot mode must be pg or mpl')
            exit()
        if self.plotmode == 'pg':
            self.pg_setup(title=title)
        elif self.plotmode == 'mpl':
            self.mpl_setup(title=title)

    def pg_setup(self, title=None):
        self.app = pg.mkQApp()
        if title is None:
            wintitle = 'NEURON run plots'
        else:
            wintitle = title
        self.view = pg.GraphicsView()
        self.lo = pg.GraphicsLayout()  #(border=(100,100,100))
        self.view.setCentralItem(self.lo)
        self.view.resize(800, 600)
        self.view.show()
        #self.win = pg.GraphicsWindow(title=wintitle)
        #self.win.show()
        #self.win.resize(800,600)
        #self.win.nextRow()
        #self.GL = pg.GraphicsLayoutWidget(parent=self.win)
        self.lo.addLabel(f"Cell: {title:s}", colspan=9, size='12pt')
        self.plots = {}
        for i in range(1, 6):
            self.plots['p%d' % i] = None
        nr1 = 6
        for i in range(1,nr1+1):
            self.lo.addLayout(row=i, col=10)
        for i in range(1,11):
            self.lo.addLayout(row=nr1+2, col = i)
        self.plots['p1'] = self.lo.addPlot(title="Vsoma", row=1, col=1, rowspan=nr1-2, colspan=9)
        self.plots['Dendrite'] = self.lo.addPlot(title="Vdend", row=nr1-1, col=1, rowspan=1, colspan=9)
        self.plots['Iinj'] = self.lo.addPlot(title="IInj", row=nr1, col=1, rowspan=1, colspan=9)
        # self.plots['p4'] = self.lo.addPlot(title="", row=6, col=0, rowspan=1, colspan=1)

    def mpl_setup(self, title):
        if title is None:
            title = "NEURON run plots"
        
        # self.mpl_P = PH.regular_grid(3, 1, order='columns', figsize=(6.0, 4.0), showgrid=False,
        #         verticalspacing=0.08, horizontalspacing=0.08,
        #         margins={'leftmargin': 0.07, 'rightmargin': 0.05, 'topmargin': 0.03, 'bottommargin': 0.1},
        #         labelposition=(0., 0.), parent_figure=None, panel_labels=['Soma', 'Dendrite', 'Iinj'],
        #         title=title)
        # define positions for each panel in Figure coordinages (0, 1, 0, 1)
        # you don't have to use an ordered dict for this, I just prefer it when debugging
        x = 0
        y = 1.05
        sizer = {'Soma': {'pos': [0.1, 0.8, 0.60, 0.35], 'labelpos': (x,y), 'noaxes': True},
                 'Dendrite': {'pos': [0.1, 0.8, 0.20, 0.35], 'labelpos': (x,y), 'noaxes': True},
                 'Iinj': {'pos': [0.1, 0.8, 0.05, 0.10], 'labelpos': (x,y), 'noaxes': True},
                }
        # dict pos elements are [left, width, bottom, height] for the axes in the plot.
        gr = [(a, a+1, 0, 1) for a in range(0, 3)]   # just generate subplots - shape does not matter
        axmap = OrderedDict(zip(sizer.keys(), gr))
        self.mpl_P = PH.Plotter((3, 1), axmap=axmap, label=True, labelalignment='left', figsize=(6., 4.))
        # PH.show_figure_grid(self.mpl_P.figure_handle)
        self.mpl_P.resize(sizer)  # perform positioning magic
        self.plots = self.mpl_P.axdict
        self.plots['p4'] = None


    def show(self):
        if self.plotmode == 'pg':
            QtGui.QApplication.instance().exec_()
        else:
            mpl.show()

    def plotResults(self, res, runInfo, somasite=['postsynapticV', 'postsynapticI', 'dendriteV']):
        clist={'axon': 'r', 'heminode': 'g', 'stalk':'y', 'branch': 'g', 'neck': 'b',
                'swelling': 'm', 'tip': 'k', 'parentaxon': '', r'synapse': 'c', 'Soma': 'k',
                'dendrite': 'c', 'dend': 'c'}
        dx = np.array([x for k,x in res['distanceMap'].items()])
        dlen = res['monitor']['postsynapticV'].shape[0]
        dxmax = np.max(dx)

        
        if self.plotmode == 'pg':
            self.plots['Soma'].setLabel('left', 'V')
            self.plots['Soma'].plot(res['monitor']['time'][:dlen], res['monitor']['postsynapticV'],
                pen=pg.mkPen(clist['Soma'], width=0.5), )
            if 'dendriteV' in somasite and self.plots['Dendrite'] is not None and len(res['monitor']['dendriteV']) > 0:
                self.plots['Dendrite'].plot(res['monitor']['time'][:dlen], 
                    res['monitor']['dendriteV'], pen=pg.mkPen(clist['dendrite'], width=0.5), )
                self.plots['Dendrite'].setLabel('left', 'V (mV)')
            if 'vec' in res.keys()and self.plots['p4'] is not None:
                for v in res['vec']:
                    self.plots['p4'].plot(res['monitor']['time'], res['vec'][v],
                                 pen=pg.mkPen(pg.intColor(int(255.*res['distanceMap'][v]/dxmax))),
                                 width=1.5)
            if 'postsynapticI' in somasite and self.plots['Iinj'] is not None:
                vlen = len(res['monitor']['postsynapticI'])
                self.plots['Iinj'].plot(res['monitor']['time'][0:vlen], res['monitor']['postsynapticI'],
                    pen = pg.mkPen('b', width=0.5))
                self.plots['Iinj'].setLabel('left', 'I (inj, nA)')
            if 'vec' in res.keys() and self.plots['p4'] is not None:
                for v in res['vec']:
                    self.plots['p4'].plot(res['monitor']['time'], res['vec'][v],
                                 pen=pg.mkPen(pg.intColor(int(255.*res['distanceMap'][v]/dxmax))),
                                 width=0.5)
            #for c in res['ICa']:
            #     self.plots['Dendrite'].plot(res['monitor']['time'], res['ICa'][c]*1e12, pen=pg.mkPen('b', width=1.5))
            #     self.plots['Iinj'].plot(res['vec']['time'], res['Cai'][c]*1e6, pen=pg.mkPen('g', width=1.5))

                    #p2.set_ylim(-5e-12, 1e-12)
            self.plots['Soma'].setXRange(0., 120., padding=0.2)
            self.plots['Dendrite'].setXRange(0., 120., padding=0.2)
            if self.plots['Iinj'] is not None:
                self.plots['Iinj'].setXRange(0., 120., padding=0.2)
            
            pgPH.cleanAxes([self.plots['Soma'], self.plots['Dendrite'], self.plots['Iinj'], self.plots['p4']])
            pgPH.nice_plot(self.plots['Soma'])
            pgPH.calbar(self.plots['Soma'], [110, -50, 20, 50])
           # pgPH.nice_plot(self.plots['Dendrite'])
        #    pgPH.calbar(self.plots['Dendrite'], [110, 0.1, 20, 1])
            if self.plots['Iinj'] is not None:
                pgPH.nice_plot(self.plots['Iinj'])
                pgPH.calbar(self.plots['Iinj'], [110, -0.1, 20, 1])


        elif self.plotmode == 'mpl':
            self.plots['Soma'].set_ylabel('V')
            self.plots['Soma'].plot(res['monitor']['time'][:dlen], res['monitor']['postsynapticV'],
            color=clist['Soma'], linewidth=0.5)
            if 'dendriteV' in somasite and self.plots['Dendrite'] is not None and len(res['monitor']['dendriteV']) > 0:
                self.plots['Dendrite'].plot(res['monitor']['time'][:dlen], 
                    res['monitor']['dendriteV'], color=clist['dendrite'], linewidth=0.5)
                self.plots['Dendrite'].set_ylabel('V (mV)')
            if 'vec' in res.keys()and self.plots['p4'] is not None:
                for v in res['vec']:
                    self.plots['p4'].plot(res['monitor']['time'], res['vec'][v],
                                 color=[[int(255.*res['distanceMap'][v]/dxmax)]*3],
                                 linewidth=0.75)
            if 'postsynapticI' in somasite and self.plots['Iinj'] is not None:
                vlen = len(res['monitor']['postsynapticI'])
                self.plots['Iinj'].plot(res['monitor']['time'][0:vlen], res['monitor']['postsynapticI'],
                    color='b', linewidth=0.5)
                self.plots['Iinj'].set_ylabel('I (inj, nA)')
            if 'vec' in res.keys() and self.plots['p4'] is not None:
                for v in res['vec']:
                    self.plots['p4'].plot(res['monitor']['time'], res['vec'][v],
                                 color=[[int(255.*res['distanceMap'][v]/dxmax)]*3],
                                 linewidth=0.75)
            #for c in res['ICa']:
            #     self.plots['Dendrite'].plot(res['monitor']['time'], res['ICa'][c]*1e12, pen=pg.mkPen('b', width=1.5))
            #     self.plots['Iinj'].plot(res['vec']['time'], res['Cai'][c]*1e6, pen=pg.mkPen('g', width=1.5))

                    #p2.set_ylim(-5e-12, 1e-12)
            self.plots['Soma'].set_xlim(0., 120.)
            self.plots['Dendrite'].set_xlim(0., 120.)
            if self.plots['Iinj'] is not None:
                self.plots['Iinj'].set_xlim(0., 120.)
            
            PH.cleanAxes([self.plots['Soma'], self.plots['Dendrite'], self.plots['Iinj']])
            PH.nice_plot(self.plots['Soma'])
            PH.calbar(self.plots['Soma'], [110, -50, 20, 50])
            PH.nice_plot(self.plots['Dendrite'])
            PH.calbar(self.plots['Dendrite'], [110, -50, 20, 50])
           # pgPH.nice_plot(self.plots['Dendrite'])
        #    pgPH.calbar(self.plots['Dendrite'], [110, 0.1, 20, 1])
            if self.plots['Iinj'] is not None:
                PH.nice_plot(self.plots['Iinj'])
                PH.calbar(self.plots['Iinj'], [110, -0.1, 20, 1])        


    def plotFit(self, panel, x, y, c='g'):
        # p = eval("self.plots[\'p%d\' % panel]")
        for j in y.keys():
            if x[j] is None:
                continue
            if self.plotmode == 'pg':
                self.plots[panel].plot(np.array(x[j]), np.array(y[j]), pen=pg.mkPen(c, width=0.5))
            else:
                self.plots[panel].plot(np.array(x[j]), np.array(y[j]), color=c, linewidth=0.35)


if __name__ == "__main__":
    # old stuff probably doesn't work... 
    app = pg.mkQApp()

    cp = IVPlots()
    path = 'Simulations/MNTB_Cell2_cleaned_14.03.09-07.49.58_mrun_ID9999.p'
    path='Simulations/MNTB_Cell2_cleaned_14.03.09-10.43.04_mrun_ID9999.p'
    fo = open(path, 'r')
    d = pickle.load(fo)
    fo.close()
    if isinstance(d['Results'], list):
        for i in range(len(d['Results'])):
            for k in d['Results'][i].keys():
                cp.plotResults(d['Results'][i][k], d['runInfo'])
    else:
        cp.plotResults(d['Results'], d['runinfo'])
    if isinstance(d['Results'], list):
        for i in range(len(d['Results'])):
            arun = AR.AnalyzeRun(d['Results'][i]) # create an instance of the class with the data
            arun.IV()  # compute the IV on the data
            cp.plotFit('Soma', arun.IVResult['taufit'][0], arun.IVResult['taufit'][1], c='r')
            cp.plotFit('Soma', arun.IVResult['ihfit'][0], arun.IVResult['ihfit'][1], c='b')
    cp.show()
