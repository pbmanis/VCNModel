__author__ = 'pbmanis'

from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
import numpy as np
import pickle
import pylibrary.pyqtgraphPlotHelpers as PH
import analyze_run as ar
import pprint
import os, sys
#import wx
pg.setConfigOption('background', 'w')  # set background to white
pg.setConfigOption('foreground', 'k')

class CalyxPlots():
    def __init__(self, title=None):
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
        self.plots = {}
        for i in range(1, 6):
            self.plots['p%d' % i] = None
        nr1 = 6
        for i in range(1,nr1+1):
            self.lo.addLayout(row=i, col=10)
        for i in range(1,11):
            self.lo.addLayout(row=nr1+2, col = i)
        self.plots['p1'] = self.lo.addPlot(title="Vsoma", row=1, col=1, rowspan=nr1-2, colspan=9)
        self.plots['p2'] = self.lo.addPlot(title="Vdend", row=nr1-1, col=1, rowspan=1, colspan=9)
        self.plots['p3'] = self.lo.addPlot(title="IInj", row=nr1, col=1, rowspan=1, colspan=9)
        # self.plots['p4'] = self.lo.addPlot(title="", row=6, col=0, rowspan=1, colspan=1)



    def show(self):
        QtGui.QApplication.instance().exec_()


    def plotResults(self, res, runInfo, somasite=['postsynapticV', 'postsynapticI', 'dendriteV']):
        clist={'axon': 'r', 'heminode': 'g', 'stalk':'y', 'branch': 'g', 'neck': 'b',
                'swelling': 'm', 'tip': 'k', 'parentaxon': '', r'synapse': 'c', 'soma': 'k',
                'dendrite': 'c', 'dend': 'c'}
        dx = np.array([x for k,x in res['distanceMap'].items()])
        self.plots['p1'].setLabel('left', 'V')
        dlen = res['monitor']['postsynapticV'].shape[0]
        self.plots['p1'].plot(res['monitor']['time'][:dlen], res['monitor']['postsynapticV'],
            pen=pg.mkPen(clist['soma'], width=1.5), )
        if 'dendriteV' in somasite and self.plots['p2'] is not None:
            self.plots['p1'].plot(res['monitor']['time'][:dlen], 
                res['monitor']['dendriteV'], pen=pg.mkPen(clist['dendrite'], width=1.5), )
            dxmax = np.max(dx)
            
            self.plots['p2'].setLabel('left', 'V (mV)')
        if 'vec' in res.keys()and self.plots['p4'] is not None:
            for v in res['vec']:
                self.plots['p4'].plot(res['monitor']['time'], res['vec'][v],
                             pen=pg.mkPen(pg.intColor(int(255.*res['distanceMap'][v]/dxmax))),
                             width=2.0)
        if 'postsynapticI' in somasite and self.plots['p3'] is not None:
            vlen = len(res['monitor']['postsynapticI'])
            self.plots['p3'].plot(res['monitor']['time'][0:vlen], res['monitor']['postsynapticI'],
                pen = pg.mkPen('b', width=1.5))
            self.plots['p3'].setLabel('left', 'I (inj, nA)')
        if 'vec' in res.keys() and self.plots['p4'] is not None:
            for v in res['vec']:
                self.plots['p4'].plot(res['monitor']['time'], res['vec'][v],
                             pen=pg.mkPen(pg.intColor(int(255.*res['distanceMap'][v]/dxmax))),
                             width=2.0)
        #for c in res['ICa']:
        #     self.plots['p2'].plot(res['monitor']['time'], res['ICa'][c]*1e12, pen=pg.mkPen('b', width=1.5))
        #     self.plots['p3'].plot(res['vec']['time'], res['Cai'][c]*1e6, pen=pg.mkPen('g', width=1.5))

                #p2.set_ylim(-5e-12, 1e-12)
        self.plots['p1'].setXRange(0., 120., padding=0.2)
        self.plots['p2'].setXRange(0., 120., padding=0.2)
        if self.plots['p3'] is not None:
            self.plots['p3'].setXRange(0., 120., padding=0.2)
            
        PH.cleanAxes([self.plots['p1'], self.plots['p2'], self.plots['p3'], self.plots['p4']])
        PH.nice_plot(self.plots['p1'])
        PH.calbar(self.plots['p1'], [110, -50, 20, 50])
       # PH.nice_plot(self.plots['p2'])
    #    PH.calbar(self.plots['p2'], [110, 0.1, 20, 1])
        if self.plots['p3'] is not None:
            PH.nice_plot(self.plots['p3'])
            PH.calbar(self.plots['p3'], [110, -0.1, 20, 1])
            


    def plotFit(self, panel, x, y, c='g'):
        p = eval("self.plots[\'p%d\' % panel]")
        for j in y.keys():
            p.plot(np.array(x[j]), np.array(y[j]), pen=pg.mkPen(c, width=1.5))


    def show(self):
        QtGui.QApplication.instance().exec_()


# class MyFrame(wx.Frame):
#     def __init__(self, parent, id, title):
#       wx.Frame.__init__(self, parent, id, title)
#
#       self.CreateStatusBar()
#       menuBar = wx.MenuBar()
#       menu = wx.Menu()
#       menu.Append(101, "&File Dialog", "Shows a File Dialog")
#       menuBar.Append(menu, "&Dialogs")
#       self.SetMenuBar(menuBar)
#
#       self.Bind(wx.EVT_MENU, self.openfile, id=101)
#
#     def openfile(self, event):
#         dlg = wx.FileDialog(self, "Choose a file", os.getcwd()+'/Canonical', "", "*.p", wx.OPEN)
#         if dlg.ShowModal() == wx.ID_OK:
#             path = dlg.GetPath()
#             mypath = os.path.basename(path)
#             self.SetStatusText("You selected: %s" % mypath)
#             fo = open(path, 'r')
#             d = pickle.load(fo)
#             fo.close()
#             print d['runInfo'].keys()
#             cp.plotResults(d['Results'], d['runInfo'])
#         dlg.Destroy()
#
# class MyApp(wx.App):
#     def OnInit(self):
#         myframe = MyFrame(None, -1, "calyxPlots")
#         myframe.CenterOnScreen()
#         myframe.Show(True)
#         return True

if __name__ == "__main__":
#    tkapp = MyApp(0)
    app = pg.mkQApp()

    cp = CalyxPlots()
    #path = 'Simulations/Normal14.02.08-16.34.41.p'
    path = 'Simulations/MNTB_Cell2_cleaned_14.03.09-07.49.58_mrun_ID9999.p'
    path='Simulations/MNTB_Cell2_cleaned_14.03.09-10.43.04_mrun_ID9999.p'
    fo = open(path, 'r')
    d = pickle.load(fo)
    fo.close()
    print 'got file'
    if isinstance(d['Results'], list):
#        print 'doing it listway'
        for i in range(len(d['Results'])):
#            print d['Results'][i]
            for k in d['Results'][i].keys():
#                print 'k = ', k
                cp.plotResults(d['Results'][i][k], d['runInfo'])
    else:
        cp.plotResults(d['Results'], d['runinfo'])
    if isinstance(d['Results'], list):
        for i in range(len(d['Results'])):
            arun = ar.AnalyzeRun(d['Results'][i]) # create an instance of the class with the data
            arun.IV()  # compute the IV on the data
            cp.plotFit(1, arun.IVResult['taufit'][0], arun.IVResult['taufit'][1], c='r')
            cp.plotFit(1, arun.IVResult['ihfit'][0], arun.IVResult['ihfit'][1], c='b')
    cp.show()
