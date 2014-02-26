__author__ = 'pbmanis'

from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
import numpy as np
import pickle
import pylibrary.pyqtgraphPlotHelpers as PH
import pprint
import os, sys
import wx
pg.setConfigOption('background', 'w')  # set background to white
pg.setConfigOption('foreground', 'k')

class CalyxPlots():
    def __init__(self, title=None):
        self.app = pg.mkQApp()
        if title is None:
            wintitle = 'NEURON run plots'
        else:
            wintitle = title

        self.win = pg.GraphicsWindow(title=wintitle)
        self.win.show()
        self.win.resize(800,600)

        self.p1 = self.win.addPlot(title="")
        self.win.nextRow()
        self.p2 = self.win.addPlot(title="")
        self.win.nextRow()
        self.p3 = self.win.addPlot(title="")
        self.win.nextRow()
        self.p4 = self.win.addPlot(title="")


    def show(self):
        QtGui.QApplication.instance().exec_()


    
    def plotResults(self, res, runInfo, somasite=['postsynapticV', 'postsynapticI']):
        clist={'axon': 'r', 'heminode': 'g', 'stalk':'y', 'branch': 'g', 'neck': 'b',
                'swelling': 'm', 'tip': 'k', 'parentaxon': '', r'synapse': 'c', 'soma': 'k'}
        dx = np.array([x for k,x in res['distanceMap'].items()])
        self.p1.plot(res['monitor']['time'], res['monitor'][somasite[0]], pen=pg.mkPen(clist['soma'], width=1.5), )
        dxmax = np.max(dx)
        if 'vec' in res.keys():
            for v in res['vec']:
                self.p3.plot(res['monitor']['time'], res['vec'][v],
                             pen=pg.mkPen(pg.intColor(int(255.*res['distanceMap'][v]/dxmax))),
                             width=2.0)
        self.p1.setLabel('left', 'V')
        vlen = len(res['monitor'][somasite[1]])
        self.p2.plot(res['monitor']['time'][0:vlen], res['monitor'][somasite[1]], pen = pg.mkPen('b', width=1.5))
        tlen = len(res['monitor']['time'])
        #self.p3.plot(res['monitor']['time'], res['monitor']['cmd'][0:tlen])
        #for c in res['ICa']:
        #     self.p2.plot(res['monitor']['time'], res['ICa'][c]*1e12, pen=pg.mkPen('b', width=1.5))
        #     self.p3.plot(res['vec']['time'], res['Cai'][c]*1e6, pen=pg.mkPen('g', width=1.5))
        self.p2.setLabel('left', 'I (inj, nA)')
        self.p3.setLabel('left', 'V (mV)')

                #p2.set_ylim(-5e-12, 1e-12)
        PH.cleanAxes([self.p1, self.p2, self.p3, self.p4])


    def show(self):
        QtGui.QApplication.instance().exec_()


class MyFrame(wx.Frame):
    def __init__(self, parent, id, title):
      wx.Frame.__init__(self, parent, id, title)

      self.CreateStatusBar()
      menuBar = wx.MenuBar()
      menu = wx.Menu()
      menu.Append(101, "&File Dialog", "Shows a File Dialog")
      menuBar.Append(menu, "&Dialogs")
      self.SetMenuBar(menuBar)

      self.Bind(wx.EVT_MENU, self.openfile, id=101)

    def openfile(self, event):
        dlg = wx.FileDialog(self, "Choose a file", os.getcwd()+'/Canonical', "", "*.p", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            mypath = os.path.basename(path)
            self.SetStatusText("You selected: %s" % mypath)
            fo = open(path, 'r')
            d = pickle.load(fo)
            fo.close()
            print d['runInfo'].keys()
            cp.plotResults(d['Results'], d['runInfo'])
        dlg.Destroy()

class MyApp(wx.App):
    def OnInit(self):
        myframe = MyFrame(None, -1, "calyxPlots")
        myframe.CenterOnScreen()
        myframe.Show(True)
        return True

if __name__ == "__main__":
#    tkapp = MyApp(0)
    app = pg.mkQApp()

    cp = CalyxPlots()
    path = 'Simulations/Normal14.02.08-16.34.41.p'
    fo = open(path, 'r')
    d = pickle.load(fo)
    fo.close()
    print 'got file'
    print d.keys()
    cp.plotResults(d['Results'], d['runInfo'])
    cp.show()
