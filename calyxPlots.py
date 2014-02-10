__author__ = 'pbmanis'

from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
import pickle
import pylibrary.pyqtgraphPlotHelpers as PH
import pprint
import os, sys
import wx
pg.setConfigOption('background', 'w')  # set background to white
pg.setConfigOption('foreground', 'k')

class CalyxPlots():
    def __init__(self):
        self.app = pg.mkQApp()
        self.win = pg.GraphicsWindow(title="calyx_Plots")
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

    
    def plotResults(self, res, runInfo):
        clist={'axon': 'r', 'heminode': 'g', 'stalk':'y', 'branch': 'g', 'neck': 'b',
                'swelling': 'm', 'tip': 'k', 'parentaxon': '', r'synapse': 'c'}

        self.p1.plot(res['vec']['time'], res['vec']['axon'], pen=pg.mkPen(clist['axon'], width=1.5), )
        for v in res['Voltages']:
            self.p1.plot(res['vec']['time'], res['Voltages'][v],  pen=pg.mkPen('k', width=1.5))
        self.p1.setLabel('left', 'V')
        for c in res['ICa']:
            self.p2.plot(res['vec']['time'], res['ICa'][c]*1e12, pen=pg.mkPen('b', width=1.5))
            self.p3.plot(res['vec']['time'], res['Cai'][c]*1e6, pen=pg.mkPen('g', width=1.5))
        self.p2.setLabel('left', 'I_{Ca}')
        self.p3.setLabel('left', '[Ca]_i')
        self.p4.plot(res['vec']['time'], res['vec']['postsynaptic'], pen = pg.mkPen('b', width=1.5))

                #p2.set_ylim(-5e-12, 1e-12)
     #   PH.cleanAxes([p1, p2, p3, p4])

        # print '='*80
        # print 'runInfo:'
        # pprint.pprint(runInfo)
        # print '='*80
        # fig=MP.figure(101)
        # p21=fig.add_subplot(3,1,1)
        # p22=fig.add_subplot(3,1,2)
        # p23=fig.add_subplot(3,1,3)
        # t = res['vec']['time']
        # tx = len(t)
        # tl = len(res['vec']['i_stim0'])
        # print tl, len(res['vec']['i_stim0'])
        # if tx < tl:
        #     tl = tx
        #
        # p21.plot(t[0:tl], res['vec']['i_stim0'][0:tl])
        # if 'i_stim1' in res['vec'].keys() and len(res['vec']['i_stim1']) == tl:
        #     p22.plot(t[0:tl], res['vec']['i_stim1'])
        # if 'i_stim2' in res['vec'].keys() and len(res['vec']['i_stim2']) == tl:
        #     p23.plot(t[0:tl], res['vec']['i_stim2'])
        # MP.show()

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
