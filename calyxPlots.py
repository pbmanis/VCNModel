__author__ = 'pbmanis'

import matplotlib.pylab as MP
import pickle
import pylibrary.PlotHelpers as PH
import pprint
import os, sys
import wx


def plotResults(res, runInfo):
    fig = MP.figure(100)
    p1 = fig.add_subplot(4,1,1)
    p2 = fig.add_subplot(4,1,2)
    p3 = fig.add_subplot(4,1,3)
    p4 = fig.add_subplot(4,1,4)
    p1.plot(res['vec']['time'], res['vec']['axon'], color=runInfo['clist'][0])
    for v in res['Voltages']:
        p1.plot(res['vec']['time'], res['Voltages'][v], color=runInfo['clist'][v])
    p1.set_ylabel('V')
    for c in res['ICa']:
        p2.plot(res['vec']['time'], res['ICa'][c]*1e12, color=runInfo['clist'][c])
        p3.plot(res['vec']['time'], res['Cai'][c]*1e6, color=runInfo['clist'][c])
    p2.set_ylabel('I_{Ca}')
    p3.set_ylabel('[Ca]_i')
    p4.plot(res['vec']['time'], res['vec']['postsynaptic'], color = 'k')
    #p2.set_ylim(-5e-12, 1e-12)
    PH.cleanAxes([p1, p2, p3, p4])

    print '='*80
    print 'runInfo:'
    pprint.pprint(runInfo)
    print '='*80
    fig=MP.figure(101)
    p21=fig.add_subplot(3,1,1)
    p22=fig.add_subplot(3,1,2)
    p23=fig.add_subplot(3,1,3)
    t = res['vec']['time']
    tx = len(t)
    tl = len(res['vec']['i_stim0'])
    print tl, len(res['vec']['i_stim0'])
    if tx < tl:
        tl = tx

    p21.plot(t[0:tl], res['vec']['i_stim0'][0:tl])
    if 'i_stim1' in res['vec'].keys() and len(res['vec']['i_stim1']) == tl:
        p22.plot(t[0:tl], res['vec']['i_stim1'])
    if 'i_stim2' in res['vec'].keys() and len(res['vec']['i_stim2']) == tl:
        p23.plot(t[0:tl], res['vec']['i_stim2'])
    MP.show()

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
        dlg = wx.FileDialog(self, "Choose a file", os.getcwd(), "", "*.p", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            mypath = os.path.basename(path)
            self.SetStatusText("You selected: %s" % mypath)
            fo = open(path, 'r')
            d = pickle.load(fo)
            fo.close()
            print d['runInfo'].keys()
            plotResults(d['Results'], d['runInfo'])
        dlg.Destroy()

class MyApp(wx.App):
    def OnInit(self):
        myframe = MyFrame(None, -1, "calyxPlots")
        myframe.CenterOnScreen()
        myframe.Show(True)
        return True

if __name__ == "__main__":
    app = MyApp(0)
    app.MainLoop()

# if __name__ == "__main__":
#     fo = open('Canonical/Normal.p')
#     d = pickle.load(fo)
#     fo.close()
#     print d['runInfo'].keys()
#     plotResults(d['Results'], d['runInfo'])