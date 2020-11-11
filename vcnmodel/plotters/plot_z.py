import numpy as np
import pickle
import pyqtgraph as pg
import matplotlib.pyplot as mpl
import pylibrary.plotting.plothelpers as PH
import pylibrary.plotting.styler as STY

# import seaborn as sns

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stixsans'
# matplotlib.rcParams['font.family'] = 'sans-serif'
# matplotlib.rcParams['font.sans-serif'] = ['stixsans'] #, 'Tahoma', 'DejaVu Sans',
#                                #'Lucida Grande', 'Verdana']
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['text.usetex'] = False


pg = False

f = ['VCN_c09_Full_Z.pkl', 'VCN_c09_NoUninnervated_Z.pkl', 'VCN_c09_NoDend_Z.pkl']
fi = [17] # [2, 5, 6, 9, 10, 11, 13, 17, 30]
#fi = [9, 11, 13, 30]
# fi = [2, 5, 6, 10, 17]
f = []
for fin in fi:
    f.append(f"VCN_c{fin:02d}_Full_Z.pkl")

# cols = ['w', 'm', 'c', 'y', 'g', 'r', 'b', pg.mkPen()]
syms = ['s', 'o', 'x', 's', 'o', 'x', 's', 'o', 'x']

if pg:
    pg.mkQApp()
    win = pg.GraphicsWindow()
    win.setGeometry(100, 100, 500, 1200) 
    p1 = win.addPlot()
    win.nextRow()
    p2 = win.addPlot()
    win.nextRow()
    p3 = win.addPlot()
    p1.setLogMode(x=True, y=False)
    p2.setLogMode(x=True, y=False)
    p1.setLabels(title="Zin", left='Zin (Mohm)', bottom='Frequency(Hz)')
    p2.setLabels(title="Phase", left='Phase (pi radians)', bottom='Frequency(Hz)')
    p3.setLabels(title="Impedance Locus", left='-Zi (MOhm)', bottom='Zr (MOhm)')
    # legend = pg.LegendItem()
    # legend.setParentItem(p1)
    legend = p2.addLegend((80, 50), offset=(-1, 1))
    # print(dir(legend))
    # legend.setProperty({textsize:'8pt'})
    print('f: ', f)
    for i, filename in enumerate(f):
        with open(filename, 'rb') as fh:
            d = pickle.load(fh)
        col = pg.intColor(i, hues=len(f))
        pz = p1.plot(
            d['f'], d['zin'], pen=col, symbol=syms[i], symbolSize=3,
            name=filename)
        pp = p2.plot(
            d['f'], d['phase'], pen=col, symbol=syms[i], symbolSize=3,
            name=filename)
        zr = d['zin']*np.sqrt(1./(1+(np.tan(d['phase'])**2.)))
        zi = np.tan(d['phase'])*zr
    
        p3.plot(zr, -zi, pen=col, symbol=syms[i], symbolSize=3)
    # legendLabelStyle = {'color': '#FFF', 'size': '12pt', 'bold': True, 'italic': False}
    # for item in legend.items:
    #     for single_item in item:
    #         if isinstance(single_item, pg.graphicsItems.LabelItem.LabelItem):
    #             single_item.setText(single_item.text, **legendLabelStyle)

    pg.QtGui.QApplication.instance().exec_()


else:
    # matplotlib etc.
    # style = STY.styler('JNeurosci', "single")
    P = PH.regular_grid(3, 2, order='rowsfirst', 
        panel_labels=['A', 'B', 'C', 'D', 'E', 'F'],
        labelposition=(-0.05, 1.05))
    
    def plot_col(col, f, P):
        ax = P.axarr
        for i, a in enumerate(ax.ravel()):
            if i < 4:
                a.set_xscale('log')
        for i, filename in enumerate(f):
            print('col: ', col, i, 'file: ', filename)
            label = filename[:].replace('_', '\_')
            with open(filename, 'rb') as fh:
                d = pickle.load(fh)
            print('dkeys: ', d.keys()
            )# col = pg.intColor(i, hues=len(f))
            pz = ax[0, col].plot(
                d['f'], d['zin'], marker=syms[i], markersize=3,
                label=label[:8])
            ax[0, col].set_ylim(0, 100.)
            ax[0, col].set_ylabel("R (M$\Omega$)")
            ax[0, col].set_xlabel("Frequency (Hz)")
            pp = ax[1, col].plot(
                d['f'], d['phase'], marker=syms[i], markersize=3,
                #label=filename
                )
            ax[1, col].set_ylim(-1.5, 0.25)
            ax[1, col].set_ylabel("$\phi$ (radians)")
            ax[1, col].set_xlabel("Frequency (Hz)")
            
            zr = d['zin']*np.sqrt(1./(1+(np.tan(d['phase'])**2.)))
            zi = np.tan(d['phase'])*zr
    
            ax[2, col].plot(zr, -zi, marker=syms[i], markersize=3,)
            ax[2, col].set_ylim(-10., 60.)
            ax[2, col].set_xlim(5, 100)
            ax[2, col].set_ylabel("-Im(Z) (M$\Omega$)")
            ax[2, col].set_xlabel("Re(Z) (M$\Omega$)")
            
        if col == 0:
            ax[0, col].legend(fontsize='small', frameon=False,
                fancybox=False)
    # style.style_apply()
    for i, cond in enumerate(['Full', 'NoDend']):
        f = []
        for fin in fi:
            f.append(f"VCN_c{fin:02d}_{cond:s}_Z.pkl")
        print(len(f))
        plot_col(col=i, f=f, P=P)
    
    mpl.show()