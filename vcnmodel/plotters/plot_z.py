import numpy as np
import pickle
import pyqtgraph as pg

f = ['VCN_c09_Full_Z.pkl', 'VCN_c09_NoUninnervated_Z.pkl', 'VCN_c09_NoDend_Z.pkl']
fi = [2, 5, 6, 9, 10, 11, 13, 17, 30]
f = []
for fin in fi:
    f.append(f"VCN_c{fin:02d}_Full_Z.pkl")

# cols = ['w', 'm', 'c', 'y', 'g', 'r', 'b', pg.mkPen()]
syms = ['s', 'o', 'x', 's', 'o', 'x', 's', 'o', 'x']
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