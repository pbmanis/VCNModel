__author__ = 'pbmanis'

from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import numpy as np

class testcyl():
    def __init__(self):
        self.app = QtGui.QApplication([])
        self.w = gl.GLViewWidget()
        self.w.show()
        self.w.setWindowTitle('hoc Render')
        self.w.setCameraPosition(distance=1000.)

        self.g = gl.GLGridItem()
        self.g.scale(2,2,1)
        self.w.addItem(self.g)
        self.draw_model(mode = 'cylinder')
        self.show()

    def show(self):
        QtGui.QApplication.instance().exec_()


    def draw_model(self, mode = "cylinder"):
        """
        modified from:neuronvisio
        Draw the model.
        Params:
        controls - the main gui obj."""


        x,y,z,d = [], [], [], []
        voltage = []
        connections = []
        x = np.array([0., 1., 1])
        y = np.array([0., 1., 2.])
        z = np.array([0., 5., 10.])
        d = np.array([1, 0.5, 0.25])
        d = np.array(d) # Transforming for easy division
        connections.append((1,0))
        connections.append((2,1))
        connections.append((2,1))
        lines = np.vstack(connections)
        print 'lines: ', lines
        self.drawMeshes(x,y,z,d,lines, mode)
        self.drawMeshes(x,y,z,d,lines, 'line')


    def drawMeshes(self, x, y, z, d, lines, mode):
        print 'x: ', len(x)
        dmax = np.max(d)
        #wmax = 20.
        for c in range(len(lines)-1):
            i = lines[c, 0] # from
            j = lines[c, 1] # to
            if i < 3:
                print 'xyzd: %6.1f %6.1f %6.1f %6.2f' % (x[j], y[j], z[j], d[j])
            pts = np.vstack([[x[j], x[i]],[y[j], y[i]],[z[j],z[i]]]).transpose()
            print 'pts: ', pts
            if mode == "line":
                plt = gl.GLLinePlotItem(pos=pts, width =(d[i]+d[j])/(2.), color=pg.glColor((int(255.*d[i]/dmax), 128)), connected=True)
                self.w.addItem(plt)
            elif mode == 'sphere':
                md = gl.MeshData.sphere(rows=10, cols=20, radius=d[i]/2.0) # , length=d(i))
                colors = np.ones((md.faceCount(), 4), dtype=float)
                colors[::2,0] = 0
                colors[:,1] = np.linspace(0, 1, colors.shape[0])
                md.setFaceColors(colors)
                m5 = gl.GLMeshItem(meshdata=md, smooth=False, drawEdges=False,)
                m5.translate(x[j],y[j],z[j])
            elif mode == "cylinder":
                cyllen = np.sqrt((x[j]-x[i])**2.0 + (y[j]-y[i])**2.0 + (z[j]-z[i])**2.0)
                #print 'cyllen: ', cyllen
                md = gl.MeshData.cylinder(rows=10, cols=20, radius=[d[j]/2., d[i]/2.], length=cyllen)
                colors = np.ones((md.faceCount(), 4), dtype=float)
                colors[::2,0] = 0
                colors[:,1] = np.linspace(0, 1, colors.shape[0])
                #colors[:,1] = (int(255.*d[i]/dmax), 128)
                md.setFaceColors(colors)
                m5 = gl.GLMeshItem(meshdata=md, smooth=False, drawEdges=False)


                m5.translate(x[j],y[j],z[j]+cyllen/2.0) # move into position

                self.w.addItem(m5)
        #self.draw_mayavi(x, y, z, d, self.edges)

if __name__ == '__main__':
    T = testcyl()