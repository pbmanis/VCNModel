from neuron import h
from cnmodel.cells import Bushy, SGC
import pyqtgraph as pg

# create a spiral ganglion cell and force it to spike four times
sgc = SGC.create(model='dummy', species='mouse')
sgc.set_spiketrain([51.0, 54.5, 60.0, 68.1])

# create a single bushy cell
bushy = Bushy.create(species='mouse', modelType='II')

# create a synapse from the spiral ganglion cell to the bushy cell
sgc.connect(bushy)

# initialize the cell
bushy.cell_initialize()

# perform a run
h.t = 0.
h.tstop = 100.
t = h.Vector()
v = h.Vector()
v.record(bushy.soma(0.5)._ref_v, sec=bushy.soma)
t.record(h._ref_t)
h.run()

# plot the data
pg.plot(t, v)
