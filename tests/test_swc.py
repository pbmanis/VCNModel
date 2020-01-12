from __future__ import print_function
import neuronvis.swc_to_hoc as swc

print(dir(swc))
SWC = swc.SWC
soma = SWC('MorphologyFiles/cellbody.swc', types={1:'soma', 2:'axon', 3:'dendrite'})
soma.set_type(1)
axon = SWC('MorphologyFiles/axonnonscaled.swc')
axon.set_type(2)
dend = SWC('MorphologyFiles/dendnonscaled.swc')
dend.set_type(3)

soma.connect(57, axon)
soma.connect(39, dend)

soma.scale(0.11, 0.11, 0.06, 0.11)

soma.write_hoc('MorphologyFiles/test.hoc')
soma.topology()