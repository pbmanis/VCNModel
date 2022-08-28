from __future__ import print_function
from pathlib import Path
import neuronvis.swc_to_hoc as swc

SWC = swc.SWC
soma = SWC('tests/test_data/cellbody.swc', types={1:'soma', 2:'axon', 3:'dendrite'})
soma.set_type(1)
axon = SWC('tests/test_data/axonnonscaled.swc')
axon.set_type(2)
dend = SWC('tests/test_data/dendnonscaled.swc')
dend.set_type(3)

soma.connect(57, axon)
soma.connect(39, dend)

soma.rescale(0.11, 0.11, 0.06, 0.11)

soma.write_hoc('tests/test_data/test.hoc')
soma.show_topology()
Path('tests/test_data/test.hoc').unlink()
Path('tests/test_data/test.segmap').unlink()