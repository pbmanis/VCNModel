"""
Test whether neuron works with small diameter dendrites (<0.1 um)

While this normally would not be a problem, these are encountered
in new high-resolution reconstructions from SBEM material.

"""


from neuron import h
from pprint import pprint


hocstruct = """
objref soma
soma = new SectionList()
objref axon
axon = new SectionList()
objref dendrite
dendrite = new SectionList()

create sections[5]

access sections[0]
soma.append()
sections[0] {
    pt3dadd(33.391190, 39.891441, -32.180580, 2.302326)  // connect dendrite 0.060241)
    pt3dadd(33.687325, 38.943192, -32.854713, 3.608622)
    pt3dadd(33.775993, 37.043259, -34.277031, 5.400022)
    pt3dadd(34.079544, 32.780701, -36.625359, 9.923010)
    pt3dadd(34.928024, 28.653046, -37.699585, 14.286328)
    pt3dadd(37.634308, 23.339478, -37.915794, 15.517240)
    pt3dadd(37.634308, 23.339478, -37.915794, 15.517240)
    pt3dadd(38.589157, 20.617432, -37.847702, 11.964212)
    pt3dadd(40.379379, 16.413605, -37.723038, 7.282564)
    pt3dadd(40.443516, 16.042000, -37.662342, 4.953794)
    pt3dadd(41.534229, 14.799618, -38.033463, 0.357617)
}

access sections[1]
dendrite.append()
connect sections[1](0), sections[0](0)
sections[1] {
  pt3dadd(-64.937929, 65.431734, 23.493839, 2.302326)
}

access sections[2]
dendrite.append()
connect sections[2](0), sections[1](1)
sections[2] {
  pt3dadd(-66.472813, 63.896850, 23.493839, 1.588886)
  pt3dadd(-67.496068, 62.617781, 23.633374, 1.191287)
  pt3dadd(-68.007696, 61.338711, 23.912444, 0.1)
  pt3dadd(-68.263510, 60.827083, 24.051979, 0.1)
  pt3dadd(-68.007696, 60.059641, 24.470583, 0.1)
  pt3dadd(-67.496068, 58.013129, 24.889188, 1.790698)
  pt3dadd(-67.240254, 56.222432, 24.889188, 1.790698)
}

access sections[3]
dendrite.append()
connect sections[3](0), sections[2](1)
sections[3] {
  pt3dadd(-66.984441, 51.873594, 24.889188, 0.01)
  pt3dadd(-67.984673, 49.556943, 24.889188, 0.01)
}

access sections[4]
dendrite.append()
connect sections[4](0), sections[3](1)
sections[4] {
  pt3dadd(-68.775138, 48.292199, 24.889188, 2.302326)
  pt3dadd(-69.030952, 47.268943, 25.447327, 2.032255)
  pt3dadd(-69.542580, 46.245688, 26.005467, 2.856291)
  pt3dadd(-69.542580, 45.478246, 26.005467, 2.385567)
}

"""

h(hocstruct)
h.topology()

for s in h.allsec():
    s.Ra = 150.
soma = h.sections[0]
distdend = h.sections[4]
soma.insert('hh')
ic = h.IClamp(soma(0.5))

ic.delay = 5.0
ic.dur = 2.0
ic.amp = 0.05
v = h.Vector().record(soma(0.5)._ref_v)
v2 = h.Vector().record(distdend(0.5)._ref_v)
t = h.Vector().record(h._ref_t)
h.finitialize(-65.0)
h.dt = 0.025

h.t = 0.0
tstop = 100.0
while h.t < tstop:
    h.fadvance()

import matplotlib.pyplot as mpl
mpl.plot(t, v)
mpl.plot(t, v2)
mpl.show()


