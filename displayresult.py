"""
displayresult shows the results of a model_run. Just enter the filename in the fn field
"""

import matplotlib.pyplot as plt
import numpy as np
import pickle as p
from cnmodel.util import sound

fn = 'VCN_Cells/VCN_c18/Simulations/AN/AN_Result_VCN_c18_delays_N020_040dB_4000.0_HS.p'
h = open(fn)
d = p.load(h)

# d has keys: ['inputSpikeTimes', 'somaVoltage', 'spikeTimes', 'time', 'dendriteVoltage', 'stimInfo', 'stimWaveform']

#print len(d['time'])
w = d['stimWaveform'][0][0].generate()
t = d['stimWaveform'][0][0].time
fig, ax = plt.subplots(3, 1)
ax[1,].plot(t, w)
print d['spikeTimes']
for i, st in enumerate(d['spikeTimes'].keys()):
    ax[0,].plot(d['spikeTimes'][st]/1000., i*np.ones(len( d['spikeTimes'][st])), 'o', color='b')
    
plt.show()