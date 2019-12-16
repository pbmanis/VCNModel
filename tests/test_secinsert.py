import cnmodel.cells as cells
from channel_decorate import ChannelDecorate
import cell_config
cellID = 'VCN_c09'
b = cells.Bushy.create(morphology='VCN_Cells/VCN_c09/Morphology/VCN_c09.hoc', decorator=ChannelDecorate, modelType='mGBC')

#b.list_sections()

synapseConfig, celltype = cell_config.makeDict(cellID)

preCell = []
synapse = []

for i, syn in enumerate(synapseConfig):
    preCell.append(cells.DummySGC(cf=4000., sr=3))
    synapse.append(preCell[-1].connect(b, pre_opts={'nzones':syn['nSyn'], 'delay':syn['delay2']},
        post_opts={'AMPAScale': 1.0, 'postlocation': syn['postlocation']}))
for i, s in enumerate(synapse):
    s.terminal.relsite.Dep_Flag = 0  # turn off depression computation

b.print_connections()


