# cell_config generates the configuration of synaptic inputs and devines a cell
#
# the synaptic input structure consists of a list of tuples.
# Each element in the list corresponds to one terminal and all of it's active zones
# each tuple consists of N sites (calculated from area * average synapses/um2)
# delay, and SR
# the 4th and 5th entries in the tuple are the length of axon from the edge of the block (?)
# the 6th entry is E or I (guess)
# and the diameter. The brances are not included. distances are in microns
#
# The make_dict routine converts these into the lists needed by model_run to population the cell
# connections. Note that we also use this to define the cell type, which determines the ion channels
# that population the cell.
#
import numpy as np
import json
from collections import OrderedDict

synperum2 = 0.65 # average density of synapses, synapses per micron sequared

# Measurements:
# distances are in microns
# size is measured as radii (NOT diameter)
#    [(ASA), nsyn(calculated), delay, SRgroup, delay2, axonlength, branch length, e or i]
    
# VCN_c18 = [ [(216.66), 0., 2, np.nan, np.nan, 49.2, 1.222, 'e' ],
#             [(122.16), 0., 2, 82.7, 1.417, np.nan, np.nan, 'e'],
#             [(46.865), 0., 2, 67.3 , 1.309, 62.1, 0.717,  'e' ],
#             [(84.045), 0., 2, 22.4 , 1.416, 90.3, 0.924, 'e'],
#             [(80.27),  0., 2, np.nan, np.nan, 120.3,  0.687, 'e'],
#             ]

# based on new data 12/17/2016, from Spirou (table in Final7_Somatic_Input_ASAs.ods)
VCN_c18 = [ [(249.94), 0., 2, np.nan, np.nan, 49.2, 1.222, 'e' ],
            [(155.8), 0., 2, 82.7, 1.417, np.nan, np.nan, 'e'], 
            [(115.99),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            [(98.37), 0., 2, 22.4 , 1.416, 90.3, 0.924, 'e'],
            [(64.3),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            [(63.27),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            [(61.3), 0., 2, 67.3 , 1.309, 62.1, 0.717,  'e' ],
            [(34.03),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            [(32.49),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            [(26.14),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            # [(19.77),  0., 2, np.nan, np.nan, 1,  1, 'e'],  # George says not to consider these yet - are of unknown origin
            # [(19.72),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            # [(17.08),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            # [(16.71),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            # [(15.66),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            # [(14.74),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            # [(13.45),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            # [(12.71),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            # [(7.11),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            # [(6.86),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            # [(6.82),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            # [(6.71),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            # [(6.55),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            # [(5.74),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            # [(5.41),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            # [(5.38),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            # [(5.02),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            # [(3.59),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            # [(1.08),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            # [(0.87),  0., 2, np.nan, np.nan, 1,  1, 'e'],
            ]

# George says cell 7 is not a GBC
VCN_c07 = [ [(15.08), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [(12.95), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [(7.93), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [(9.13), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [(5.15), 0., 2., np.nan, np.nan,np.nan, np.nan, 'e'],
            [(8.77), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
        ]
        
VCN_c08 = [ [(66.26), 0., 2., 11.1, np.nan, 16.1, np.nan, 'e'],
            [(6.43), 0., 2., np.nan, np.nan, 25.5, np.nan, 'e'],
            [(174.6), 0., 2., 67.7, np.nan, 14.1, np.nan, 'e'],
            [(43.54), 0., 2., np.nan, np.nan, 20.5, np.nan, 'e'],
            [(89.99), 0., 2., 36.9, np.nan, 94.0, np.nan, 'e'],
            [(127.44), 0., 2., 82.2, np.nan, 160.3, np.nan, 'e'],
            [(61.45), 0., 2., 36.2, np.nan, 22.5, np.nan, 'e'],
            [(4.56), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [(1.37), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [(61.83),0., 2., np.nan, np.nan, 131.3, np.nan, 'e'],
            ]


VCN_c09 = [
    [(302.55), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(168.63), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(143.97), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(78.81), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(53.46), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(53.10), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    ]

VCN_c10 = [
    [(190.61), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(189.65), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(184.59), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(170.10), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(147.10), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(145.07), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(114.02), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(90.13), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(89.52), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(79.91), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(64.21), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(42.41), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    ]

VCN_c11 = [
    [(114.10), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(61.70), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(55.88), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    ]

VCN_c13 = [
    [(226.54), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(209.74), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(180.73), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(159.34), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(147.89), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(134.78), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(128.53), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(31.99), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    ]

VCN_c14 = [
    [(242.50), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(128.20), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(104.90), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(89.10), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(61.90), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(6.90), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    ]

VCN_c15 = [
    [(236.43), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(221.16), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(193.71), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(121.09), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(48.80), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(42.83), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(26.92), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(24.86), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(24.50), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(20.77), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    ]

VCN_c16 = [
    [(385.57), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(287.97), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(247.41), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(192.02), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(152.74), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(92.84), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(55.07), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(32.79), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(22.96), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    ]

VCN_c17 = [
    [(159.93), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(135.32), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(127.19), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(67.20), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    [(16.62), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e' ],
    ]



        
VCN_Inputs = {'VCN_c18': ['bushy', VCN_c18], 'VCN_c07': ['notbushy', VCN_c07],
        'VCN_c08': ['bushy', VCN_c08], 'VCN_c09': ['bushy', VCN_c09],
        'VCN_c10': ['bushy', VCN_c10], 'VCN_c11': ['bushy', VCN_c11],
        'VCN_c13': ['bushy', VCN_c13],
        'VCN_c14': ['bushy', VCN_c14], 'VCN_c15': ['bushy', VCN_c15],
        'VCN_c16': ['bushy', VCN_c16], 'VCN_c17': ['bushy', VCN_c17],
        }
        
def makeDict(cell, velocity=3.0):
    assert cell in VCN_Inputs
    indata = VCN_Inputs[cell][1]
    celltype = VCN_Inputs[cell][0]
    r = [None]*len(indata)
    for j, i in enumerate(indata):
        r[j] = OrderedDict([('input', j+1), ('asa', i[0]),
            ('nSyn', int(np.around(i[0]*synperum2))),
            ('delay', i[1]), ('SR', i[2]),
            ("delay2", np.nansum([i[3],i[5]])*0.001/velocity),
            ('axonLen', i[3]), ('axonR', i[4]),
            ('branchLen', i[5]), ('branchR', i[6]), ('type', i[7])])
    
    return r, celltype


def printCellInputs(r):
    print(json.dumps(r, indent=4))




if __name__ == '__main__':
    # Check the formatting and display the results  
# Check the formatting and display the results
    r, celltype = makeDict('VCN_c18')
    print 'Cell Type: ', celltype
    printCellInputs(r)