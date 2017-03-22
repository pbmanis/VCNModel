# cell_config generates the configuration of synaptic inputs and devines a cell
#
# the synaptic input structure consists of a list of tuples.
# Each element in the list corresponds to one terminal and all of it's active zones
# each tuple consists of N sites (calculated from area * average synapses/um2)
# delay, and SR
# the 4th and 5th entries in the tuple are the length of axon from the edge of the block (?)
# the 6th entry is E or I (guess)
# the 7th entry is the location to insert the synapse, defaulting to 'soma'
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
    
# VCN_c18 = [ [(216.66), 0., 2, np.nan, np.nan, 49.2, 1.222, 'e', 'soma' ],
#             [(122.16), 0., 2, 82.7, 1.417, np.nan, np.nan, 'e', 'soma'],
#             [(46.865), 0., 2, 67.3 , 1.309, 62.1, 0.717,  'e' ],
#             [(84.045), 0., 2, 22.4 , 1.416, 90.3, 0.924, 'e', 'soma'],
#             [(80.27),  0., 2, np.nan, np.nan, 120.3,  0.687, 'e', 'soma'],
#             ]

# based on new data 12/17/2016, from Spirou (table in Final7_Somatic_Input_ASAs.ods)
VCN_c18 = [ 
    
            [(249.94), 0., 2, np.nan, np.nan, 49.2, 1.222, 'e', 'soma' ],
            [(155.8), 0., 2, 82.7, 1.417, np.nan, np.nan, 'e', 'soma'], 
            [(115.99),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            [(98.37), 0., 2, 22.4 , 1.416, 90.3, 0.924, 'e', 'soma'],
            [(64.3),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            [(63.27),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            [(61.3), 0., 2, 67.3 , 1.309, 62.1, 0.717,  'e' ],
            [(34.03),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            [(32.49),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            [(26.14),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            # [(19.77),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],  # George says not to consider these yet - are of unknown origin
            # [(19.72),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            # [(17.08),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            # [(16.71),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            # [(15.66),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            # [(14.74),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            # [(13.45),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            # [(12.71),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            # [(7.11),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            # [(6.86),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            # [(6.82),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            # [(6.71),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            # [(6.55),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            # [(5.74),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            # [(5.41),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            # [(5.38),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            # [(5.02),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            # [(3.59),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            # [(1.08),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            # [(0.87),  0., 2, np.nan, np.nan, 1,  1, 'e', 'soma'],
            ]

# George says cell 7 is not a GBC
VCN_c07 = [ [(15.08), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [(12.95), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [(7.93), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [(9.13), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [(5.15), 0., 2., np.nan, np.nan,np.nan, np.nan, 'e', 'soma'],
            [(8.77), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
        ]
        
VCN_c08 = [
            [(174.6), 0., 2., 67.7, np.nan, 14.1, np.nan, 'e', 'soma'],
            [(127.44), 0., 2., 82.2, np.nan, 160.3, np.nan, 'e', 'soma'],
            [(89.99), 0., 2., 36.9, np.nan, 94.0, np.nan, 'e', 'soma'],
            [(66.26), 0., 2., 11.1, np.nan, 16.1, np.nan, 'e', 'soma'],
            [(61.83),0., 2., np.nan, np.nan, 131.3, np.nan, 'e', 'soma'],
            [(61.45), 0., 2., 36.2, np.nan, 22.5, np.nan, 'e', 'soma'],
            [(43.54), 0., 2., np.nan, np.nan, 20.5, np.nan, 'e', 'soma'],
            [(6.43), 0., 2., np.nan, np.nan, 25.5, np.nan, 'e', 'soma'],
            [(4.56), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [(1.37), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            ]


VCN_c09 = [
            [(302.55), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(168.63), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(143.97), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(89.05), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(78.81), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(53.46), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(53.10), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            ]

VCN_c99 = [  # VCN_c09 without any dendrites
            [(302.55), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(168.63), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(143.97), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(89.05), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(78.81), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(53.46), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(53.10), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            ]

VCN_c10 = [
            [(190.61), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(189.65), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(184.59), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(170.10), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(147.10), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(145.07), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(114.02), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(90.13), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(89.52), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(79.91), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(64.21), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(42.41), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            ]

VCN_c11 = [
            [(114.10), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(61.70), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(55.88), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            ]

VCN_c13 = [
            [(226.54), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(209.74), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(180.73), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(159.34), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(147.89), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(134.78), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(128.53), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(31.99), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            ]

VCN_c14 = [
            [(242.50), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(128.20), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(104.90), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(89.10), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(61.90), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(6.90), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            ]

VCN_c15 = [
            [(236.43), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(221.16), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(193.71), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(121.09), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(48.80), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(42.83), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(26.92), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(24.86), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(24.50), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(20.77), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
    ]

VCN_c16 = [
            [(385.57), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(287.97), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(247.41), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(192.02), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(152.74), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(92.84), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(55.07), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(32.79), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(22.96), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            ]

VCN_c17 = [
            [(159.93), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(135.32), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(127.19), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(67.20), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            [(16.62), 0., 2, np.nan, np.nan, np.nan, np.nan, 'e', 'soma' ],
            ]

VCN_c19= [  [ (188.97), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (149.21), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (97.28), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (89.92), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (68.79), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (60.29), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (51.17), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (45.03), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (43.70), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (31.65), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (30.52), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (26.54), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (24.65), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (23.09), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (18.77), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            ]

VCN_c20= [ 
            [ (150.80), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (90.40), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (81.21), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (65.35), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (61.36), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (54.67), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (52.14), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (44.32), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (13.94), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            ]

VCN_c21= [
            [ (111.92), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (108.54), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (80.79), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (75.98), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (68.72), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (66.46), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (60.10), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (49.46), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (46.32), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (12.88), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            ]

VCN_c22= [ 
            [ (141.62), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (134.79), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (134.56), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (132.22), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (121.92), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (101.31), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            [ (53.04), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'soma'],
            ]

        
VCN_Inputs = {'VCN_c07': ['notbushy', VCN_c07],
        'VCN_c08': ['bushy', VCN_c08], 'VCN_c09': ['bushy', VCN_c09],
        'VCN_c99': ['bushy', VCN_c99],
        'VCN_c10': ['bushy', VCN_c10], 'VCN_c11': ['bushy', VCN_c11],
        'VCN_c13': ['bushy', VCN_c13],
        'VCN_c14': ['bushy', VCN_c14], 'VCN_c15': ['bushy', VCN_c15],
        'VCN_c16': ['bushy', VCN_c16], 'VCN_c17': ['bushy', VCN_c17],
        'VCN_c18': ['bushy', VCN_c18], 
        'VCN_c19': ['bushy', VCN_c19], 'VCN_c20': ['bushy', VCN_c20],
        'VCN_c21': ['bushy', VCN_c21], 'VCN_c22': ['bushy', VCN_c22],
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