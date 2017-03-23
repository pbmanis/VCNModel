"""
NOTE: FIle superseeded by cell_config.py
"""
    # synconfig consists of a list of tuples.
    # Each element in the list corresponds to one terminal and all of it's active zones
    # each tuple consists of N sites (calculated from area * average synapses/um2)
    # delay, and SR
    # the 4th and 5th entries in the tuple are the length of axon from the edge of the block (?)
    # the 6th entry is E or I (guess)
    # and the diameter. The brances are not included. distances are in microns
import numpy as np
import json
from collections import OrderedDict

synperum2 = 0.65 # average density of synapses, synapses per micron sequared


# Measurements:
# distances are in microns
# size is measured as radii (NOT diameter)
             # sites  del SR  axlen   axrad, branchlen, branchradius ,type(e or i)

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

VCN_c09= [  [ (53.46), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (78.81), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],  #partially on hillock
            [ (168.63), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (302.55), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (143.97), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (53.10), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (89.05), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (18.40), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            ]

VCN_c17= [  [ (149.80), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (173.82), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (82.66), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (261.75), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (42.97), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (44.77), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (108.71), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (212.03), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            ]

# VCN_c18 = [ [(216.66), 0., 2, np.nan, np.nan, 49.2, 1.222, 'e' ],
#             [(122.16), 0., 2, 82.7, 1.417, np.nan, np.nan, 'e'],
#             [(46.865), 0., 2, 67.3 , 1.309, 62.1, 0.717,  'e' ],
#             [(84.045), 0., 2, 22.4 , 1.416, 90.3, 0.924, 'e'],
#             [(80.27),  0., 2, np.nan, np.nan, 120.3,  0.687, 'e'],

VCN_c18= [  [ (249.94), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (155.80), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (61.30), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (98.37), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (5.41), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (5.02), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (19.72), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (63.27), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (5.38), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (5.74), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (34.03), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (14.74), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (13.45), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (1.08), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (3.59), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (32.49), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (0.87), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (115.99), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (64.30), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (6.71), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (17.08), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (15.66), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (19.77), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (26.14), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (7.11), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (16.71), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (6.55), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (12.71), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (6.82), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (6.96), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            ]

VCN_c19= [  [ (149.21), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (45.03), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (188.97), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (89.92), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (24.65), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (26.54), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (43.70), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (23.09), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (68.79), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (30.52), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (51.17), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (97.28), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (60.29), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (18.77), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (31.65), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            ]

VCN_c20= [  [ (81.21), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (65.35), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (13.94), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (90.40), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (54.67), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (52.14), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (61.36), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (44.32), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (150.80), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            ]

VCN_c21= [  [ (66.46), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (49.46), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (68.72), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (111.92), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (60.10), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (12.88), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (80.79), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (108.54), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (75.98), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (46.32), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            ]

VCN_c22= [  [ (132.22), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (134.56), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (121.92), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (53.04), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (134.79), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (141.62), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            [ (101.31), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],
            ]
        
VCN_Inputs = {'VCN_c07': VCN_c07, 'VCN_c08': VCN_c08,
              'VCN_c09': VCN_c09, 'VCN_c17': VCN_c17,
              'VCN_c18': VCN_c18, 'VCN_c19': VCN_c19,
              'VCN_c20': VCN_c20,
              'VCN_c21': VCN_c21, 'VCN_c22': VCN_c22
    }
        
def makeDict(cell, velocity=3.0):
    print (cell, VCN_Inputs.keys())
    assert cell in VCN_Inputs
    indata = VCN_Inputs[cell]
    r = [None]*len(indata)
    for j, i in enumerate(indata):
        r[j] = OrderedDict([('input', j+1), ('asa', i[0]),
            ('nSyn', int(np.around(i[0]*synperum2))),
            ('delay', i[1]), ('SR', i[2]),
            ("delay2", np.nansum([i[3],i[5]])*0.001/velocity),
            ('axonLen', i[3]), ('axonR', i[4]),
            ('branchLen', i[5]), ('branchR', i[6]), ('type', i[7])])
    
    return r


def printCellInputs(r):
    print(json.dumps(r, indent=4)) 
    

if __name__ == '__main__':
    # Check the formatting and display the results  
    r = makeDict('VCN_c18')
    printCellInputs(r)
    