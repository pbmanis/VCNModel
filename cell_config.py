"""cell_config generates the configuration of synaptic inputs and devines a cell

the synaptic input structure consists of a list of tuples.
Each element in the list corresponds to one terminal and all of it's active zones
each tuple consists of N sites (calculated from area * average synapses/um2)
delay, and SR
the 4th and 5th entries in the tuple are the length of axon from the edge of the block (?)
the 6th entry is E or I (guess)
the 7th entry is the location to insert the synapse, defaulting to {'soma': [0, 0.5, 1.0]}
and the diameter. The brances are not included. distances are in microns

The make_dict routine converts these into the lists needed by model_run to population the cell
connections. Note that we also use this to define the cell type, which determines the ion channels
that population the cell.

Measurements:
distances are in microns
size is measured as radii (NOT diameter)
   [(ASA), nsyn(calculated), delay, SRgroup, delay2, axonlength, branch length, syntype, postlocationdescriptor]

syntype is "AN" for auditory nerve, "AMPA", "AMPA+NMDA", "glycine", "GABA", or can be more specific as to actual mechanism

location descriptor is as follows:
{{'soma': [0, 0.5, 1.0]}: [0, 0.5, 1.0]}
The synapse might be split across sections as follows:
{'nameofsection': {nrnsection1#: [location, fractionofgmax], nrnsection2#: [location, fractionofgmax]}}
if fraction of gmax is -1, then it is computed as the residual of the remaining gmax.
(this allows things like a big ending that crosses anatomically defined boundaries)
"""

"""
Data sets are from George Spirou & Co. (Michael Morehead, Nathan Spencer)

Latest revsion 4/1/2018 from data sent 12/22/2017 by Nathan Spencer (ASA_12Cells.ods)
The following cells were updated:
Cell09 Cell11 Cell14 Cell16 Cell17
Cell18 Cell19 Cell20 Cell21 Cell22
Cell24 Cell29

"""
import numpy as np
import json
from collections import OrderedDict

synperum2 = 0.65 # average density of synapses, synapses per micron squared 
                 # Original value, from Spriou measurements in MNTB.

# VCN_c18 = [ [(216.66), 0., 2, np.nan, np.nan, 49.2, 1.222, 'AN', {'soma': [0, 0.5, 1.0]} ],
#             [(122.16), 0., 2, 82.7, 1.417, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
#             [(46.865), 0., 2, 67.3 , 1.309, 62.1, 0.717,  'AN' ],
#             [(84.045), 0., 2, 22.4 , 1.416, 90.3, 0.924, 'AN', {'soma': [0, 0.5, 1.0]}],
#             [(80.27),  0., 2, np.nan, np.nan, 120.3,  0.687, 'AN', {'soma': [0, 0.5, 1.0]}],
#             ]

# LM test cell from Luke Campagnola (mouse, CBA): Simple 3 input, 100 each (no exp. data on inputs)
LC_bushy = [ [(50./synperum2), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
             [(50./synperum2), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
             [(50./synperum2), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
        ]
# LM test cell from Luke Campagnola (mouse, CBA): Simple 3 input, 100 each (no exp. data on inputs)
LC_neuromantic_scaled = [ [(50./synperum2), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
             [(50./synperum2), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
             [(50./synperum2), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
        ]

# based on new data 12/17/2016, from Spirou (Table in Final7_Somatic_Input_ASAs.ods)

# George says cell 7 is not a GBC
VCN_c07 = [ [(15.08), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [(12.95), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [(7.93), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [(9.13), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [(5.15), 0., 2., np.nan, np.nan,np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [(8.77), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
        ]
        
VCN_c08 = [
            [(174.6), 0., 2., 67.7, np.nan, 14.1, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [(127.44), 0., 2., 82.2, np.nan, 160.3, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [(89.99), 0., 2., 36.9, np.nan, 94.0, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [(66.26), 0., 2., 11.1, np.nan, 16.1, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [(61.83),0., 2., np.nan, np.nan, 131.3, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [(61.45), 0., 2., 36.2, np.nan, 22.5, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [(43.54), 0., 2., np.nan, np.nan, 20.5, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [(6.43), 0., 2., np.nan, np.nan, 25.5, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [(4.56), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [(1.37), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            ]

VCN_c09 = [
            [(302.55), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(168.63), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(143.97), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(89.05), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(78.81), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(53.46), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1]} ],
            [(53.10), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            ]

VCN_c09_h = [  # input 6 split on hillock
            [(302.55), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(168.63), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(143.97), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(89.05), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(78.81), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(53.46), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 0.0], 'hillock': [210, 0.5, 1]}],
            [(53.10), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            ]

VCN_c09_nd = [  # VCN_c09 without any dendrites
            [(302.55), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(168.63), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(143.97), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(89.05), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(78.81), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(53.46), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(53.10), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            ]

VCN_c10 = [
            [(190.61), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(189.65), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(184.59), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(170.10), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(147.10), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(145.07), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(114.02), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(90.13), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(89.52), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(79.91), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(64.21), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(42.41), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            ]

VCN_c11 = [ 
            [(263.75), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e' 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(236.42), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e' 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(181.60), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e' 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(94.02), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e', 'AN', {'soma': [0, 0.5, 1.0]}],
            [(87.25), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e' 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(49.85), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e' 'AN', {'soma': [0, 0.5, 1.0]} ],
#            [(41.69), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e' 'AN', {'soma': [0, 0.5, 1.0]} ],  # Missing from "ASA_12 Cells" from Nathan Spencer, Dec 2017
            ]

VCN_c13 = [
            [(226.54), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(209.74), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(180.73), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(159.34), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(147.89), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(134.78), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(128.53), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(31.99), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            ]

VCN_c14 = [ [(309.83), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e' 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(121.64), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e' 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(150.36), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e' 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(19.74), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e' 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(70.32), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e' 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(80.39), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e' 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(83.04), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e' 'AN', {'soma': [0, 0.5, 1.0]} ],
#            [(30.87), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e' 'AN', {'soma': [0, 0.5, 1.0]} ],  # Missing from "ASA_12 Cells" from Nathan Spencer, Dec 2017
            ]

VCN_c15 = [
            [(236.43), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(221.16), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(193.71), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(121.09), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(48.80), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(42.83), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(26.92), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(24.86), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(24.50), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(20.77), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
        ]

# VCN_c16 = [  # Cell revised in "ASA_12 Cells" from Nathan Spencer, Dec 2017
#               Original based on new data 12/17/2016, from Spirou (table in Final7_Somatic_Input_ASAs.ods)
#             [(385.57), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],  
#             [(287.97), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
#             [(247.41), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
#             [(192.02), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
#             [(152.74), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
#             [(92.84), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
#             [(55.07), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
#             [(32.79), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
#             [(22.96), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
#             ]

VCN_c16 = [  # revised version
            [(218.437), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(198.641), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(101.207), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(100.633), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(65.3895), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            ]        

# VCN_c17 = [  # Cell revised in "ASA_12 Cells" from Nathan Spencer, Dec 2017
#               Original based on new data 12/17/2016, from Spirou (table in Final7_Somatic_Input_ASAs.ods)
#             [( 159.93), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
#             [( 135.32), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
#             [( 127.19), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
#             [( 67.20), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
#             [( 16.62), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
#             ]

VCN_c17 = [  # revised version
            [(261.749), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(212.031), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(173.817), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(149.802), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(108.713), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(82.663), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]} ],
            ]  
            
VCN_c18 = [ 
            [(249.94), 0., 2, np.nan, np.nan, 49.2, 1.222, 'AN', {'soma': [0, 0.5, 1.0]} ],
            [(155.8), 0., 2, 82.7, 1.417, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}], 
            [(115.99),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            [(98.37), 0., 2, 22.4 , 1.416, 90.3, 0.924, 'AN', {'soma': [0, 0.5, 1.0]}],
            [(64.3),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            [(63.27),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            [(61.3), 0., 2, 67.3 , 1.309, 62.1, 0.717,  'AN', {'soma': [0, 0.5, 1.0]}],
#               Original based on new data 12/17/2016, from Spirou (table in Final7_Somatic_Input_ASAs.ods)
            # [(34.03),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}], # Cell revised in "ASA_12 Cells" from Nathan Spencer, Dec 2017
            # [(32.49),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            # [(26.14),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            # [(19.77),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],  # George says not to consider these yet - are of unknown origin
            # [(19.72),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            # [(17.08),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            # [(16.71),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            # [(15.66),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            # [(14.74),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            # [(13.45),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            # [(12.71),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            # [(7.11),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            # [(6.86),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            # [(6.82),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            # [(6.71),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            # [(6.55),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            # [(5.74),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            # [(5.41),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            # [(5.38),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            # [(5.02),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            # [(3.59),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            # [(1.08),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            # [(0.87),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': [0, 0.5, 1.0]}],
            ]

# VCN_c19= [  Original version
#               Original based on new data 12/17/2016, from Spirou (table in Final7_Somatic_Input_ASAs.ods)
#             [ (   188.97), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
#             [ (   149.21), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
#             [ (   97.28), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
#             [ (   89.92), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
#             [ (   68.79), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
#             [ (   60.29), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
#             [ (   51.17), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
#             [ (   45.03), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
#             [ (   43.70), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
#             [ (   31.65), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
#             [ (   30.52), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
#             [ (   26.54), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
#             [ (   24.65), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
#             [ (   23.09), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
#             [ (   18.77), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
#             ]

VCN_c19= [ # Cell revised in "ASA_12 Cells" from Nathan Spencer, Dec 2017
            [ (188.972), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (149.212), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ ( 97.278), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ ( 89.920), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ ( 68.795), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ ( 60.292), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ ( 51.173), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],

            ]     

VCN_c20= [ 
            [ (150.80), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (90.40), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (81.21), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (65.35), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (61.36), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (54.67), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (52.14), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            # [ (44.32), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],  Cell revised in "ASA_12 Cells" from Nathan Spencer, Dec 2017
            # [ (13.94), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            ]

VCN_c21= [
            [ (111.92), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (108.54), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (80.79), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (75.98), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (68.72), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (66.46), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (60.10), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (49.46), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (46.32), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            # [ (12.88), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],  Cell revised in "ASA_12 Cells" from Nathan Spencer, Dec 2017
            ]

VCN_c22= [ 
            [ (141.62), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (134.79), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (134.56), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (132.22), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (121.92), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (101.31), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (53.04), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            ]

VCN_c24= [
            [ (281.332), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (136.014), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (131.227), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (123.791), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (118.848), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (106.619), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (72.043), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (57.603), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (47.749), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            ]         

VCN_c29= [
            [ (191.376), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (160.929), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (147.643), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (132.020), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (130.695), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            [ (101.560), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': [0, 0.5, 1.0]}],
            ]


VCN_Inputs = {'VCN_c07': ['notbushy', VCN_c07],
        'LC_bushy': ['bushy', LC_bushy], 
        'LC_neuromantic_scaled': ['bushy', LC_bushy], 
        'VCN_c08': ['bushy', VCN_c08],
        'VCN_c09': ['bushy', VCN_c09], 'VCN_c09h': ['bushy', VCN_c09_h], 'VCN_c09nd': ['bushy', VCN_c09_nd],  # variants
        'VCN_c10': ['bushy', VCN_c10],
        'VCN_c11': ['bushy', VCN_c11],
        'VCN_c13': ['bushy', VCN_c13],
        'VCN_c14': ['bushy', VCN_c14],
        'VCN_c15': ['bushy', VCN_c15],
        'VCN_c16': ['bushy', VCN_c16],
        'VCN_c17': ['bushy', VCN_c17],
        'VCN_c18': ['bushy', VCN_c18], 
        'VCN_c19': ['bushy', VCN_c19],
        'VCN_c20': ['bushy', VCN_c20],
        'VCN_c21': ['bushy', VCN_c21],
        'VCN_c22': ['bushy', VCN_c22],
        'VCN_c24': ['bushy', VCN_c24],
        'VCN_c29': ['bushy', VCN_c29],
        }

VCN_Inputs_Clean = {  # Spriou SBEM VCN Bushy cells only, no variants
        'VCN_c08': ['bushy', VCN_c08],
        'VCN_c09': ['bushy', VCN_c09],
        'VCN_c10': ['bushy', VCN_c10],
        'VCN_c11': ['bushy', VCN_c11],
        'VCN_c13': ['bushy', VCN_c13],
        'VCN_c14': ['bushy', VCN_c14],
        'VCN_c15': ['bushy', VCN_c15],
        'VCN_c16': ['bushy', VCN_c16],
        'VCN_c17': ['bushy', VCN_c17],
        'VCN_c18': ['bushy', VCN_c18], 
        'VCN_c19': ['bushy', VCN_c19],
        'VCN_c20': ['bushy', VCN_c20],
        'VCN_c21': ['bushy', VCN_c21],
        'VCN_c22': ['bushy', VCN_c22],
        'VCN_c24': ['bushy', VCN_c24],
        'VCN_c29': ['bushy', VCN_c29],
        }
# convert the table for a given cell to the dictionary used by model_run to decroate
# cell with synapses.

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
            ('branchLen', i[5]), ('branchR', i[6]), ('type', i[7]),
            ('postlocations', i[8]),
        ])
    return r, celltype


def printCellInputs(r):
    print(json.dumps(r, indent=4))


def summarize_inputs():
    import matplotlib.pyplot as mpl
    import pylibrary.PlotHelpers as PH
    P = PH.regular_grid(2, 3, order='columns', figsize=(10., 8), showgrid=False)
    ax = [P.axdict[x] for x in P.axdict.keys()]
    print(ax)
 #   PH.nice_plot(ax)
    allendings = []
    cellendings = {}
    for cell, v in iter(VCN_Inputs_Clean.items()):
        cellendings[cell] = []
        for s in v[1]:
            # print cell, s[0]
            allendings.append(s[0])
            cellendings[cell].append(s[0])
    ax[0].hist(allendings, bins=20)
    ax[0].set_title('All ending areas')
    ax[0].set_ylim((0, 25))
    ax[0].set_ylabel('N')
    ax[0].set_xlim((0,350))
    ax[0].set_xlabel('Area (um^2)')
    
    normd = []
    ratio1 = []
    ratio2 = []
    meansize = []
    maxsize = []
    convergence = []
    
    for cell in cellendings.keys():
            normd.extend(cellendings[cell]/np.max(cellendings[cell]))
            ratio1.append(cellendings[cell][1]/cellendings[cell][0])
            ratio2.append(np.mean(cellendings[cell])/cellendings[cell][0])
            meansize.append(np.mean(cellendings[cell]))
            maxsize.append(np.max(cellendings[cell]))
            convergence.append(len(cellendings[cell]))
    print('convergence: ', convergence)
    ax[1].set_title('Normaized by largest')
    ax[1].hist(normd, bins=20, range=(0,1.0), align='mid')
    ax[1].set_xlabel('Area (um^2)')
    ax[1].set_ylabel('N')

    ax[2].set_title('ratio largest to next largest')
    ax[2].hist(ratio1, bins=10, range=(0,1.0), align='mid')
    ax[2].set_xlabel('Area (um^2)')
    ax[2].set_title('ratio mean to largest')

    ax[3].set_title('Ratio of mean to largest')
    ax[3].hist(ratio2, bins=10, range=(0,1.0), align='mid')
    ax[3].set_xlabel('Area (um^2)')

    ax[4].set_xlim((0., 200.))
    ax[4].set_xlabel('Area (um^2)')
    ax[4].set_ylim((0., 15.))  
    ax[4].set_ylabel('Convergence')
    ax[4].set_title('Convergence vs. mean size')
    fit = np.polyfit(meansize, convergence, 1)
    fit_fn = np.poly1d(fit) 
    ax[4].scatter(meansize, convergence)
    ax[4].plot(meansize, fit_fn(meansize), '--k')


    ax[5].set_title('Convergence vs max size')
    fit = np.polyfit(maxsize, convergence, 1)
    fit_fn = np.poly1d(fit) 
    ax[5].scatter(maxsize, convergence)
    ax[5].plot(maxsize, fit_fn(maxsize), '--k')
    ax[5].set_xlim((0., 350.))
    ax[5].set_xlabel('Area (um^2)')
    ax[5].set_ylim((0., 15.))
    ax[5].set_xlabel('Convergence')

    

    mpl.show()
        


if __name__ == '__main__':
    # Check the formatting and display the results  
    summarize_inputs()
    exit(0)
    r, celltype = makeDict('VCN_c09')
    print 'Cell Type: ', celltype
    printCellInputs(r)