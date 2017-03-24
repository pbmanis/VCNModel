# cell_config generates the configuration of synaptic inputs and devines a cell
#
# the synaptic input structure consists of a list of tuples.
# Each element in the list corresponds to one terminal and all of it's active zones
# each tuple consists of N sites (calculated from area * average synapses/um2)
# delay, and SR
# the 4th and 5th entries in the tuple are the length of axon from the edge of the block (?)
# the 6th entry is E or I (guess)
# the 7th entry is the location to insert the synapse, defaulting to {'soma': {0: [0.5, 1.0]}}
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
#    [(ASA), nsyn(calculated), delay, SRgroup, delay2, axonlength, branch length, syntype, postlocationdescriptor]

# syntype is "AN" for auditory nerve, "AMPA", "AMPA+NMDA", "glycine", "GABA", or can be more specific as to actual mechanism

# location descriptor is as follows:
# {{'soma': {0: [0.5, 1.0]}}: {0: [0.5, 1.0]}}
# The synapse might be split across sections as follows:
# {'nameofsection': {nrnsection1#: [location, fractionofgmax], nrnsection2#: [location, fractionofgmax]}}
# if fraction of gmax is -1, then it is computed as the residual of the remaining gmax.
# (this allows things like a big ending that crosses anatomically defined boundaries)
    
# VCN_c18 = [ [(216.66), 0., 2, np.nan, np.nan, 49.2, 1.222, 'AN', {'soma': {0: [0.5, 1.0]}} ],
#             [(122.16), 0., 2, 82.7, 1.417, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
#             [(46.865), 0., 2, 67.3 , 1.309, 62.1, 0.717,  'AN' ],
#             [(84.045), 0., 2, 22.4 , 1.416, 90.3, 0.924, 'AN', {'soma': {0: [0.5, 1.0]}}],
#             [(80.27),  0., 2, np.nan, np.nan, 120.3,  0.687, 'AN', {'soma': {0: [0.5, 1.0]}}],
#             ]

# based on new data 12/17/2016, from Spirou (table in Final7_Somatic_Input_ASAs.ods)
VCN_c18 = [ 
    
            [(249.94), 0., 2, np.nan, np.nan, 49.2, 1.222, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(155.8), 0., 2, 82.7, 1.417, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}], 
            [(115.99),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [(98.37), 0., 2, 22.4 , 1.416, 90.3, 0.924, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [(64.3),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [(63.27),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [(61.3), 0., 2, 67.3 , 1.309, 62.1, 0.717,  'AN' ],
            [(34.03),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [(32.49),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [(26.14),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            # [(19.77),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],  # George says not to consider these yet - are of unknown origin
            # [(19.72),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            # [(17.08),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            # [(16.71),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            # [(15.66),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            # [(14.74),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            # [(13.45),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            # [(12.71),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            # [(7.11),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            # [(6.86),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            # [(6.82),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            # [(6.71),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            # [(6.55),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            # [(5.74),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            # [(5.41),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            # [(5.38),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            # [(5.02),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            # [(3.59),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            # [(1.08),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            # [(0.87),  0., 2, np.nan, np.nan, 1,  1, 'AN', {'soma': {0: [0.5, 1.0]}}],
            ]

# George says cell 7 is not a GBC
VCN_c07 = [ [(15.08), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [(12.95), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [(7.93), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [(9.13), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [(5.15), 0., 2., np.nan, np.nan,np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [(8.77), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
        ]
        
VCN_c08 = [
            [(174.6), 0., 2., 67.7, np.nan, 14.1, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [(127.44), 0., 2., 82.2, np.nan, 160.3, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [(89.99), 0., 2., 36.9, np.nan, 94.0, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [(66.26), 0., 2., 11.1, np.nan, 16.1, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [(61.83),0., 2., np.nan, np.nan, 131.3, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [(61.45), 0., 2., 36.2, np.nan, 22.5, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [(43.54), 0., 2., np.nan, np.nan, 20.5, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [(6.43), 0., 2., np.nan, np.nan, 25.5, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [(4.56), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [(1.37), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            ]


VCN_c09 = [
            [(302.55), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(168.63), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(143.97), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(89.05), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(78.81), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(53.46), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(53.10), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            ]

VCN_c99 = [  # VCN_c09 without any dendrites
            [(302.55), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(168.63), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(143.97), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(89.05), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(78.81), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(53.46), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(53.10), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            ]

VCN_c10 = [
            [(190.61), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(189.65), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(184.59), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(170.10), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(147.10), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(145.07), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(114.02), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(90.13), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(89.52), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(79.91), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(64.21), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(42.41), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            ]

VCN_c11 = [
            [(114.10), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(61.70), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(55.88), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            ]

VCN_c13 = [
            [(226.54), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(209.74), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(180.73), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(159.34), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(147.89), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(134.78), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(128.53), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(31.99), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            ]

VCN_c14 = [
            [(242.50), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(128.20), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(104.90), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(89.10), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(61.90), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(6.90), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            ]

VCN_c15 = [
            [(236.43), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(221.16), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(193.71), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(121.09), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(48.80), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(42.83), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(26.92), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(24.86), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(24.50), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(20.77), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
    ]

VCN_c16 = [
            [(385.57), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(287.97), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(247.41), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(192.02), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(152.74), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(92.84), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(55.07), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(32.79), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(22.96), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            ]

VCN_c17 = [
            [(159.93), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(135.32), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(127.19), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(67.20), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            [(16.62), 0., 2, np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}} ],
            ]

VCN_c19= [  [ (188.97), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (149.21), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (97.28), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (89.92), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (68.79), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (60.29), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (51.17), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (45.03), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (43.70), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (31.65), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (30.52), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (26.54), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (24.65), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (23.09), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (18.77), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            ]

VCN_c20= [ 
            [ (150.80), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (90.40), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (81.21), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (65.35), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (61.36), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (54.67), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (52.14), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (44.32), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (13.94), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            ]

VCN_c21= [
            [ (111.92), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (108.54), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (80.79), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (75.98), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (68.72), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (66.46), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (60.10), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (49.46), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (46.32), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (12.88), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            ]

VCN_c22= [ 
            [ (141.62), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (134.79), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (134.56), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (132.22), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (121.92), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (101.31), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
            [ (53.04), 0., 2., np.nan, np.nan, np.nan, np.nan, 'AN', {'soma': {0: [0.5, 1.0]}}],
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
            ('branchLen', i[5]), ('branchR', i[6]), ('type', i[7]),
            ('postlocation', i[8]),
        ])
    return r, celltype


def printCellInputs(r):
    print(json.dumps(r, indent=4))




if __name__ == '__main__':
    # Check the formatting and display the results  
# Check the formatting and display the results
    r, celltype = makeDict('VCN_c18')
    print 'Cell Type: ', celltype
    printCellInputs(r)