"""
Figure data: a dict of which simulations are in which figures

The dictionaries all should have the same basic structure: top level keys are
the cell number/ID. Second level keys relate to the experimental conditions. The
data in the second level key is either the name of the directory that holds the
simulation data file, or the name of a pickled file that holds some intermediate
results

This file is used both for plotting, and for extracting the original data for
plots into a new directory structure for general release. Plotting is controlled
by plot_sims.py and figures.py. Data extraction is controlled by
util/make_shared_dataset.py.

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2017-2022 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 
"""

import VS_datasets_15dB
import VS_datasets_15dB_BC09_NoUninnervated
import VS_datasets_30dB

BC_name = "GBC"  # or "BC" in the original version

figure_IV = {
    17: {
        "normal": "runIV-all-2021-07-27.12-18-24",
        "pasdend": "runIV-all-2021-07-27.12-23-27",
        "actdend": "runIV-all-2021-07-27.12-28-28",
        "Z_normal": "VCN_c17_Full_MeshInflate_normal_Z.pkl",
        "Z_passive": "VCN_c17_Full_MeshInflate_pasdend_Z.pkl",
        "Z_active": "VCN_c17_Full_MeshInflate_actdend_Z.pkl",
    }
}

"""
Updated 8/1/2021 with runs from 7/27/2021, 12:14-12:29pm
"""
figure_AllIVs = {
    # using standard axon ******************
    2: {
        "normal": "runIV-all-2021-08-01.13-47-29",
        "pasdend": "runIV-all-2021-08-01.13-52-43",
        "actdend": "runIV-all-2021-08-01.13-57-57",
    },
    # using standard axon ******************
    5: {
        "normal": "runIV-all-2021-08-01.13-48-01",
        "pasdend": "runIV-all-2021-08-01.13-53-16",
        "actdend": "runIV-all-2021-08-01.13-58-29",
    },
    6: {
        "normal": "runIV-all-2021-08-01.13-48-34",
        "pasdend": "runIV-all-2021-08-01.13-53-48",
        "actdend": "runIV-all-2021-08-01.13-59-02",
    },
    9: {
        "normal": "runIV-all-2021-08-01.13-49-02",
        "pasdend": "runIV-all-2021-08-01.13-54-15",
        "actdend": "runIV-all-2021-08-01.13-59-30",
    },
    10: {
        "normal": "runIV-all-2021-08-01.13-49-36",
        "pasdend": "runIV-all-2021-08-01.13-54-48",
        "actdend": "runIV-all-2021-08-01.14-00-02",
    },
    11: {
        "normal": "runIV-all-2021-08-01.13-50-09",
        "pasdend": "runIV-all-2021-08-01.13-55-22",
        "actdend": "runIV-all-2021-08-01.14-00-33",
    },
    13: {
        "normal": "runIV-all-2021-08-01.13-50-35",
        "pasdend": "runIV-all-2021-08-01.13-55-49",
        "actdend": "runIV-all-2021-08-01.14-00-59",
    },
    17: {
        "normal": "runIV-all-2021-08-01.13-51-07",
        "pasdend": "runIV-all-2021-08-01.13-56-20",
        "actdend": "runIV-all-2021-08-01.14-01-31",
    },
    18: {
        "normal": "runIV-all-2021-08-01.13-51-42",
        "pasdend": "runIV-all-2021-08-01.13-56-53",
        "actdend": "runIV-all-2021-08-01.14-02-02",
    },
    30: {
        "normal": "runIV-all-2021-08-01.13-52-12",
        "pasdend": "runIV-all-2021-08-01.13-57-26",
        "actdend": "runIV-all-2021-08-01.14-02-31",
    },
}

"""
The efficacy data is taken from runs using the latest measure
of synapse density, 0.7686 syn/um2  11/15/2020
Two tables here; one at 30 dB (singles) and one at 0dB (singles)
"""
figure_efficacy_supplement_30dB = {
    2: {
        "Full": "runANSingles-all-2021-08-02.17-58-19",
        "NoDend": "runANSingles-all-2021-08-02.23-13-48",
    },
    5: {
        "Full": "runANSingles-all-2021-08-02.18-07-35",
        "NoDend": "runANSingles-all-2021-08-02.23-17-51",
    },
    6: {
        "Full": "runANSingles-all-2021-08-02.18-21-35",
        "NoDend": "runANSingles-all-2021-08-02.23-22-53",
    },
    9: {
        "Full": "runANSingles-all-2021-08-02.18-30-57",
        "NoDend": "runANSingles-all-2021-08-02.23-27-30",
    },
    10: {
        "Full": "runANSingles-all-2021-08-02.18-45-37",
        "NoDend": "runANSingles-all-2021-08-02.23-34-04",
    },
    11: {
        "Full": "runANSingles-all-2021-08-02.19-06-04",
        "NoDend": "runANSingles-all-2021-08-02.23-41-21",
    },
    13: {
        "Full": "runANSingles-all-2021-08-02.19-16-40",
        "NoDend": "runANSingles-all-2021-08-02.23-46-44",
    },
    17: {
        "Full": "runANSingles-all-2021-08-02.19-27-45",
        "NoDend": "runANSingles-all-2021-08-02.23-51-15",
    },
    18: {
        "Full": "runANSingles-all-2021-08-02.19-42-24",
        "NoDend": "runANSingles-all-2021-08-02.23-57-14",
    },
    30: {
        "Full": "runANSingles-all-2021-08-02.19-56-49",
        "NoDend": "runANSingles-all-2021-08-03.00-03-56",
    },
}

figure_efficacy_supplement_Spont = {
    2: {
        "Full": "runANSingles-all-2021-07-31.11-13-17",
    },
    5: {
        "Full": "runANSingles-all-2021-07-31.11-22-28",
    },
    6: {
        "Full": "runANSingles-all-2021-07-31.11-36-25",
    },
    9: {
        "Full": "runANSingles-all-2021-07-31.11-45-33",
    },
    10: {
        "Full": "runANSingles-all-2021-07-31.12-00-01",
    },
    11: {
        "Full": "runANSingles-all-2021-07-31.12-20-18",
    },
    13: {
        "Full": "runANSingles-all-2021-07-31.12-30-51",
    },
    17: {
        "Full": "runANSingles-all-2021-07-31.12-41-53",
    },
    18: {
        "Full": "runANSingles-all-2021-07-31.12-56-36",
    },
    30: {
        "Full": "runANSingles-all-2021-07-31.13-11-02",
    },
}




figure_revcorr_example = {
    17: {
        "Spont": "runANPSTH-all-2020-11-24.09-12-49",
        "40dB": "runANPSTH-all-2020-08-20.14-54-39",
    }
}

# figure_revcorr = {
#     2: {
#         "Spont": "runANPSTH-all-2020-11-24.06-59-18",
#         "40dB": "runANPSTH-all-2020-11-24.11-13-17",
#     },
#     5: {
#         "Spont": "runANPSTH-all-2020-11-24.07-15-41",
#         "40dB": "runANPSTH-all-2020-11-24.11-19-05",
#     },
#     6: {
#         "Spont": "runANPSTH-all-2020-11-24.07-38-28",
#         "40dB": "runANPSTH-all-2020-11-24.11-27-20",
#     },
#     9: {
#         "Spont": "runANPSTH-all-2020-11-24.07-55-04",
#         "40dB": "runANPSTH-all-2020-11-24.11-33-12",
#     },
#     10: {
#         "Spont": "runANPSTH-all-2020-11-24.08-15-28",
#         "40dB": "runANPSTH-all-2020-11-24.11-40-19",
#     },
#     11: {
#         "Spont": "runANPSTH-all-2020-11-24.08-39-48",
#         "40dB": "runANPSTH-all-2020-11-24.11-47-36",
#     },
#     13: {
#         "Spont": "runANPSTH-all-2020-11-24.08-55-43",
#         "40dB": "runANPSTH-all-2020-11-24.11-52-44",
#     },
#     17: {
#         "Spont": "runANPSTH-all-2020-11-24.09-12-49",
#         "40dB": "runANPSTH-all-2020-11-24.11-58-11",
#         "50dB": "runANPSTH-all-2020-08-20.14-54-39",
#     },
#     18: {
#         "Spont": "runANPSTH-all-2020-11-24.09-37-01",
#         "40dB": "runANPSTH-all-2020-11-24.12-05-53",
#     },
#     30: {
#         "Spont": "runANPSTH-all-2020-11-24.09-51-06",
#         "40dB": "runANPSTH-all-2020-11-24.12-10-28",
#     },
# }

figure_revcorr = {
    2: {
        "Spont": "runANPSTH-all-2021-08-05.21-36-21",
        "30dB": "runANPSTH-all-2021-08-03.12-48-09",
    },
    5: {
        "Spont": "runANPSTH-all-2021-08-05.21-41-33",
        "30dB": "runANPSTH-all-2021-08-03.12-53-12",
    },
    6: {
        "Spont": "runANPSTH-all-2021-08-05.21-48-33",
        "30dB": "runANPSTH-all-2021-08-03.13-00-08",
    },
    9: {
        "Spont": "runANPSTH-all-2021-08-05.21-53-41",
        "30dB": "runANPSTH-all-2021-08-03.13-05-00",
    },
    10: {
        "Spont": "runANPSTH-all-2021-08-05.21-59-59",
        "30dB": "runANPSTH-all-2021-08-03.13-10-54",
    },
    11: {
        "Spont": "runANPSTH-all-2021-08-05.22-07-14",
        "30dB": "runANPSTH-all-2021-08-03.13-17-50",
    },
    13: {
        "Spont": "runANPSTH-all-2021-08-05.22-12-15",
        "30dB": "runANPSTH-all-2021-08-03.13-22-54",
    },
    17: {
        "Spont": "runANPSTH-all-2021-08-05.22-17-36",
        "30dB": "runANPSTH-all-2021-08-03.13-28-16",
        # "50dB": "runANPSTH-all-2020-08-20.14-54-39",
    },
    18: {
        "Spont": "runANPSTH-all-2021-08-05.22-25-16",
        "30dB": "runANPSTH-all-2021-08-03.13-35-32",
    },
    30: {
        "Spont": "runANPSTH-all-2021-08-05.22-31-41",
        "30dB": "runANPSTH-all-2021-08-03.13-41-59",
    },
}


figure_VClamp = {
    17: {
        "passive": "runVC-all-2020-07-29.10-36-59",
        "normal": "runVC-all-2020-07-29.10-30-30",
        "active": "runVC-all-2020-07-29.12-17-13",
    },
}

# figure_psth = {
#     2: {"40dB": "runANPSTH-all-2020-11-24.15-39-05",},
#     5: {"40dB": "runANPSTH-all-2020-11-24.15-46-45",},
#     6: {"40dB": "runANPSTH-all-2020-11-24.15-57-32",},
#     9: {"40dB": "runANPSTH-all-2020-11-24.16-05-10",},
#     10: {"40dB": "runANPSTH-all-2020-11-24.16-14-34",},
#     11: {"40dB": "runANPSTH-all-2020-11-24.16-25-29",},
#     13: {"40dB": "runANPSTH-all-2020-11-24.16-32-59",},
#     17: {"40dB": "runANPSTH-all-2020-11-24.16-40-56",},
#     18: {"40dB": "runANPSTH-all-2020-11-24.16-51-50",},
#     30: {"40dB": "runANPSTH-all-2020-11-24.16-58-24",},
# }
figure_psth = {
    2: {
        "30dB": "runANPSTH-all-2021-08-07.10-59-50",
    },
    5: {
        "30dB": "runANPSTH-all-2021-08-07.11-05-56",
    },
    6: {
        "30dB": "runANPSTH-all-2021-08-07.11-14-07",
    },
    9: {
        "30dB": "runANPSTH-all-2021-08-07.11-20-03",
    },
    10: {
        "30dB": "runANPSTH-all-2021-08-07.11-27-22",
    },
    11: {
        "30dB": "runANPSTH-all-2021-08-07.11-35-47",
    },
    13: {
        "30dB": "runANPSTH-all-2021-08-07.11-41-38",
    },
    17: {
        "30dB": "runANPSTH-all-2021-08-07.11-47-48",
    },
    18: {
        "30dB": "runANPSTH-all-2021-08-07.11-56-23",
    },
    30: {
        "30dB": "runANPSTH-all-2021-08-07.12-03-55",
    },
}

figure_SAC = {
    2: {
        "all": "runANPSTH-all-2022-01-02.19-56-54",
        "largestonly": "runANPSTH-largestonly-2022-01-02.20-02-29",
        "removelargest": "runANPSTH-removelargest-2022-01-02.20-08-14",
    },
    5: {
        "all": "runANPSTH-all-2022-01-02.20-14-23",
        "largestonly": "runANPSTH-largestonly-2022-01-02.20-21-12",
        "removelargest": "runANPSTH-removelargest-2022-01-02.20-27-52",
    },
    6: {
        "all": "runANPSTH-all-2022-01-02.20-34-56",
        "largestonly": "runANPSTH-largestonly-2022-01-02.20-39-59",
        "removelargest": "runANPSTH-removelargest-2022-01-02.20-45-00",
    },
    9: {
        "all": "runANPSTH-all-2022-01-02.20-50-30",
        "largestonly": "runANPSTH-largestonly-2022-01-02.20-56-36",
        "removelargest": "runANPSTH-removelargest-2022-01-02.21-02-45",
    },
    10: {
        "all": "runANPSTH-all-2022-01-02.21-09-25",
        "largestonly": "runANPSTH-largestonly-2022-01-02.21-16-25",
        "removelargest": "runANPSTH-removelargest-2022-01-02.21-23-22",
    },
    11: {
        "all": "runANPSTH-all-2022-01-02.21-30-35",
        "largestonly": "runANPSTH-largestonly-2022-01-02.21-35-34",
        "removelargest": "runANPSTH-removelargest-2022-01-02.21-40-28",
    },
    13: {
        "all": "runANPSTH-all-2022-01-02.21-45-49",
        "largestonly": "runANPSTH-largestonly-2022-01-02.21-51-00",
        "removelargest": "runANPSTH-removelargest-2022-01-02.21-56-14",
    },
    17: {
        "all": "runANPSTH-all-2022-01-02.22-02-02",
        "largestonly": "runANPSTH-largestonly-2022-01-02.22-09-30",
        "removelargest": "runANPSTH-removelargest-2022-01-02.22-16-35",
    },
    18: {
        "all": "runANPSTH-all-2022-01-02.22-24-05",
        "largestonly": "runANPSTH-largestonly-2022-01-02.22-30-05",
        "removelargest": "runANPSTH-removelargest-2022-01-02.22-36-07",
    },
    30: {
        "all": "runANPSTH-all-2022-01-02.22-42-35",
        "largestonly": "runANPSTH-largestonly-2022-01-02.22-48-47",
        "removelargest": "runANPSTH-removelargest-2022-01-02.22-54-58",
    },
}

figure_SAM_SAC = {
    11: {
        "SAM": "runANPSTH-all-2021-12-05.21-59-17",
        "SAC": "runANPSTH-all-2022-01-02.21-30-35",
    },
    17: {
        "SAM": "runANPSTH-all-2021-12-06.00-03-21",
        "SAC": "runANPSTH-all-2022-01-02.22-02-02",
    }
}

figure_No_Dend = {
    9: {
        "normal_noUn": "runIV-all-2022-05-23.13-52-18",
        "normal_Full": "runIV-all-2021-08-01.13-49-02", # same as control one in allIVs
    }
}

figure_cell9_nouninnervated2 = {
    9: {
        "NoUninnervated2": "runANSingles-all-2022-05-23.14-34-45",
        "Full": "runANSingles-all-2022-04-15.12-09-24",
    }
}

# we get the datasets here from the other files
figure_VS_SAM_30 = VS_datasets_30dB.samdata
figure_VS_SAM_15 = VS_datasets_15dB.samdata
figure_VS_SAM_BC09 = VS_datasets_15dB_BC09_NoUninnervated.samdata

all_figures = {
    "AN_PSTH": figure_psth,
    "VC_ex": figure_VClamp,
    "AN_revcorr": figure_revcorr,
    "AN_ex_revcorr": figure_revcorr_example,
    "AN_efficacy": figure_efficacy_supplement_30dB,
    "IV_ex": figure_IV,
    "IV_all": figure_AllIVs,
    "AN_SAC": figure_SAC,
    "AN_SAM_SAC": figure_SAM_SAC,
    "AN_BC_09_Pruned": figure_cell9_nouninnervated2,
    "AN_VS_15": figure_VS_SAM_15,
    "AN_VS_30": figure_VS_SAM_30,
    "AN_VS_15_BC09": figure_VS_SAM_BC09,
}



if __name__ == "__main__":
    pass
