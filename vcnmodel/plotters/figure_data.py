"""
Figure data: a dict of which simulations are in which figures

"""

figure_IV = {
    "Cell": 17,
    "normal": "runIV-all-2021-07-27.12-18-24",
    "pasdend": "runIV-all-2021-07-27.12-23-27",
    "actdend": "runIV-all-2021-07-27.12-28-28",
    "Z_normal": "VCN_c17_Full_MeshInflate_normal_Z.pkl",
    "Z_passive": "VCN_c17_Full_MeshInflate_pasdend_Z.pkl",
    "Z_active": "VCN_c17_Full_MeshInflate_actdend_Z.pkl",
}

"""
Updated 8/1/2021 with runs from 7/27/2021, 12:14-12:29pm
"""
figure_AllIVs = {
    # using standard axon ******
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
    # 10: {
    #     "normal": "runIV-all-2020-07-30.18-25-42",
    #     "passive": "runIV-all-2020-07-30.18-38-31",
    #     "active": "runIV-all-2020-07-30.18-52-07",
    # },
    10: {
        "normal": "runIV-all-2021-08-01.13-49-36",
        "pasdend": "runIV-all-2021-08-01.13-54-48",
        "actdend": "runIV-all-2021-08-01.14-00-02",
    },
    # unclear; cell 11 did not have all runs on 7/27 (missing pasdend)
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
"""
figure_efficacy_supplement = {
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
    17: [
        "runVC-all-2020-07-29.10-36-59",
        "runVC-all-2020-07-29.10-30-30",
        "runVC-all-2020-07-29.12-17-13",
    ],
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

all_figures = [
    figure_psth,
    figure_VClamp,
    figure_revcorr,
    figure_revcorr_example,
    figure_efficacy_supplement,
    figure_IV,
    figure_AllIVs,
]


if __name__ == '__main__':
    pass
