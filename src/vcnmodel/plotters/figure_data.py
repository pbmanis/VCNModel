"""
Figure data: a dict of which simulations are in which figures

"""

figure_IV = {
    "Cell": 17,
    "normal": "runIV-all-2021-07-20.11-04-10",
    "pasdend": "runIV-all-2021-07-20.11-08-06",
    "actdend": "runIV-all-2021-07-20.11-11-59",
    "Z_normal": "VCN_c17_Full_MeshInflate_normal_Z.pkl",
    "Z_passive": "VCN_c17_Full_MeshInflate_pasdend_Z.pkl",
    "Z_active": "VCN_c17_Full_MeshInflate_actdend_Z.pkl",
}

figure_AllIVs = {
    2: {  # using standard axon **********************
        "normal": "runIV-all-2021-07-20.11-01-30",
        "pasdend": "runIV-all-2021-07-20.11-05-23",
        "actdend": "runIV-all-2021-07-20.11-09-20",
    },
    5: {  # using standard axon **********************
        "normal": "runIV-all-2021-07-20.11-01-51",
        "pasdend": "runIV-all-2021-07-20.11-05-46",
        "actdend": "runIV-all-2021-07-20.11-09-41",
    },
    6: {
        "normal": "runIV-all-2021-07-20.11-02-18",
        "pasdend": "runIV-all-2021-07-20.11-06-12",
        "actdend": "runIV-all-2021-07-20.11-10-08",
    },
    9: {
        "normal": "runIV-all-2021-07-20.11-02-37",
        "pasdend": "runIV-all-2021-07-20.11-06-32",
        "actdend": "runIV-all-2021-07-20.11-10-27",
    },
    # 10: {
    #     "normal": "runIV-all-2020-07-30.18-25-42",
    #     "passive": "runIV-all-2020-07-30.18-38-31",
    #     "active": "runIV-all-2020-07-30.18-52-07",
    # },
    10: {
        "normal": "runIV-all-2021-07-20.11-03-01",
        "pasdend": "runIV-all-2021-07-20.11-06-55",
        "actdend": "runIV-all-2021-07-20.11-10-51",
    },
    11: {
        "normal": "runIV-all-2021-07-20.11-03-28",
        "pasdend": "runIV-all-2021-07-20.11-07-23",
        "actdend": "runIV-all-2021-07-20.11-11-17",
    },
    13: {
        "normal": "runIV-all-2021-07-20.11-03-47",
        "pasdend": "runIV-all-2021-07-20.11-07-42",
        "actdend": "runIV-all-2021-07-20.11-11-37",
    },
    17: {
        "normal": "runIV-all-2021-07-20.11-04-10",
        "pasdend": "runIV-all-2021-07-20.11-08-06",
        "actdend": "runIV-all-2021-07-20.11-11-59",
    },
    18: { # new reconstruction June 15 2021 - fixed swc->hoc
        "normal": "runIV-all-2021-07-20.11-04-36",
        "pasdend": "runIV-all-2021-07-20.11-08-33",
        "actdend": "runIV-all-2021-07-20.11-12-24",
    },
    30: {
        "normal": "runIV-all-2021-07-20.11-04-59",
        "pasdend": "runIV-all-2021-07-20.11-08-55",
        "actdend": "runIV-all-2021-07-20.11-12-47",
    },
}

"""
The efficacy data is taken from runs using the latest measure
of synapse density, 0.7686 syn/um2  11/15/2020
"""
figure_efficacy_supplement = {
    2: {
        "NoDend": "runANSingles-all-2020-11-16.17-04-23",
        #"Full": "runANSingles-all-2020-11-16.17-08-55",
        "Full": "runANSingles-all-2021-07-19.17-41-17",
    },
    5: {
        "NoDend": "runANSingles-all-2020-11-16.17-19-30",
        #"Full": "runANSingles-all-2020-11-16.17-25-11",
        "Full": "runANSingles-all-2021-07-19.17-51-57",
    },
    6: {
        "NoDend": "runANSingles-all-2020-11-16.17-40-50",
        # "Full": "runANSingles-all-2020-11-16.17-46-05",
         "Full": "runANSingles-all-2021-07-19.18-07-39",
    },
    9: {
        "NoDend": "runANSingles-all-2020-11-16.17-56-43",
        # "Full": "runANSingles-all-2020-11-16.18-04-06",
         "Full": "runANSingles-all-2021-07-19.18-18-08",
    },
    10: {
        "NoDend": "runANSingles-all-2020-11-16.18-20-31",
        # "Full": "runANSingles-all-2020-11-16.18-28-43",
         "Full": "runANSingles-all-2021-07-19.18-34-51",
    },
    11: {
        "NoDend": "runANSingles-all-2020-11-16.18-51-40",
        # "Full": "runANSingles-all-2020-11-16.18-57-43",
        "Full": "runANSingles-all-2021-07-19.18-57-54",
    },
    13: {
        "NoDend": "runANSingles-all-2020-11-16.19-09-30",
        # "Full": "runANSingles-all-2020-11-16.19-14-35",
        "Full": "runANSingles-all-2021-07-19.19-09-53",
    },
    17: {
        "NoDend": "runANSingles-all-2020-11-16.19-27-02",
        # "Full": "runANSingles-all-2020-11-16.19-33-47",
        "Full": "runANSingles-all-2021-07-19.19-22-32",
    },
    18: {
        "NoDend": "runANSingles-all-2020-11-16.19-50-22",
        # "Full": "runANSingles-all-2021-05-18.14-43-58",  # new run
        # "Full": "runANSingles-all-2020-11-16.19-57-52",
        "Full": "runANSingles-all-2021-07-19.19-39-15",
    },
    30: {
        "NoDend": "runANSingles-all-2020-11-16.20-09-36",
        # "Full": "runANSingles-all-2020-11-16.20-16-56",
        "Full": "runANSingles-all-2021-07-19.19-55-51",
    },
}


figure_revcorr_example = {
    17: {
        "Spont": "runANPSTH-all-2020-11-24.09-12-49",
        "40dB": "runANPSTH-all-2020-08-20.14-54-39",
    }
}

figure_revcorr = {
    2: {
        "Spont": "runANPSTH-all-2020-11-24.06-59-18",
        "40dB": "runANPSTH-all-2020-11-24.11-13-17",
    },
    5: {
        "Spont": "runANPSTH-all-2020-11-24.07-15-41",
        "40dB": "runANPSTH-all-2020-11-24.11-19-05",
    },
    6: {
        "Spont": "runANPSTH-all-2020-11-24.07-38-28",
        "40dB": "runANPSTH-all-2020-11-24.11-27-20",
    },
    9: {
        "Spont": "runANPSTH-all-2020-11-24.07-55-04",
        "40dB": "runANPSTH-all-2020-11-24.11-33-12",
    },
    10: {
        "Spont": "runANPSTH-all-2020-11-24.08-15-28",
        "40dB": "runANPSTH-all-2020-11-24.11-40-19",
    },
    11: {
        "Spont": "runANPSTH-all-2020-11-24.08-39-48",
        "40dB": "runANPSTH-all-2020-11-24.11-47-36",
    },
    13: {
        "Spont": "runANPSTH-all-2020-11-24.08-55-43",
        "40dB": "runANPSTH-all-2020-11-24.11-52-44",
    },
    17: {
        "Spont": "runANPSTH-all-2020-11-24.09-12-49",
        "40dB": "runANPSTH-all-2020-11-24.11-58-11",
        "50dB": "runANPSTH-all-2020-08-20.14-54-39",
    },
    18: {
        "Spont": "runANPSTH-all-2020-11-24.09-37-01",
        "40dB": "runANPSTH-all-2020-11-24.12-05-53",
    },
    30: {
        "Spont": "runANPSTH-all-2020-11-24.09-51-06",
        "40dB": "runANPSTH-all-2020-11-24.12-10-28",
    },
}


figure_VClamp = {
    17: [
        "runVC-all-2020-07-29.10-36-59",
        "runVC-all-2020-07-29.10-30-30",
        "runVC-all-2020-07-29.12-17-13",
    ],
}

figure_psth = {
    2: {"40dB": "runANPSTH-all-2020-11-24.15-39-05",},
    5: {"40dB": "runANPSTH-all-2020-11-24.15-46-45",},
    6: {"40dB": "runANPSTH-all-2020-11-24.15-57-32",},
    9: {"40dB": "runANPSTH-all-2020-11-24.16-05-10",},
    10: {"40dB": "runANPSTH-all-2020-11-24.16-14-34",},
    11: {"40dB": "runANPSTH-all-2020-11-24.16-25-29",},
    13: {"40dB": "runANPSTH-all-2020-11-24.16-32-59",},
    17: {"40dB": "runANPSTH-all-2020-11-24.16-40-56",},
    18: {"40dB": "runANPSTH-all-2020-11-24.16-51-50",},
    30: {"40dB": "runANPSTH-all-2020-11-24.16-58-24",},
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
