"""
Figure data: a dict of which simulations are in which figures

"""

figure_IV = {
    "Cell": 17,
    "normal": "runIV-all-2020-07-30.18-29-53",
    "passive": "runIV-all-2020-07-30.18-43-29",
    "active": "runIV-all-2020-07-30.18-56-13",
    "Z_normal": "VCN_c17_Full_normal_Z.pkl",
    "Z_passive": "VCN_c17_Full_pasdend_Z.pkl",
    "Z_active": "VCN_c17_Full_actdend_Z.pkl",
}

figure_AllIVs = {
    2: {
        "normal": "runIV-all-2020-07-30.18-20-16",
        "passive": "runIV-all-2020-07-30.18-33-07",
        "active": "runIV-all-2020-07-30.18-46-46",
    },
    5: {
        "normal": "runIV-all-2020-07-30.18-21-32",
        "passive": "runIV-all-2020-07-30.18-34-20",
        "active": "runIV-all-2020-07-30.18-47-59",
    },
    6: {
        "normal": "runIV-all-2020-07-30.18-23-12",
        "passive": "runIV-all-2020-07-30.18-36-04",
        "active": "runIV-all-2020-07-30.18-49-38",
    },
    9: {
        "normal": "runIV-all-2020-07-30.18-24-19",
        "passive": "runIV-all-2020-07-30.18-37-10",
        "active": "runIV-all-2020-07-30.18-50-44",
    },
    # 10: {
    #     "normal": "runIV-all-2020-07-30.18-25-42",
    #     "passive": "runIV-all-2020-07-30.18-38-31",
    #     "active": "runIV-all-2020-07-30.18-52-07",
    # },
    10: {
        "normal": "runIV-all-2021-07-14.14-13-45",
        "passive": "runIV-all-2021-07-14.14-14-12",
        "active": "runIV-all-2021-07-14.14-14-40",
    },
    11: {
        "normal": "runIV-all-2020-07-30.18-27-24",
        "passive": "runIV-all-2020-07-30.18-40-45",
        "active": "runIV-all-2020-07-30.18-53-50",
    },
    13: {
        "normal": "runIV-all-2020-07-30.18-28-30",
        "passive": "runIV-all-2020-07-30.18-42-00",
        "active": "runIV-all-2020-07-30.18-54-51",
    },
    17: {
        "normal": "runIV-all-2020-07-30.18-29-53",
        "passive": "runIV-all-2020-07-30.18-43-29",
        "active": "runIV-all-2020-07-30.18-56-13",
    },
    18: {
        "normal": "runIV-all-2021-06-16.13-11-52",
        "passive": "runIV-all-2021-06-16.13-12-14",
        "active": "runIV-all-2021-06-16.13-12-36",
    },
    # 18: {  # new reconstruction June 15 2021 - fixed swc->hoc
    #     'normal': 'runIV-all-2021-05-18.12-37-05' ,
    #     'passive': 'runIV-all-2021-05-18.12-37-29' ,
    #     'active': 'runIV-all-2021-05-18.12-37-53' ,
    # },
    # 18: {
    #     "normal": "runIV-all-2020-11-16.10-54-53",
    #     "active": "runIV-all-2020-11-16.10-55-23",
    #     "passive": "runIV-all-2020-11-16.10-55-08",
    # },
    30: {
        "normal": "runIV-all-2020-07-30.18-31-35",
        "passive": "runIV-all-2020-07-30.18-45-12",
        "active": "runIV-all-2020-07-30.18-57-54",
    },
}

"""
The efficacy data is taken from runs using the latest measure
of synapse density, 0.7686 syn/um2  11/15/2020
"""
figure_efficacy_supplement = {
    2: {
        "NoDend": "runANSingles-all-2020-11-16.17-04-23",
        "Full": "runANSingles-all-2020-11-16.17-08-55",
    },
    5: {
        "NoDend": "runANSingles-all-2020-11-16.17-19-30",
        "Full": "runANSingles-all-2020-11-16.17-25-11",
    },
    6: {
        "NoDend": "runANSingles-all-2020-11-16.17-40-50",
        "Full": "runANSingles-all-2020-11-16.17-46-05",
    },
    9: {
        "NoDend": "runANSingles-all-2020-11-16.17-56-43",
        "Full": "runANSingles-all-2020-11-16.18-04-06",
    },
    10: {
        "NoDend": "runANSingles-all-2020-11-16.18-20-31",
        "Full": "runANSingles-all-2020-11-16.18-28-43",
    },
    11: {
        "NoDend": "runANSingles-all-2020-11-16.18-51-40",
        "Full": "runANSingles-all-2020-11-16.18-57-43",
    },
    13: {
        "NoDend": "runANSingles-all-2020-11-16.19-09-30",
        "Full": "runANSingles-all-2020-11-16.19-14-35",
    },
    17: {
        "NoDend": "runANSingles-all-2020-11-16.19-27-02",
        "Full": "runANSingles-all-2020-11-16.19-33-47",
    },
    18: {
        "NoDend": "runANSingles-all-2020-11-16.19-50-22",
        "Full": "runANSingles-all-2021-05-18.14-43-58",  # new run
        # "Full": "runANSingles-all-2020-11-16.19-57-52",
    },
    30: {
        "NoDend": "runANSingles-all-2020-11-16.20-09-36",
        "Full": "runANSingles-all-2020-11-16.20-16-56",
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
