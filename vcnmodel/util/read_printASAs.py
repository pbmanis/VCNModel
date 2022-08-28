
"""
Build synconfig (print to terminal) based on reading the ASA file
e.g., make this::

    VCN_c18 = [ [(216.66), 0., 2, np.nan, np.nan, 49.2, 1.222, 'e' ],
            [(122.16), 0., 2, 82.7, 1.417, np.nan, np.nan, 'e'], 
            [(46.865), 0., 2, 67.3 , 1.309, 62.1, 0.717,  'e' ],
            [(84.045), 0., 2, 22.4 , 1.416, 90.3, 0.924, 'e'],
            [(80.27),  0., 2, np.nan, np.nan, 120.3,  0.687, 'e'],
            ]

Copy the output and paste it into syn_config.py 

This code is obselete.

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2022 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 
"""
import toml
from pathlib import Path
import pandas as pd
import numpy as np
import vcnmodel.cell_config as cell_config


def main():
    config = toml.load("wheres_my_data.toml")
    cf = cell_config.CellConfig(
            verbose=True,
        )
    return

    datafile = cf.dendqualfile
    with open(datafile, "rb") as fh:
        SDSummary = pd.read_excel(
            fh, config["SomaAndDendriteData"], skiprows=3
        )

    with open(datafile, "rb") as fh:
        ASA = pd.read_excel(
            fh,
            config["asaData"],
            skiprows=config["asaHeaderSkip"],
           # engine="openpyxl",
        )

   # f = open(Path(config["baseMorphologyDirectory"], 'Final7_Somatic_Input_ASAs.xls'))
   #  f = open(Path(config["baseMorphologyDirectory"], 'Dendrite Quality and Surface Areas_comparisons_9_July_2021.xlsx'))
    # d = pd.read_excel(f, skiprows=3, skipfooter=12)
    # df = ASA
    # print(df.head())
    # df.drop(df.columns[8:],inplace=True,axis=1)

    # nh = df.iloc[0]
    # df = df[1:]
    # df = df.rename(columns = nh)

    # for dx in df.columns.values:
    #     s = ''
    #     if not isinstance(dx, int):
    #         continue
    #     u = df[dx].dropna()
    #     s += 'VCN_c%02d' % dx
    #     s += '= [ '
    #     for asa in u:
    #         s += "            [ (%.2f), 0., 2., np.nan, np.nan, np.nan, np.nan, 'e'],\n" % asa
    #     s += '            ]'
    #     print('')
    #     print(s)


if __name__ == "__main__":
    main()