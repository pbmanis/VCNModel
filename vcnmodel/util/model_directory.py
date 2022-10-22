"""
Tool for making sure that the VCN cells are properly populated
with hoc files, and to list the simulations in the cell directories

Also see make_shared_datasets.py in the util subdirectory.

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2021, 2022 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 
"""
import datetime
import operator
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Union

import model_params
import toml


@dataclass
class SIM:
    name: Union[Path, str, None] = None
    date: float = 0.0
    datestr: str = ""
    dendscaling: bool = False
    somascaling: bool = False


@dataclass
class Morph:
    name: Union[Path, str, None] = None
    date: float = 0.0
    datestr: str = ""


class ModelDirectory(object):
    def __init__(self):

        where_is_data = Path("wheres_my_data.toml")
        if where_is_data.is_file():
            self.datapaths = toml.load("wheres_my_data.toml")
        else:
            self.datapaths = {
                "baseDirectory": Path("../VCN-SBEM-Data", "VCN_Cells")
            }  # "here"
        self.baseDirectory = self.datapaths["baseDirectory"]
        self.morphDirectory = "Morphology"
        self.initDirectory = "Initialization"
        self.simDirectory = "Simulations"

        # get a list of all simulations

    def get_sims(self, protocol="IV"):
        dates = []
        # print(self.baseDirectory)
        dirs = Path(self.baseDirectory).glob("VCN_c*")
        vcndirs = sorted(list(dirs))
        for d in vcndirs:
            sim_d = Path(d, self.simDirectory)
            # print('\n')
            if sim_d.is_dir:
                sim_IV_d = Path(sim_d, protocol)
                iv_sims = sorted(list(sim_IV_d.glob("*.p")))
                # print(iv_sims)
                for f in iv_sims:
                    mtime = f.stat().st_mtime
                    datestr = datetime.datetime.fromtimestamp(mtime).strftime(
                        "%Y.%m.%d-%H:%M:%S"
                    )
                    # print(mtime)
                    u = SIM(name=f, date=mtime, datestr=datestr)
                    dates.append(u)
        dates = sorted(dates, key=operator.attrgetter("date"))
        for f in range(len(dates)):
            print(self.datestr(dates[f].date), dates[f].name.name)

    def datestr(self, mtime):
        return datetime.datetime.fromtimestamp(mtime).strftime("%Y.%m.%d-%H:%M:%S")

    def backup_hoc(self):
        dates = []
        dirs = Path(self.baseDirectory).glob("VCN_c*")
        vcndirs = sorted(list(dirs))
        for d in vcndirs:
            # print(d)
            morph_d = Path(d, self.morphDirectory)
            # print('\n')
            if morph_d.is_dir:
                print(morph_d)
                morph_files = sorted(list(morph_d.glob("*.hoc")))
                # morph_files = sorted(list(morph_d.glob('*_Full.hoc')))
                morph_files.extend(sorted(list(morph_d.glob("*.hocx"))))
                morph_files.extend(sorted(list(morph_d.glob("*.swc"))))
                # print(iv_sims)
                for f in morph_files:
                    mtime = f.stat().st_mtime
                    datestr = datetime.datetime.fromtimestamp(mtime).strftime(
                        "%Y.%m.%d-%H:%M:%S"
                    )
                    # print(mtime)
                    u = Morph(name=f, date=mtime, datestr=datestr)
                    dates.append(u)
        dates = sorted(dates, key=operator.attrgetter("date"))
        for f in range(len(dates)):
            print(self.datestr(dates[f].date), dates[f].name.name)
        return morph_files

    def backup_morphology(self):
        """
        Reads through the morphology directories under VCN_Cells, and creates a
        date-named subdirectory that has a copy of the current .hoc, .hocx and
        .swc files Files retain a creation and modification date that are set to
        the modification dates in the parent directory (creation dates cannot be
        easily kept, unless we use rsync)

        """
        dates = []
        dirs = Path(self.baseDirectory).glob("VCN_c*")
        vcndirs = sorted(list(dirs))
        for d in vcndirs:
            # print(d)
            morph_d = Path(d, self.morphDirectory)
            print("morphd: ", morph_d)  # print('\n')
            if morph_d.is_dir:
                print("directory: ", morph_d)
                morph_files = sorted(list(morph_d.glob("*_Full.hoc")))
                morph_files.extend(sorted(list(morph_d.glob("*.hoc"))))
                morph_files.extend(sorted(list(morph_d.glob("*.hocx"))))
                morph_files.extend(sorted(list(morph_d.glob("*.swc"))))
            today = datetime.datetime.today()
            today_dirname = today.strftime("%Y.%m.%d")
            bkdir = Path(d, self.morphDirectory, today_dirname)
            print("bkdir: ", bkdir)
            bkdir.mkdir(parents=True, exist_ok=True)
            for f in morph_files:
                print("from: ", f)
                print("to: ", Path(d, bkdir, f.name))

            for f in morph_files:
                shutil.copy2(str(f), str(Path(d, bkdir, f.name)))


def main():
    MD = ModelDirectory()
    # MD.get_sims(protocol='IV')
    MD.backup_morphology()


if __name__ == "__main__":
    main()
