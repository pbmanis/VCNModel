"""
List all of the data in the directories of the SBEM cells This mostly is to test
the integrity of the data Note that a number of file 'types' are also handled,
based on the data structure. Some of these are quite old formats. Very old data
(well, sort of) had "params" structures embedded. These are skipped. The data
structures have subsequently moved to dataclasses and dicts, which can be read
by standard library routines.

"""

__author__ = "pbmanis"
import os
import pickle
import pprint
import sys
from pathlib import Path

import toml


def list_dirs():
    config = toml.load(open("wheres_my_data.toml", "r"))

    print("Starting at top directory: ", config["cellDataDirectory"])
    alldirs = sorted(list(Path(config["cellDataDirectory"]).glob("*")))
    for d in alldirs:
        if d.is_dir():
            rootDir = Path(d, "Simulations/")
            list_dirs_for_cell(rootDir)


def list_dirs_for_cell(celldir):

    filesanddirs = celldir.rglob("*")
    for f in list(filesanddirs):
        if f.is_dir():  # and f != celldir:
            print("Found directory: %s" % f)
            list_dirs_for_cell(f)

        if f.is_file():
            if f.name in [".DS_Store"]:
                continue
            print("found file: ", f)
            with open(f, "rb") as fh:

                try:
                    d = pickle.load(fh, encoding="latin1")
                except:
                    print("   Pylibrary.params cannot be read....")
                    break
                print("File: %s\n", f.name)
                try:
                    if "runInfo" in d:
                        pprint.pprint(d["runInfo"], indent=4)
                    else:
                        pprint.pprint(d["Params"], indent=4)
                except:
                    if isinstance(d, dict):
                        pprint.pprint(d.keys(), indent=4)
                    elif isinstance(d, list):
                        pprint.pprint(d, indent=4)
                    else:
                        try:
                            pprint.pprint(d.__dict__, indent=4)
                        # pprint.pprint(d.runInfo, indent=4)
                        except:
                            raise


if __name__ == "__main__":
    list_dirs()
