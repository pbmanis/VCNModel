"""make_data_dirs.py
Create a directory structure based on figure_data.py (in the plotters
subdirectory) and populate it with all of the simulation files and intermediate
analyses needed to regenerate the figures using DataTablesVCN.py This also
populates the Morphology directory for each cell with the .hoc and .swc
morphology files.

The directory created by this routine was compressed (zip) and uploaded
to figshare:
Private link: https://figshare.com/s/2f87300f58e1c1e9270c
(DOI: 10.6084/m9.figshare.19520041)

"""
import shutil
from pathlib import Path

import toml
from vcnmodel.plotters import (
    figure_data as FD,
)  # table of simulation runs used for plotting figures

config = toml.load(open("wheres_my_data.toml", "r"))
sourcepath = Path(config["baseDataDirectory"])

simpath = Path("/Volumes/Pegasus_002/BU_simulation_data")

intermediates = Path(simpath, "IntermediateAnalyses")
intermediates.mkdir(exist_ok=True, parents=True)

simulations = Path(simpath, "Simulations")
simulations.mkdir(exist_ok=True, parents=True)

Zin_dest_dir = Path(simulations, "Impedance_Calculations")
Zin_dest_dir.mkdir(exist_ok=True, parents=True)
Zin_source_dir = Path(sourcepath, "VCN_Cells", "Impedance_Calculations")

allcells = [2, 5, 6, 9, 10, 11, 13, 17, 18, 30]
cdirs = {}
sdirs = {}
for cn in allcells:
    cell_id = f"VCN_c{cn:02d}"
    celldir = Path(simulations, cell_id)
    cdirs[cn] = Path(simulations, celldir)
    cdirs[cn].mkdir(exist_ok=True)
    for ddir in ["Simulations/AN", "Simulations/IV", "Simulations/VC", "Morphology"]:
        Path(cdirs[cn], ddir).mkdir(exist_ok=True, parents=True)
    sdirs[cn] = Path(sourcepath, "VCN_Cells", cell_id)

# for p in list(simpath.glob("*/**")):
#     print(p)


def mk_cellid(cellN):
    return f"VCN_c{cellN:02d}"


def find_data(pkldir):
    pass


def copy_morphology(allcells):
    for cn in allcells:
        cell_id = f"VCN_c{cn:02d}"
        source_morph = list(Path(sdirs[cn], "Morphology").glob("*.hoc"))
        source_morph.extend(list(Path(sdirs[cn], "Morphology").glob("*.swc")))
        dest_morph = Path(cdirs[cn], "Morphology")
        for f in source_morph:
            dest_file = Path(dest_morph, f.name)
            print(f"File to copy: {str(dest_file):s}")
            if not dest_file.is_file():  # only copy if not already there
                print("copying ", f, " to ", dest_morph)
                shutil.copy2(f, dest_morph)
            else:
                print(f"Morphology file {cell_id:s}  {str(f):s}  already present")
                print("    dest file was: ", dest_file)


copy_morphology(allcells)

exit()

for fig in FD.all_figures:
    figd = FD.all_figures[fig]
    if fig in ["IV_ex", "IV_all", "VC_ex"]:
        for cell_id in list(figd.keys()):
            cell = mk_cellid(cell_id)
            sourcedir = Path(sdirs[cell_id], "Simulations", fig[:2])
            targetdir = Path(cdirs[cell_id], "Simulations", fig[:2])
            print("\n", "=" * 80, "\n", fig, " target dir: ", targetdir)
            for expt in figd[cell_id].keys():  # for each expt in the list
                print("   expt: ", expt)
                sourcefile = Path(sourcedir, figd[cell_id][expt])
                if expt.startswith("Z_"):  # capture impedance runs
                    sourcefile = Path(Zin_source_dir, figd[cell_id][expt])
                    if sourcefile.is_file and sourcefile.suffix == ".pkl":
                        print("   Zin pkl copying ", sourcefile, " to ", Zin_dest_dir)
                        shutil.copy2(sourcefile, Zin_dest_dir)

                elif sourcefile.is_dir():  # make and fill the subdirectory
                    targetsubdir = Path(targetdir, sourcefile.name)
                    targetsubdir.mkdir(exist_ok=True)
                    sourcefiles = sourcefile.glob("*.p")
                    for f in sourcefiles:
                        print("   data copying ", f, " to ", targetsubdir)
                        shutil.copy2(f, targetsubdir)
                elif sourcefile.is_file and sourcefile.suffix == ".pkl":
                    print("   suffix: ", sourcefile.suffix)
                    print("   pkl copying ", sourcefile, " to ", targetdir)
                    shutil.copy2(sourcefile, targetdir)
    elif fig in [
        "AN_PSTH",
        "AN_revcorr",
        "AN_ex_revcorr",
        "AN_efficacy",
        "AN_SAC",
        "AN_SAM_SAC",
    ]:
        for cell_id in list(figd.keys()):
            cell = mk_cellid(cell_id)
            sourcedir = Path(sdirs[cell_id], "Simulations", fig[:2])
            targetdir = Path(cdirs[cell_id], "Simulations", fig[:2])
            print("\n", "=" * 80, "\n", fig, " target dir: ", targetdir)
            for expt in figd[cell_id].keys():  # for each expt in the list
                print("   expt: ", expt)
                sourcefile = Path(sourcedir, figd[cell_id][expt])

                if sourcefile.is_dir():  # make and fill the subdirectory
                    targetsubdir = Path(targetdir, sourcefile.name)
                    targetsubdir.mkdir(exist_ok=True)
                    sourcefiles = sourcefile.glob("*.p")
                    for f in sourcefiles:
                        print("   data copying ", f, " to ", targetsubdir)
                        shutil.copy2(f, targetsubdir)
                elif sourcefile.is_file and sourcefile.suffix == ".pkl":
                    print("   suffix: ", sourcefile.suffix)
                    print("   pkl copying ", sourcefile, " to ", targetdir)
                    shutil.copy2(sourcefile, targetdir)
