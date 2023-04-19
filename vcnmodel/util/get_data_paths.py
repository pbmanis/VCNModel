"""
Read and prepare paths for the data and other input files.
All 
"""

from pathlib import Path
import toml


def get_data_paths():
    """get the full path for a the data files
    Reads the "wheres_my_data.toml" file, prepends the disk,
    and returns the configuration.

    (A call to this functions should replace all instances of "with open('wheres_my_data))...")

    Parameters
    ----------
        None
    Returns
    -------
    Path object
        Full path to the output file
    """    
    with open("../wheres_my_data.toml", "r") as fh:
        config = toml.load(fh)
    config["basepath"] = config["baseDataDirectory"]
    config["baseDataDirectory"] = Path(config["disk"], config["baseDataDirectory"])
    config["cellDataDirectory"] = Path(config["disk"], config["cellDataDirectory"])
    config["revcorrDataDirectory"] = Path(config["disk"], config["revcorrDataDirectory"])
    config["baseMorphologyDirectory"] = Path(config["disk"], config["baseMorphologyDirectory"])
    
    return config

def update_disk(filename, datapaths):
    pfs = str(filename)
    match = pfs.find(str(datapaths["basepath"]))
    pfs = Path(datapaths["disk"], pfs[match:])
    return pfs

if __name__ == "__main__":
    # print the paths
    config = get_data_paths()
    print(config)