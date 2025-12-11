"""
Read and prepare paths for the data and other input files.
All 
"""

from pathlib import Path
import toml
import ephys.tools.get_configuration as get_configuration


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
    print("*** get data paths")
    
    if Path("wheres_my_data.toml").is_file():
        with open("wheres_my_data.toml", "r") as fh:
            config = toml.load(fh)
        print("    Read config from toml file")
    else:   # try to read a configuration file from the parent
        dset, models = get_configuration.get_configuration("./config/models.cfg",
                                                    check_completeness=False)
        config = models[dset[0]]
        print("    Read config from models.cfg file")
    # config["basepath"] = config["baseDataDirectory"]
    # config["baseDataDirectory"] = Path(config["disk"], config["baseDataDirectory"])
    # config["cellDataDirectory"] = Path(config["disk"], config["cellDataDirectory"])
    # config["revcorrDataDirectory"] = Path(config["disk"], config["revcorrDataDirectory"])
    # config["baseMorphologyDirectory"] = Path(config["disk"], config["baseMorphologyDirectory"])
    
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