"""
Read and set paths
"""

from pathlib import Path
from typing import Union
import toml


def set_data_paths():
    """Set the full path for a figure, with a specific extension
    This is a tool to keep the figure files organized and consistently named.

    Parameters
    ----------
    fignum : int
        Figure number - for top directory "Figure#"
    filedescriptor : str
        A descriptor for the file, such as "Ephys_main_V2",
    suppnum : Union[int, None], optional
        if the figure will be a supplement, number for the supplement, by default None
        This puts the file in a subdirectory of the fignum (Figure#/Figure#_supp)
        with the name "Figure#_Supplemental#_filedescriptor"
    suffix: Union[str, None], optional
        Suffix to add, by default ".pdf"
        if the suffix is None, then no suffix is added.

    Returns
    -------
    Path object
        Full path to the output file
    """    
    with open("wheres_my_data.toml", "r") as fh:
        config = toml.load(fh)
    config.baseDataDirectory = Path(config["disk"], config["baseDataDirectory"])
    config.cellDataDirectory = Path(config.disk, config.cellDataDirectory)
    config.revcorrDataDirectory = Path(config.disk, config.revcorrDataDirectory)
    config.baseMorphologyDirectory = Path(config.disk, config.baseMorphologyDirectory)
    
    return config
