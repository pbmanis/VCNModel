"""
build output path for figures
"""

from pathlib import Path
from typing import Union
from vcnmodel.util.get_data_paths import get_data_paths

def set_figure_path(fignum:int, filedescriptor:str, suppnum:Union[int, None]=None, suffix:Union[str, None] = ".pdf"):
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

    config = get_data_paths()
        
    figpath = Path(config["disk"], config["baseDataDirectory"], config["figureDirectory"], f"Figure{fignum:d}")
    if suppnum is None:
        figpath = Path(figpath, f"Figure{fignum:d}_{filedescriptor:s}")
    else:
        figpath = Path(figpath, f"Figure{fignum:d}_supp", 
            f"Figure{fignum:d}_Supplemental{suppnum:d}_{filedescriptor:s}")
    if suffix is not None:
        figpath = figpath.with_suffix(suffix)
    return figpath
