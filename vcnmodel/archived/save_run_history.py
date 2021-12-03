from pathlib import Path
import time
import datetime
from dataclasses import dataclass, field
from typing import Union, List
import operator
import model_params
import toml

@dataclass
class SIM:
    name:Union[Path, str, None] = None
    date:float = 0.0
    datestr:str = ''
    dendscaling: bool=False
    somascaling: bool=False