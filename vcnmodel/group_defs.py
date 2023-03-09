"""data definitions that are widely used

"""
from typing import List, Union

# all the best "grade A" cells
gradeACells = [2, 5, 6, 9, 10, 11, 13, 17, 18, 30]
MixedMode = [9, 11, 17, 18]  # according to Figure 2, supplemental figure 2
Coincidence = [2, 5, 6, 10, 13, 30]  # according to figure 2, supplemental figure 3

# The fit groups are categorized differently than the mm/coincidence groups.
# These are based on efficacy, specifically for the weaker inputs that define group 2
Fit_Group1 = [9, 11, 17]
Fit_Group2 = [2, 5, 6, 10, 13, 18, 30]

# convience function
def grAList() -> List:
    return gradeACells

def get_cells_in_fit_group(group:int):
    if group == 1:
        return Fit_Group1
    elif group == 2:
        return Fit_Group2
    else:
        raise ValueError("Group must be either 1 or 2")

def get_fit_group(cellnumber):
    if cellnumber in Fit_Group1:
        return 1
    elif cellnumber in Fit_Group2:
        return 2

# Groups as defined in Figure 5E : mm or coincidence
def get_BC_Group(cellnumber):
    if cellnumber in MixedMode:
        return 2
    else:
        return 1

def get_cells_in_group(group:int):
    if group == 1:
        return Coincidence
    elif group == 2:
        return MixedMode
    else:
        raise ValueError("Group must be either 1 or 2")

def get_cells_by_group(group:str):
    if group == "Coincidence":
        return Coincidence
    elif group in ["MixedMode", "First-in"]:
        return MixedMode
    else:
        raise ValueError("Group must be either Mixed Mode/First-in or Coiincidence")

# remap the grouping
def remap_Group(cell):
    if cell in MixedMode:
        return "Mixed-mode"
    if cell in Coincidence:
        return "Coincidence"
    else:
        raise ValueError("Encoutnered cell not in the group definitions")

group_symbols = ['o', 'D'] 

def get_group_symbol(identity: Union[str, int]) -> str:
    if isinstance(identity, int):
        g = get_BC_Group(identity)
        if g not in [1,2]:
            raise ValueError("Group must be in 1 or 2")
        return group_symbols[g-1]

    elif isinstance(identity, str):
        groupsyms = {'Coincidence': group_symbols[0], 'First-in': group_symbols[1], "MixedMode": group_symbols[1]}
        return groupsyms[identity]        

#Standard map of colors to use for cells

sns_colors = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
    (1.0, 0.4980392156862745, 0.054901960784313725),
    (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
    (0.8392156862745098, 0.15294117647058825, 0.1568627450980392),
    (0.5803921568627451, 0.403921568627451, 0.7411764705882353),
    (0.5490196078431373, 0.33725490196078434, 0.29411764705882354),
    (0.8901960784313725, 0.4666666666666667, 0.7607843137254902),
    (0.4980392156862745, 0.4980392156862745, 0.4980392156862745),
    (0.7372549019607844, 0.7411764705882353, 0.13333333333333333),
    (0.09019607843137255, 0.7450980392156863, 0.8117647058823529),
]

