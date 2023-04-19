"""
Rename import paths for pickle load after restructuring data 

Approach taken From stack overflow:
https://stackoferflow.com/questions/2121874/python-pickling-after-changing-a-modules-directory

This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01 DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Copyright 2020 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 
"""
import io
import pickle

import dill
import vcnmodel.analyzers.reverse_correlation
import vcnmodel.plotters.plot_sims



class RenameUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        renamed_module = module
        #  Previously had modules one deeper (src/vcnmodel); now just at top
        if module in [
            "src.vcnmodel",
            "src.vcnmodel.model_params",
            "src.vcnmodel.table_manager",
            "src.vcnmodel.plotters.plot_sims",
            "src.vcnmodel.plotters.figures",
            "src.vcnmodel.util",
            "src.vcnmodel.analyzers.reverse_correlation",
        ]:
            renamed_module = module  # used to be "src." + module
        if module.startswith("src."):
            renamed_module = module[4:]  # so strip src. from module name
        # print("renamed module: ", renamed_module)
        return super(RenameUnpickler, self).find_class(
            renamed_module,
            name,
        )


def pickle_load(file_obj):
    return RenameUnpickler(file_obj).load()


def pickle_loads(pickled_bytes):
    # file_obj = io.BytesIO(pickled_bytes)
    file_obj = dill.loads(pickled_bytes)
    return renamed_load(file_obj)
