import io
import pickle

"""
Rename import paths for pickle load after restructuring data 

From stack overflow:
https://stackoverflow.com/questions/2121874/python-pickling-after-changing-a-modules-directory

"""
class RenameUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        renamed_module = module
        if module in ["vcnmodel","vcnmodel.model_params", 'vcnmodel.table_manager',
            "vcnmodel.plotters.plot_sims"]:
            renamed_module = "src." + module
        return super(RenameUnpickler, self).find_class(renamed_module, name)

def pickle_load(file_obj):
    return RenameUnpickler(file_obj).load()


def pickle_loads(pickled_bytes):
    file_obj = io.BytesIO(pickled_bytes)
    return renamed_load(file_obj)
