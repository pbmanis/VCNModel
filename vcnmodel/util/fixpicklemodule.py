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
        #  Previously had modules one deeper (src/vcnmodel); now just at top
        if module in ["src.vcnmodel","src.vcnmodel.model_params", 'src.vcnmodel.table_manager',
            "src.vcnmodel.plotters.plot_sims", "src.vcnmodel.plotters.figures", 
            "src.vcnmodel.util"]:
            renamed_module =  module  # used to be "src." + module
        if module.startswith("src."):
            renamed_module = module[4:]  # so strip src. from module name
        # print("renamed module: ", renamed_module)
        return super(RenameUnpickler, self).find_class(renamed_module, name,)


def pickle_load(file_obj):
    return RenameUnpickler(file_obj).load()


def pickle_loads(pickled_bytes):
    #file_obj = io.BytesIO(pickled_bytes)
    file_obj = dill.loads(pickled_bytes)
    return renamed_load(file_obj)
