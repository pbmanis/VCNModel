import io
import pickle
# import pandas
"""
Rename import paths for pickle load after restructuring data 

From stack overflow:
https://stackoverflow.com/questions/2121874/python-pickling-after-changing-a-modules-directory

"""
class RenameUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        renamed_module = module
        # if module.startswith('vcnmodel'):
        #     print('module: ', module)
        if module in ["vcnmodel","vcnmodel.model_params", 'vcnmodel.table_manager',
            "vcnmodel.plotters.plot_sims"]:
            renamed_module = "src." + module
        # print(renamed_module)
        return super(RenameUnpickler, self).find_class(renamed_module, name)

# class RenameUnpickler_pandas(pandas.read_pickle):
#     def find_class(self, module, name):
#         renamed_module = module
#         if module.startswith('vcnmodel'):
#             print('module: ', module)
#         if module == "vcnmodel":
#             renamed_module = "src.vcnmodel"
#         if module == "vcnmodel.model_params":
#             renamed_module = "src.vcnmodel.model_params"
#         elif module == 'vcnmodel.table_manager':
#             renamed_module = "src.vcnmodel.table_manager"
#
#         return super(RenameUnpickler_pandas, self).find_class(renamed_module, name)
#

def pickle_load(file_obj):
    return RenameUnpickler(file_obj).load()


def pickle_loads(pickled_bytes):
    file_obj = io.BytesIO(pickled_bytes)
    return renamed_load(file_obj)

# def pd_read_pickle(file_obj):
#     return RenamUnpickler_pandas(file_obj).read_pickle
#