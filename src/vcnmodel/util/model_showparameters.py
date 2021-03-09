from collections import OrderedDict
import pickle
import pylibrary as PL
import pprint
import src.vcnmodel.model_params as model_params


""""
Show the parameters for a given model run
"""

MP = model_params.ModelParams()
pp = pprint.PrettyPrinter()

fn = 'VCN_Cells/VCN_c09/Simulations/IV/VCN_c09_pulse_XM13_nacncoop_II_soma=1.514_dend=1.514_ASA=1.514_monitor.p'
def show_model_params(filename: Union[Path, str]):
    fn = Path(filename)
    with open(fn, 'rb'):
        d = pickle.load(open(fn, 'rb'))
    print(d.keys())
    if 'Paramsx' in list(d.keys()):
        pp.pprint(d['Params'])

    else:
        # fill params from filename and modelpars
        print('\nmodelPars: ')
        pp.pprint(d['modelPars'])
        print()
        print('\nrunInfo: ')
        pp.pprint(d['runInfo'])
        print()
        print('\nDefault Params: ')
        pp.pprint(MP.Params)
    
    Par = MP.Params.copy()
    modelPars = d['modelPars']
    runInfo = d['runInfo']

    pk = list(Par.keys())
    for mpk in list(modelPars.keys()):
        if mpk in pk:
            if Par[mpk] == Par[mpk]:
                print(f"modelpar Duplicate key/value: ", mpk, Par[mpk], Par[mpk])
            else:
                print(f"modelpar Matched key, different value: ", mpk, Par[mpk], Par[mpk])
        else:
            print(f"Shoud add modelpar key {mpk:s} to Par keys")

    for mpk in list(runInfo.keys()):
        if mpk in pk:
            if Par[mpk] == Par[mpk]:
                print(f"runInfo Duplicate key/value: ", mpk, Par[mpk], Par[mpk])
            else:
                print(f"runInfo Matched key, different value: ", mpk, Par[mpk], Par[mpk])
        else:
            print(f"Shoud add runInfo key {mpk:s} to Par keys")
