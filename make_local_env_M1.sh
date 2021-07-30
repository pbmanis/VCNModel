set -e # force failure if anyting fails in script - ensuring completion
set -o errexit
ENVNAME="vcn_venv"
if [ -d $ENVNAME ]
then
    echo "Removing previous environment: $ENVNAME"
    set +e
    rsync -aR --remove-source-files $ENVNAME ~/.Trash/ || exit 1
    set -e
    rm -R $ENVNAME
else
    echo "No previous environment - ok to proceed"
fi

python3.8 -m venv $ENVNAME || exit 1

source $ENVNAME/bin/activate || exit 1

pip3 install --upgrade pip  # be sure pip is up to date in the new env.
pip3 install wheel  # seems to be missing (note singular)
pip3 install cython
# # if requirements.txt is not present, create:
# # pip install pipreqs
# # pipreqs
#
# #Then:
#
pip3 install -r requirements_local_M1.txt || exit 1
source $ENVNAME/bin/activate
pip3 install mayavi
# build the mechanisms
# this may equire a separate install of the standard NEURON package
# with the same version as we have provided
nrnivmodl ../../cnmodel/mechanisms

# these has to be done separately as the compilation depends on what happens above
pip3 install -e git+https://github.com/pbmanis/cochlea-1.git@c2e8c9612481ebe397ba0e7762b8f2772500388d#egg=cochlea
rm -rf cochleae
pip3 install -e git+https://github.com/pbmanis/neuronvis.git@Python3#egg=neuronvis
source $ENVNAME/bin/activate
python setup.py develop

echo "Success in installing acq4 environment!"
