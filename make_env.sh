# Takes one command line argument "local".
# if "local", then we use the local sources rather than calling down
# the sources from github.
# This does not apply to cochlea, which needs to be pulled from the upstream
# on pbmanis with the correct commit. 

# define the "local" environment path. The required local sources
# must be under this subdirectory. Modify as needed

LOCAL="/Users/pbmanis/Desktop/Python/"

set -e # force failure if anyting fails in script - ensuring completion
set -o errexit
ENVNAME="vcn_venv"
if [ -d $ENVNAME ]
then
    echo "Removing previous environment: $ENVNAME"
    rm -Rf $ENVNAME
    # set +e
    # rsync -aR --remove-source-files $ENVNAME ~/.Trash/ || exit 1
    # set -e
    # rm -R $ENVNAME
else
    echo "No previous environment - ok to proceed"
fi

python3.10 -m venv $ENVNAME || exit 1

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
pip3 install -r requirements.txt || exit 1
source $ENVNAME/bin/activate
pip3 install mayavi
# build the mechanisms
# this may require a separate install of the standard NEURON package
# with the same version as we have provided
nrnivmodl ../cnmodel/cnmodel/mechanisms

# these has to be done separately as the compilation depends on what happens above
pip3 install -e git+https://github.com/pbmanis/cochlea-1.git@8152f032e6e619d8632548b6b632fcda5f0638ed#egg=cochlea
rm -rf cochleae

if [[$0 == "local"]]
then
    pip install -e $LOCAL"cnmodel"
    pip install -e $LOCAL"pylibrary"
    pip install -e $LOCAL"ephys"
    pip install -e $LOCAL"montage"
    pip install -e $LOCAL"neuronvis"
else
    pip install -e git+https://github.com/pbmanis/cnmodel.git#egg=cnmodel
    pip install -e git+https://github.com/pbmanis/pylibrary.git@main#egg=pylibrary
    pip install -e git+https://github.com/pbmanis/ephys.git@022edf3885bb64328a84b30b1d83aa8f1871cae0#egg=ephys
    pip install -e git+https://github.com/pbmanis/montager.git#egg=montage
    pip install -e git+https://github.com/pbmanis/neuronvis@179cef80314ce88eb8e352200be8fb6f02ffa6e1#egg=neuronvis
 fi
source $ENVNAME/bin/activate
python setup.py develop

echo "Success in installing vcnmodel environment!"
