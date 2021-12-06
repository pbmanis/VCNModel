Installation
============

First, clone the vcnmodel repository from github.

In order to control the environment in which *vcnmodel* is run, an installation script is provided. 
This script will run under Unix-based (Posix) systems; it may need to be modified as a batch
file to run under Windows; a better solution might be to run under the Linux subsystem on a Windows machine.

The script (`make_local.sh`) will create a local environment, named `vcn_venv`. This environment
will have Python3.8, cnmodel, cochlea, and a host of standard Python modules, plus a few specific others,
needed to run the various
programs. Note that the environment also includes some packages that are only used during development, such as
black, isort, flake8. The script also compiles *cochlea* and the mechanisms used by *cnmodel*. If there are any
red lines in the print out during this setup, you will have to figure out how to fix the problem, clean the
vcn_venv directory (remove it), and run the script again. 

The install script also calls setup.py, which installs some console shortcuts::

  'console_scripts': [
       'model_run=vcnmodel.model_run2:main',
       'allgbcivs=vcnmodel.all_gbc_iv:main',
       'show_swc=scnmodel.util.show_swc:main',
       'render=vcnmodel.util.render:main',
       'plot_sims=vcnmodel.plotters.plot_sims:main',
       'datatable=vcnmodel.DataTablesVCN:main',
       'hocswcmap = vcnmodel.util.hoc_swc_sectionmap:main',
       ],

The environment should be activated with: `source vcn_venv/bin/activate`. In zsh, I use a couple of aliases to simplify matters::

    # clean deactivation - no message if there is no deactivate
    # command, otherwise, just do it.
    deact() {
        if [ "$(command -v deactivate)" ]; then
            deactivate
        fi
    }
    alias vcn="deact; cd ~/Desktop/Python/vcnmodel; source vcn_venv/bin/activate"  # do vcn model, is variant of cnmodel

In this case, typing "vcn" at the command line jumps into the directory, activates the environment, and it is ready to go.

Next, the directory where the results will be stored needs to be set up and populated with the morphology files. The
location of this directory is specified in the "wheres_my_data.toml" file at the top level of the repository. 

6 December 2021 pbm.

