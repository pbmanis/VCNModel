VCNModel
========

VCNModel is a set of Python routines in a repository for running simulations of fully reconstructed neurons from the mouse cochlear nucleus.
Neurons, including their dendrites and the portion of the axon that remained withing the volume, afferent fibers, 
and afferent terminals were reconstructed from serial blockface electron microscopic images. 
The simulations include accounting for the size (number of release sites) of multiple convergent endbulbs onto the bushy cells.
Simulation modes include voltage-clamp evaluation of currents (used to set current densities in the model), current-clamp evaluation of firing patterns to current steps, responses to acoustic stimuli with different patterns of input (all inputs active; only one input active, all but the largest one or two inputs active) in response to tone pips and sinusoidal amplitude-modulated tones. In some simulations, parts of the cell were removed (e.g., all dendrites, unninvevated dendrites, etc). In some simulations, the consequences of variations of the dendritic ion channel densities were explored.

Used in Spirou et al. (in preparation, 2022).

Installation
------------

First, clone the vcnmodel repository from github.

**Installation Script**

In order to control the environment in which vcnmodel is run, an installation script is provided. The script (make_local_env.sh) will create a local environment, named vcn_venv, in the top level of the repository. This script will run under Unix-based (Posix) systems; it may need to be modified as a batch file to run under Windows; a better solution might be to run under the Linux subsystem on a Windows machine.

The environment that is created will have Python3.9, cnmodel, cochlea, and a host of standard Python modules, plus a few specific others, needed to run the various programs. Note that the environment also includes some packages that are only used during development, such as black, isort, flake8. The script also compiles cochlea and the mechanisms used by cnmodel. If there are any red lines in the print out during this setup, you will have to figure out how to fix the problem, clean the vcn_venv directory (remove it), and run the script again. The install script is for use on macOS and Linux systems. No Windows-based batch file is provided at this time.

The install script also calls setup.py, which installs some console shortcuts:

```
'console_scripts': [
     'model_run=vcnmodel.model_run2:main',
     'allgbcivs=vcnmodel.all_gbc_iv:main',
     'show_swc=scnmodel.util.show_swc:main',
     'render=vcnmodel.util.render:main',
     'plot_sims=vcnmodel.plotters.plot_sims:main',
     'datatable=vcnmodel.DataTablesVCN:main',
     'hocswcmap = vcnmodel.util.hoc_swc_sectionmap:main',
     ],
```


Environment
-----------
The environment should be activated with: source vcn_venv/bin/activate. In zsh, a couple of aliases will simplify matters:

  ```
  # clean deactivation - no message is printed if there is no deactivate
  # command, otherwise, just do it.

  deact() {
      if [ "$(command -v deactivate)" ]; then
          deactivate
      fi
  }
  alias vcn="deact; cd ~/Desktop/Python/vcnmodel; source vcn_venv/bin/activate"  # do vcn model, is variant of cnmodel.
  ```
  
  Once this has been set up, typing “vcn” at the command line will change into the vcnmodel directory and activate the environment.
  At that point, it is ready to run the shell scripts or Python scripts.

Documentation
-------------
At this point, it is recommended to generate and read the documentation:

  ```
  cd /docs
  make html
  open html/index.html
  ```

Plese also read the document Workflow.md in this directory.



