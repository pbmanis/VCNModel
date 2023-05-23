VCNmodel: Simulations from SBEM reconstructions of mouse bushy neurons
----------------------------------------------------------------------

*vcnmodel* is a set of Python routines in a repository for running simulations of fully reconstructed
neurons from the mouse cochlear nucleus. Neurons, including their dendrites and the portion of the axon
that remained withing the volume, afferent fibers, and afferent terminals
were reconstructed from serial blockface electron
microscopic images. The simulations include accounting for the size (number of release sites)
of multiple convergent endbulbs onto the bushy cells. Simulation modes include voltage-clamp evaluation
of currents (used to set current densities in the model), current-clamp evaluation of firing patterns
to current steps, responses to acoustic stimuli with
different patterns of input (all inputs active; only one input active, all but the largest one or two inputs
active) in response to tone pips and sinusoidal amplitude-modulated tones. In some simulations, parts of the
cell were removed (e.g., all dendrites, unninvevated dendrites, etc). In some simulations, the consequences
of variations of the dendritic ion channel densities were explored.

Used in Spirou, et al 2022/2023 (under review).

.. toctree::
   :maxdepth: 2
   :caption: Contents
   
   overview
   install
   workflow
   directories
   simfile
   scripts
   datatable
   figures
   modules

.. toctree::
   :maxdepth: 1
   :caption: Modules
   
   vcnmodel
   analysis
   morphology
   plotting
   miscanalysis
   utility
   api

Indices and Modules

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


