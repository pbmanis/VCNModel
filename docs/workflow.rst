Workflow for vcnmodel simulations
=================================

The simulations that can be run under vcnmodel require consideration of a large number
of parameters and stimulus conditions. Most simulations Workflow for vcnmodel simulations
will be run using the scripts
provided (in the scripts directory), which together with specific parameter files (in the toml
directory) and the Excel sheet describing the pattern for the inputs (pointed to
by the filename in the 'wheres_my_data.toml' file as the `dendriteQualityFile`). See the
`directories` entry for the structure of the various aspects of the data.

Basic steps: Morphology
-----------------------

First, the morphology of the cells must be translated from the swc format to the hoc format,
the surfaces ares adjusted, and if needed, the axons also adjusted. 

The swc file (converted to a hoc format for NEURON by `swctohoc`) for each cell must be generated
first. The hoc file that corresponds directly to the swc file will have "_Full" appended to the base
file name. 

Next, the hoc soma and dendrite areas must be "inflated" to match the surface area of the 
reconstructed cell as determined from the meshing process applied to the traced EM data. These
surface areas are stored in the dendriteQualityFile (Excel file). The inflation process requires
running the cells through the `adjust_areas` script, which will be run from the command line.
See that script for details. The script will generate an hoc file with the added name "MeshInflated". 

It is possible to replace the axons of the cells with a "standard axon". Axon measures are
summarized by the `get_axon_lengths` script. The "standard axon" is
hard coded into the `make_standardized_axon.py` script, which will take the mesh-inflated
hoc files, and replace the axon structures with the averaged axon, writing a new hoc file
with `standardized_axon` inserted into the file name.

Thus, in most of the indivicual cell Morphology directories, there will be a single SWC file,
plus several hoc files representing the successive application of these modifications. Most
simulations will be run using the "MeshInflated" version of the files, as these are the
reference morphologies.

The sizes of the auditory nerve synapses are defined in the Excel file, in the worksheet "Input_ASA_FromautoASA.py". The
script `cell_config.py` directly reads these from the Excel table, and together with the average # of synapses
per square micron (in the Synapse_Counting worksheet in the same file), is used to compute the number of active
zones to use in the `multisite` synapse mechanism that is used for each endbulb.



Basic steps: Establishing Channel Densities
-------------------------------------------

The ion channel densities are established with reference to Rothman and Manis (2003), and
Cao and Oertel (2007). The method for doing this, and the results, are shown in the manuscript.
The ion channel distributions and decoration (or assignment) to different parts of the cells
are set in the text-formatted tables in the `vcnmodel/model_data` directory. If you wish to 
explore different channel densities, these tables can be modified. The tables used for a given
simulation are called out on the command line to model_run2.py (see the scripts for examples). These tables are very
sensitive to spacing and column alignment, so they are easy to break.

Running simulations
-------------------

Simulations are directly run using `model_run2.py`, which accesses
the parameter data used for each simulation, and handles different kinds of simulations. See
the API documentation for `model_run2.py` for details of all of the command line parameters
that are available. There are a large number of commands that allow one to vary certain
experimental parameters, stimulus conditions, etc. 
See the scripts directory for how this program is called to execute different
kinds of simulations.

Analysis
--------

Most analysis and plotting is provided by the DataTablesVCN GUI tool, which is described in detail elsewhere
in this documentation.


Last update: 6 December 2021 pbm


