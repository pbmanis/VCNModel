********************
Directory structures
********************

The workflow to run the simulations requires a specific directory structure that is *outside* the main repository. The top-level file "wheres_my_data.toml" is used
by the program to determine where to look for files.


Code
====

The code directory is vcnmodel::

    vcnmodel
            |
            docs  (sphinx documentation)
            nb  (jupyter lab notebooks)
            scripts  (shell scripts that perform batch runs of specific simulations)
            toml  (parameter files for specific simulations and path locations)
            (vcn_venv) : The pipenv that is created on installation
            vcnmodel  (sources)
                    |
                    analyzers  (code for performing specific analyses)
                    archived  (old code and data - delete)
                    model_data  (data tables used by cnmodel for specific simulations)
                    plotters  (routines for plotting results, including calling analyses)
                    simulators  (python routines that run simulations: obselete)
                    tests  (a test routine for cross-correlations)
                    util  (a variety of utility files, some of which are not used)

Data
====

The data directory holds the swc and hoc files, and the result files::

    VCN-SBEM-Data
                 |
                 VCN_Cells  (organized by cell ID number, with subdirectories, VCN_cnn)
                          |
                          VCN_c02
                                 |
                                 Initialization  (directory of initialization runs)
                                     > Files start with AN_Init_VCN_cnn, or IV_Init_VCN_cnn,
                                           followed by "fields" in the filename
                                     > Files are Python Pickled files (.p extension)
                                     > Fields are:
                                        inp=(self, or other cell)
                                        cell model (XM13A)_(sodium channel)
                                        HF=(hoc file name)
                                        multisite (or simple) - synapse type
                                        ##.# dB
                                        HS, MS, LS - spont rate 
                                  |
                                  Morphology
                                      > Stores the morphology files:
                                      > Includes the initial swc tracing, and the resulting hoc files.
                                      > Some hoc files are modifications (no dendrites, standardized axon,
                                        mesh inflated, etc). These are indicated in the file name.
                                  |
                                  Simulations
                                             |
                                             AN
                                             IV
                                         
                                             Simulations are held in subdirectories, with a filename that includes
                                             the name of the protocol (e.g., runANPSTH), a manipulation within that
                                             protocol (e.g., largestonly), and the date and time of the run. The file(s)
                                             within the directories are names similarly to the files in the Initialization
                                             directory.
                                             The simulation directories also hold small files that are "index" files into
                                             the simulations. These files hold (duplicate) the parameter blocks that are
                                             stored in the simulation files in the directories, but because they are small
                                             (they don't hold the results themselves), they can be read rapidly.
                                             When directed to do so (scan/update buttons), the DataTablesVCN module 
                                             examines the simulation directories and writes the index files (.pkl extension)
                                             as needed. The information in the index files is displayed by DataTablesVCN,
                                             for access to specific simulation results, plotting, analysis, and figure
                                             generation.
                            |
                            (tree structure is repeated for each cell)