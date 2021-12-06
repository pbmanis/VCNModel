********
Overview
********

`vnmodel` relies heavily on 3 other packages:
    * `cnmodel` (Manis and Campagnola, 2014), which provides the code that constructs
      cells from the hoc files, and decorates their
      different compartments with ion channels according to text-based tables.
    * The auditory nerve model of Zilany et al. (2014), accessed via the package `cochlea` (Rudnicki
      and Hemmert, 2017).
    * `NEURON` (Hines and Carnevale, 20xx); www.neuron.yale.edu.
    
`vcnmodel` is primarily about organizing the simulations, exploring parameter spaces, and performing
data analysis and plotting.

The `vcnmodel` repository contains Python routines for:
    * Configuring simulations and cell morphology (Morphology modules).
    * Controlling and running simulations (Simulation modules).
    * Evaluating and analyzing data (Analysis modules)
    * Plotting (documenting) simulation results, and generating figures (Plotting modules)
    * Miscellaneous modules for data management.

The repository does not hold the primary data that drives the simulations, nor simulation results, with the exception of the
tables that are read into *cnmodel*, which are held in the vcnmodel/data_tables subdirectory.







