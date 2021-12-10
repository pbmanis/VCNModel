*************************
DataTables
*************************

.. figure:: _static/DataTablesVCN.png
  :width: 720
  :align: center
  :alt: Image of DataTableVCN GUI
 
Th `vcnmodel` project is capable of generating a large number of simulation files, which require
an approach to keep track of which simulations go together, along with a way of 
easily showing the parameters used for different types of
simulations. DataTablesVCN (`datatables` command) is a GUI-based helper that can
provide a listing of the data from each cell, including filtering and sorting the
simulation results. It also provides access to analysis and plotting routines,
as well as figure generation.

**NOTE: DataTablesVCN does not control simulation runs!**

DataTablesVCN consists of two primary scripts,
DataTablesVCN.py and table_manager.py, which together manage the table interface,
collect and display simulation entries, and have access to the plotting and analysis
programs (largely in plot_sims.py and figures.py). DataTablesVCN is (of course) written
in Python, and makes extensive use of pyqtgraph's parametertree for user interaction,
and docking system for display of tables, graphic data and text.

GUI layout
----------

DataTables has a pane on the left hand side that provides the controls for the program (titled "Params").
THe main tabular display on the right ("Simulation Table") displays information about the (selected) simulation runs 
of a given type (AN, IV) for a given cell selected in the Cells drop-down list. Below the main table is a "Reporting"
window that displays some of the output from the program (indicating what the program is doing, errors, and sometimes
text that can be copied to Excel, Prism, or a Python script for further analysis or manual display). 

The main pane is actually a "dock", with two windows. The second window can be accessed by clicking on the dock title bar that 
says "Traces", and the main table can be viewed again by clicking on the dock title "Simulation Table". The "Traces"
dock window is only used by a specific part of the program, as indicated below.

The Params control window is organized into subgroups:

Models
^^^^^^

The `Models` group has two buttons.

    `Scans`
        Clicking this button causes the program to scan the existing simulation runs to populate the Simulation Table.
        Only runs that have already been abstracted are shown.
        Selection of a new cell (see `Selections`) also runs the scan to populate the table

    `Update Runs`
        The second button, "Update Runs" also populates the table, but in addition deterimes 
        which runs have not been abstracted yet, and will generate the file with that 
        information before updating the table.

Selections
^^^^^^^^^^

The `Selections` group helps to select which data sets will be displayed in the table, through drop-down lists. The following
actions are provided:

    `Run Type`
        Selects a run type either AN (auditory nerve; acoustic stimulation), IV (current clamp), or VC (voltage clamp).

    `Cells`
        This drop-down allows the user to select a specific cell for display in the table.
        
    `Start Date`
        Selects a starting date for the data that will be displayed. If *None* is selected,
        then there is no defined start date and the oldest simulation will be included in the table.
        Note that all runs are tagged by date and time.

    `End Date`
        Sets the end date for the data that will be displayed. If *None* is selected, the there is no end date, and
        the most recent simulation will be included in the table.

    `ModelType`
        Selects the *modelName* for models to be included in the Simulation Table. These are usually differing decoration
        patterns for ion channels.

    `Mode`
        Not sure what this does; I never use it.

    `Experiment`
        Selects the *Experiment* type that will be displayed. Experiments correspond to manipulations of either the 
        pattern of inputs, or changes to the patterns of the inputs (e.g., using inputs that are all the same
        strength is selected by choosing `mean`).

    `Analysis`
    
    `Dendrites`
        Displays only simulations with a specific type of dendrite decoration.

Analysis
^^^^^^^^

The `Analysis` group performs specific analyses on the data that is selected in the table. 

    `Traces`
        Plots all of the traces in the selected run on top of each other

    `IV` or `VC`
        If the selected dataset is an IV or VC run, computes and plots an IV or VC Figure. Otherwise it just plots the traces.

    `Singles`
        If the selected dataset is a "Singles" run (e.g., individual inputs were tested independently), this will plot
        a stacked set of traces with action potentials marked.
        This will also print out (in the Reporting pane) the efficacy (output spikes/input spikes) for each synapse in the
        run, along with some ancillary information. This table can be copied.

    `Revcorr X`
        The 3 revcorr buttons compute the reverse correlation using slightly different methods. RevcorrSPKS uses Brian SImulator
        revcorr calculation. RevcorrSimple does a simple (and probably inefficient) revcorr computation. RevcorrSTTC would
        compute the spike time tiling correlation, but is not implemented.

    `PSTH`
        For runs with acoustic stimuli, this will display a figure with a voltage trace, a raster plot of the postsynaptic
        spikes, a PSTH bar plot, the stimulus waveform, a phase histogram, the ANF spike trains for each input (one trial),
        the ANF PSTH, and finally, a first/second spike latency histogram.

Filters
^^^^^^^

The `Filters` group permits filtering the displayed data set based on a number of different simulation
parmeters. Multiple filters are applied in sequence. Fiters are not applited until the `Apply` button in the `Filter Actions`
line is clicked. Filters are removed by clicking the `Clear` button.

Options
^^^^^^^

Figures
^^^^^^^

Tools
^^^^^
These are various tools that no other home.
    
    `Reload`
        Forces reloading of the Python code from a select set of modules. Normally used only during program development.

    `View IndexFile`
        Prints (to the terminal) nicely formatted, and way too much, information from the selected simulation. This is
        basically a text output from the .pkl file that is generated from the directory that holds the simulation.

    `Print File Info`
        This prints, to the "Reporting" pane, a short text that represents the python dict that points to this file (for
        use in the figure_data program).

    `Delete Selected Sim`
        Sometimes you just don't need to keep a simulation - the data was corrupted, or the run is not useful, or it is a 
        duplication. This button lets you actually remove the simulation data and folder from the disk. Use with care.


7 December 2021 pbm



