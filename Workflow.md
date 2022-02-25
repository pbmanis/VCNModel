Workflow for vcnmodel simulations =================================

Last update: 16 Aug 2021 pbm

Simulations are controlled by shell scripts in the scripts subdirectory. These scripts set up the parameters for the
runs, including the cnmodel datatables and various stimulus parameters. The results are stored on a by-cell basis in
the target directory as "VCN-SBEM-Data/VCN_Cells/VCN_Cnm/IV (or AN)" with a long filename. The scripts invoke
model_run2.py, which handles all of the simulations and "experiments", wrapped around cnmodel. The data files contain
all of the information regarding the model run, as stored in the dataclasses Params and RunInfo (there should be no
control parameters that are not described in these tables). This includes the full text of the "data\_table" that is
used ton control cell decoration of channels in cnmodel.

Analysis is handled through "DataTablesVCN.py", which has a graphical interface. This is described in more detail below.

Figures -------

For some figures, the data to be plotted is collected manualy using the "print file info" in DataTablesVCN.py, under
the Tools bar. The results of the analysis or the files selected are printed into the lower "reporting" window on the
right, where they may be copied out to the appropriate tables in either "vcnmodel/plotting/figure_data.py" or
"vs_datasets.py"

Figure files are sent in PDF format to the "VCN-SBEM-Data/Figures" directory. The numbering is sequential, and parts of
these individual files are reassembled in Illustrator to generate the main figures; others are summaries that are
presented as supplemental data.
