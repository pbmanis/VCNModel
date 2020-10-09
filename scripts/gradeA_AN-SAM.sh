# Example:
# scripts/process_gbcIV.sh run all 
#

#######################################################
# Full models are from data/reconstuctions Matthew Kersting sent on
# March 6, 2020. 
# Note we do not have a full reconstruction for cell 18
# in that dataset.
#######################################################
CELLNAMES="05 06 10" # "02 05 06 09 10 11 13 17 30"
#CONFIG="noscale.toml" #"autoscale.toml"
CONFIG="autoscale_multisite_SAM.toml"
RUNTEXT="running the individual initialization and running AN SAM protocols"
WORKERS="16"
REPS="100"
echo $RUNTEXT

for f in $CELLNAMES
do
    echo $f
#   echo $f
    python vcnmodel/model_run2.py VCN_c$f  -F -P initAN --configfile $CONFIG  --datatable data_XM13A_nacncoop
    python vcnmodel/model_run2.py VCN_c$f  -F -P runANPSTH -r $REPS --dB 20 --Spirou all --workers $WORKERS --configfile $CONFIG --datatable data_XM13A_nacncoop
    if [ $? -ne 0 ]; then
        exit 1
    fi
    python vcnmodel/model_run2.py VCN_c$f  -F -P runANPSTH -r $REPS --dB 20 --Spirou largestonly --workers $WORKERS --configfile $CONFIG --datatable data_XM13A_nacncoop
	if [ $? -ne 0 ]; then
	    exit 1
	fi
	python vcnmodel/model_run2.py VCN_c$f  -F -P runANPSTH -r $REPS --dB 20 --Spirou removelargest --workers $WORKERS --configfile $CONFIG --datatable data_XM13A_nacncoop
	if [ $? -ne 0 ]; then
	    exit 1
	fi
done


wait

echo AN SAM runs complete
# with "A", we do all cells in grade A
# python vcnmodel/plotters/plot_gbc_ivresults_2.py "A" -s
