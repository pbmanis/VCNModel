# Example:
# scripts/process_gbcIV.sh run all 
#

#######################################################
# Full models are from data/reconstuctions Matthew Kersting sent on
# March 6, 2020. 
# Note we do not have a full reconstruction for cell 18
# in that dataset.
# This script runs to get the thresholds for current injection for one
# spike. 
#######################################################
CELLNAMES="02 05 06 09 10 11 13 17 18 30" # 
CELLNAMES="18"
#CONFIG="noscale.toml" #"autoscale.toml"
CONFIG="autoscale_multisite_0dB_parallel_synthr.toml"
RUNTEXT="running the individual initialization and running AN PSTH protocols"
AXON="default" # or "standardized"
echo $RUNTEXT
for f in $CELLNAMES
do
    echo $f

    python src/vcnmodel/model_run2.py VCN_c$f  -D Full -A $AXON -P initIV -r 1 --configfile $CONFIG --datatable data_XM13A_nacncoop
    python src/vcnmodel/model_run2.py VCN_c$f  -D Full -A $AXON -P runIVSpikeThreshold -r 1 --configfile $CONFIG --datatable data_XM13A_nacncoop
    if [ $? -ne 0 ]; then
        exit 1
    fi
done

wait

echo IV threshold runs complete
# with "A", we do all cells in grade A
# python vcnmodel/plotters/plot_gbc_ivresults_2.py "A" -s