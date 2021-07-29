# Example:
# scripts/test_IVs_thr.sh run all 
#

#######################################################
# Full models are from data/reconstuctions Matthew Kersting sent on
# March 6, 2020. 
# Note we do not have a full reconstruction for cell 18
# in that dataset.
# This script runs to get the thresholds for current injection for one
# spike. 
#######################################################
CELLNAMES="02 05 06 09 10 11 13 17 18 30"
#CELLNAMES="05"

CONFIG="--configfile xm13a_multisite_testSpikeThr.toml"
DATATABLE="--datatable data_XM13A_nacncoop_normal"
RUNTEXT="running the individual initialization and running AN PSTH protocols"
# AXON="default" # or "standardized"
echo $RUNTEXT
for f in $CELLNAMES
do
    echo $f
    AXON=""
    case $f in
        02 | 05)
            AXON="-A standardized"
            ;;
    esac
    python src/vcnmodel/model_run2.py VCN_c$f  -D Full $AXON -P initIV -r 1 $CONFIG $DATATABLE
    python src/vcnmodel/model_run2.py VCN_c$f  -D Full $AXON -P runIVSpikeThreshold -r 1 $CONFIG $DATATABLE
    if [ $? -ne 0 ]; then
        exit 1
    fi
done

wait

echo IV threshold runs complete
# with "A", we do all cells in grade A
# python vcnmodel/plotters/plot_gbc_ivresults_2.py "A" -s