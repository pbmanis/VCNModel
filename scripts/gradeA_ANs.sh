# Example:
# scripts/process_gbcIV.sh run all 
#

#######################################################
# Full models are from data/reconstuctions Matthew Kersting sent on
# March 6, 2020. 
# 
#######################################################
PROTO="runANPSTH"
#######################################################
# Full models are from data/reconstuctions Matthew Kersting sent on
# March 6, 2020. 
# Note we do not have a full reconstruction for cell 18
# in that dataset.
#######################################################
CELLNAMES="02 05 06 09 10 11 13 17 18 30"
CONFIG="--configfile xm13a_multisite_parallel.toml"
DATATABLE="--datatable data_XM13A_nacncoop_normal"
DATATABLEA="--datatable data_XM13A_nacncoop_actdend"
DATATABLEP="--datatable data_XM13A_nacncoop_pasdend"
CONFIG="autoscale_xm13a_multisite_parallel.toml"
RUNTEXT="running the individual initialization and running AN PSTH protocols"
WORKERS="16"
REPS="100"
CHECK="--check"  # or "" to run

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
    echo $f
    echo $AXON
    python vcnmodel/model_run2.py VCN_c$f  -F -P initAN $CONFIG $AXON $DATATABLE $CHECK
    python vcnmodel/model_run2.py VCN_c$f  -F $PROTO -r $REPS --dB 10 --Spirou all $CONFIG $AXON $DATATABLE $CHECK
    if [ $? -ne 0 ]; then
        exit 1
    fi
    python vcnmodel/model_run2.py VCN_c$f  -F $PROTO -r $REPS --dB 10 --Spirou largestonly $CONFIG $AXON $DATATABLE $CHECK
    if [ $? -ne 0 ]; then
        exit 1
    fi
    python vcnmodel/model_run2.py VCN_c$f  -F $PROTO -r $REPS --dB 10 --Spirou removelargest $CONFIG $AXON $DATATABLE $CHECK
   if [ $? -ne 0 ]; then
        exit 1
    fi

done

echo AN runs complete
