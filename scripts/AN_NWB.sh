# Example:
# scripts/process_gbcIV.sh run all 
#

#######################################################
# Full models are from data/reconstuctions Matthew Kersting sent on
# March 6, 2020. 
# Updated 30 July 2021 
#######################################################
PROTO="-P runANPSTH"
#######################################################
# Full models are from data/reconstuctions Matthew Kersting sent on
# March 6, 2020. 
# Note we do not have a full reconstruction for cell 18
# in that dataset.
#######################################################
CELLNAMES="17"
CONFIG="--configfile xm13a_multisite_parallel.toml"
DATATABLE="--datatable data_XM13A_nacncoop_normal"
DATATABLEA="--datatable data_XM13A_nacncoop_actdend"
DATATABLEP="--datatable data_XM13A_nacncoop_pasdend"
RUNTEXT="running the individual initialization and running AN PSTH protocols"
WORKERS="16"
REPS="1"
TEST="" # "--testsetup"
SAVEALL="--saveall"
CHECK="" #"--check"  # or "" to run

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
    python vcnmodel/model_run2.py VCN_c$f -D Full -P initAN --dendritemode normal $CONFIG $AXON $DATATABLE $CHECK $TEST
    python vcnmodel/model_run2.py VCN_c$f -D Full $PROTO -r $REPS --dB 10 --saveall --Spirou all --dendritemode normal $CONFIG $AXON $DATATABLE $CHECK $TEST
    if [ $? -ne 0 ]; then
        exit 1
    fi
   #  python vcnmodel/model_run2.py VCN_c$f -D Full $PROTO -r $REPS --dB 10 --Spirou largestonly --dendritemode normal $CONFIG $AXON $DATATABLE $CHECK $TEST
   #  if [ $? -ne 0 ]; then
   #      exit 1
   #  fi
   #  python vcnmodel/model_run2.py VCN_c$f -D Full $PROTO -r $REPS --dB 10 --Spirou removelargest --dendritemode normal $CONFIG $AXON $DATATABLE $CHECK $TEST
   # if [ $? -ne 0 ]; then
   #      exit 1
   #  fi

done

echo AN runs complete
