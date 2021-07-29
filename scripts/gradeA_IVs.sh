# Example:
# scripts/process_gbcIV.sh run all 
#
proto="runIV"
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

echo "running the individual initialization and/or running of IV protocols"
for f in $CELLNAMES
do
    AXON=""
    case $f in
        02 | 05)
            AXON="-A standardized"
            ;;
    esac
    echo $f
    echo $AXON
    python src/vcnmodel/model_run2.py VCN_c$f -P initIV $CONFIG -D Full $AXON $DATATABLE --dendritemode normal
    python src/vcnmodel/model_run2.py VCN_c$f -P runIV  $CONFIG -D Full $AXON $DATATABLE --dendritemode normal
done
wait

for f in $CELLNAMES
do
    AXON=""
    case $f in
        02 | 05)
            AXON="-A standardized"
            ;;
    esac
    echo $f
    echo $AXON
    python src/vcnmodel/model_run2.py VCN_c$f -P initIV $CONFIG -D Full $AXON $DATATABLEP --dendritemode passive
    python src/vcnmodel/model_run2.py VCN_c$f -P runIV  $CONFIG -D Full $AXON $DATATABLEP --dendritemode passive
done
wait

for f in $CELLNAMES
do
    AXON=""
    case $f in
        02 | 05)
            AXON="-A standardized"
            ;;
    esac
    echo $f
    echo $AXON
    python src/vcnmodel/model_run2.py VCN_c$f -P initIV $CONFIG -D Full $AXON $DATATABLEA --dendritemode active
    python src/vcnmodel/model_run2.py VCN_c$f -P runIV  $CONFIG -D Full $AXON $DATATABLEA --dendritemode active
done
wait

echo IV runs complete
