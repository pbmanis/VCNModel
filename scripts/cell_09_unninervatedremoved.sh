# Example:
# scripts/process_gbcIV.sh run all 
#
proto="testIV"
#######################################################
# Full models are from data/reconstuctions Matthew Kersting sent on
# March 6, 2020. 
# This just computes for cell
#######################################################
CELLNAMES="09" # "02 05 06 09 10 11 13 17 18 30"
CONFIG="--configfile xm13a_multisite_parallel.toml"

MODETABLEN="--dendritemode normal --datatable data_XM13A_nacncoop_normal"
MODETABLEA="--dendritemode active --datatable data_XM13A_nacncoop_actdend"
MODETABLEP="--dendritemode passive --datatable data_XM13A_nacncoop_pasdend"
DENDMODE="NoUninnervated"

# echo "computing Zin for each configuration"
#
#
# for f in $CELLNAMES
# do
#     echo $f
#     AXON=""
#     case $f in
#         02 | 05)
#             AXON="-A standardized"
#             ;;
#     esac
#     python vcnmodel/model_run2.py VCN_c$f -P initIV $CONFIG -D $DENDMODE $AXON $MODETABLEP
#     python vcnmodel/model_run2.py VCN_c$f -P Zin  $CONFIG -D $DENDMODE  $AXON $MODETABLEP
# done
# wait
#
# for f in $CELLNAMES
# do
#     echo $f
#     AXON=""
#     case $f in
#         02 | 05)
#             AXON="-A standardized"
#             ;;
#     esac
#     python vcnmodel/model_run2.py VCN_c$f -P initIV $CONFIG -D $DENDMODE $AXON $MODETABLEA
#     python vcnmodel/model_run2.py VCN_c$f -P Zin $CONFIG -D $DENDMODE $AXON $MODETABLEA
# done
# wait
#
# for f in $CELLNAMES
# do
#     echo $f
#     AXON=""
#     case $f in
#         02 | 05)
#             AXON="-A standardized"
#             ;;
#     esac
#     python vcnmodel/model_run2.py VCN_c$f -P initIV $CONFIG -D $DENDMODE $AXON $MODETABLEN
#     python vcnmodel/model_run2.py VCN_c$f -P Zin $CONFIG -D $DENDMODE $AXON $MODETABLEN
# done
# wait
#
#
# echo Zin runs complete


PROTO="runIV"
#######################################################
# Full models are from data/reconstuctions Matthew Kersting sent on
# March 6, 2020.
#
#######################################################
CONFIG="--configfile xm13a_multisite_parallel.toml"
DATATABLE="--datatable data_XM13A_nacncoop_normal"
DATATABLEA="--datatable data_XM13A_nacncoop_actdend"
DATATABLEP="--datatable data_XM13A_nacncoop_pasdend"

# echo "running the individual initialization and/or running of IV protocols"
# for f in $CELLNAMES
# do
#     AXON=""
#     case $f in
#         02 | 05)
#             AXON="-A standardized"
#             ;;
#     esac
#     echo $f
#     echo $AXON
#     python vcnmodel/model_run2.py VCN_c$f -P initIV $CONFIG -D $DENDMODE $AXON $DATATABLE --dendritemode normal
#     python vcnmodel/model_run2.py VCN_c$f -P runIV  $CONFIG -D $DENDMODE $AXON $DATATABLE --dendritemode normal
# done
# wait
#
# for f in $CELLNAMES
# do
#     AXON=""
#     case $f in
#         02 | 05)
#             AXON="-A standardized"
#             ;;
#     esac
#     echo $f
#     echo $AXON
#     python vcnmodel/model_run2.py VCN_c$f -P initIV $CONFIG -D $DENDMODE $AXON $DATATABLEP --dendritemode passive
#     python vcnmodel/model_run2.py VCN_c$f -P runIV  $CONFIG -D $DENDMODE $AXON $DATATABLEP --dendritemode passive
# done
# wait
#
# for f in $CELLNAMES
# do
#     AXON=""
#     case $f in
#         02 | 05)
#             AXON="-A standardized"
#             ;;
#     esac
#     echo $f
#     echo $AXON
#     python vcnmodel/model_run2.py VCN_c$f -P initIV $CONFIG -D $DENDMODE $AXON $DATATABLEA --dendritemode active
#     python vcnmodel/model_run2.py VCN_c$f -P runIV  $CONFIG -D $DENDMODE $AXON $DATATABLEA --dendritemode active
# done
# wait
#
# echo IV runs complete


#############################################################
# Spike threshold test
#############################################################

CONFIG="--configfile xm13a_multisite_testSpikeThr.toml"
DATATABLE="--datatable data_XM13A_nacncoop_normal"
RUNTEXT="running the individual initialization and measuring thresholds"
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
    python vcnmodel/model_run2.py VCN_c$f -D $DENDMODE $AXON -P initIV -r 1 $CONFIG $DATATABLE
    python vcnmodel/model_run2.py VCN_c$f -D $DENDMODE $AXON -P runIVSpikeThreshold -r 1 $CONFIG $DATATABLE
    if [ $? -ne 0 ]; then
        exit 1
    fi
done

wait

echo IV threshold runs complete
mv thrrun.txt ../VCN-SBEM-Data/Cell09_Threshold_Comparisons.txt
