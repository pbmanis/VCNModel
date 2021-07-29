# Example:
# scripts/process_gbcIV.sh run all 
#
proto="testIV"
#######################################################
# Full models are from data/reconstuctions Matthew Kersting sent on
# March 6, 2020. 
# Note we do not have a full reconstruction for cell 18
# in that dataset.
#######################################################
CELLNAMES="02 05 06 09 10 11 13 17 18 30" # "02 05 06 09 10 11 13 17 18 30"
CONFIG="--configfile xm13a_multisite_parallel.toml"

MODETABLEN="--dendritemode normal --datatable data_XM13A_nacncoop"
MODETABLEA="--dendritemode active --datatable data_XM13A_nacncoop_actdend"
MODETABLEP="--dendritemode passive --datatable data_XM13A_nacncoop_pasdend"

echo "computing Zin for each gradeA Cell"


# for f in $CELLNAMES
# do
#     echo $f
#     AXON=""
#     case $f in
#         02 | 05)
#             AXON="-A standardized"
#             ;;
#     esac
#     python src/vcnmodel/model_run2.py VCN_c$f -P initIV $CONFIG -D Full $AXON $MODETABLEN
#     python src/vcnmodel/model_run2.py VCN_c$f -P Zin $CONFIG -D Full $AXON $MODETABLEN
# done
# wait

# for f in $CELLNAMES
# do
#     echo $f
#     AXON=""
#     case $f in
#         02 | 05)
#             AXON="-A standardized"
#             ;;
#     esac
#     python src/vcnmodel/model_run2.py VCN_c$f -P initIV $CONFIG -D Full $AXON $MODETABLEP
#     python src/vcnmodel/model_run2.py VCN_c$f -P Zin  $CONFIG -D Full  $AXON $MODETABLEP
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
#     python src/vcnmodel/model_run2.py VCN_c$f -P initIV $CONFIG -D Full $AXON $MODETABLEA
#     python src/vcnmodel/model_run2.py VCN_c$f -P Zin $CONFIG -D Full $AXON $MODETABLEA
# done
# wait

for f in $CELLNAMES
do
    echo $f
    AXON=""
    case $f in
        02 | 05)
            AXON="-A standardized"
            ;;
    esac
    python src/vcnmodel/model_run2.py VCN_c$f -P initIV $CONFIG -D NoDend $AXON $MODETABLEN
    python src/vcnmodel/model_run2.py VCN_c$f -P Zin $CONFIG -D NoDend $AXON $MODETABLEN
done
wait

for f in $CELLNAMES
do
    echo $f
    AXON=""
    case $f in
        02 | 05)
            AXON="-A standardized"
            ;;
    esac
    python src/vcnmodel/model_run2.py VCN_c$f -P initIV $CONFIG -D AxonOnly $AXON $MODETABLEN
    python src/vcnmodel/model_run2.py VCN_c$f -P Zin $CONFIG -D AxonOnly $AXON $MODETABLEN
done
wait

echo Zin runs complete