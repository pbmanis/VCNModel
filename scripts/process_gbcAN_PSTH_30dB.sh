#!/bin/bash
CONFIG="--configfile singles_multisite_parallel_PSTH.toml"
RUNTEXT="running the PSTH for RevCorr data at 30 dB"
CELLNAMES="02 05 06 09 10 11 13 17 18 30"
TABLES="data_XM13A_nacncoop_normal" # " data_XM13A_nacncoop_pasdend data_XM13A_nacncoop_actdend"
EXPERIMENT="all" # all=mean max=mean removelargest largestonly twolargest"
REPS=100
TEST="" # "--testsetup"  # or "--testsetup"
echo $RUNTEXT
for f in $CELLNAMES
do

    AXON=""
    case $f in
        02 | 05)
            AXON="-A standardized"
            ;;
    esac
    echo $f
	for t in $TABLES
    do
		for e in $EXPERIMENT
		do
		    python src/vcnmodel/model_run2.py VCN_c$f -P runANPSTH -r $REPS -D Full $AXON --Spirou $e $CONFIG --datatable $t --workers 16 $TEST
            exit_status=$?
            if [ "${exit_status}" -ne 0 ];
            then
                echo "model_run2 fails with: ${exit_status}"
                echo "Loop was: cell=${f} table=${t}  experiment=${e}"
                exit
            else
                echo "model_run2 success"
            fi
		done
	done
done

wait

echo All Variations Complete



