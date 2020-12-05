#!/bin/bash
CONFIG="toml/autoscale_multisite_40dB_parallel.toml"
RUNTEXT="running the individual initialization and running AN PSTH protocols"
CELLNAMES="02 05 06 09 10 11 13 17 18 30"
TABLES="data_XM13A_nacncoop_normal" # " data_XM13A_nacncoop_pasdend data_XM13A_nacncoop_actdend"
EXPERIMENT="all" # all=mean max=mean removelargest largestonly twolargest"
REPS=25
TEST="" #"--testsetup"  # or "--testsetup"
echo $RUNTEXT
for f in $CELLNAMES
do
	for t in $TABLES
    do
		for e in $EXPERIMENT
		do
		    python vcnmodel/model_run2.py VCN_c$f -P runANPSTH -r $REPS -D Full --Spirou $e --configfile $CONFIG --datatable $t --workers 8 $TEST
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



