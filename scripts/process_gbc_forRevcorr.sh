#!/bin/bash

CELLNAMES="02 05 06 09 10 11 13 17 30"
#CONFIG="noscale.toml" #"autoscale.toml"
CONFIG="autoscale_multisite_10dB_parallel.toml"
RUNTEXT="running the AN PSTH for Revcorr"
TABLES="data_XM13A_nacncoop_pasdend data_XM13A_nacncoop_actdend"

echo $RUNTEXT
for t in $TABLES
    do
    echo $t
    for f in $CELLNAMES
    do
        echo $f
        # python vcnmodel/model_run2.py VCN_c$f  -F -P initAN --configfile $CONFIG
        python vcnmodel/model_run2.py VCN_c$f -F -P runANPSTH -r 25 --configfile $CONFIG --datatable $t
    done
done

wait
echo AN PSTH for Revcorr complete
# for f in $FILES
# do
#     echo "Cell: <$f>"
#     ls -lat VCN_Cells/$f/Simulations/AN
#     echo " "
# done


