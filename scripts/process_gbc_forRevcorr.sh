#!/bin/bash

CELLNAMES="11" # "02 05 06 09 10 11 13 17 30"
#CONFIG="noscale.toml" #"autoscale.toml"
CONFIG="autoscale_multisite_20dB_parallel.toml"
RUNTEXT="running the AN PSTH for Revcorr"
echo $RUNTEXT
for f in $CELLNAMES
do
    echo $f
    # python vcnmodel/model_run2.py VCN_c$f  -F -P initAN --configfile $CONFIG
    python vcnmodel/model_run2.py VCN_c$f -F -P runANPSTH -r 20 --configfile $CONFIG
done


wait
echo AN PSTH for Revcorr complete
# for f in $FILES
# do
#     echo "Cell: <$f>"
#     ls -lat VCN_Cells/$f/Simulations/AN
#     echo " "
# done


