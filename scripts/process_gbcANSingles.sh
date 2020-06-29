#!/bin/bash

CELLNAMES="02 05 06 09 10 11 13 17 30"
#CONFIG="noscale.toml" #"autoscale.toml"
CONFIG="singles_autoscale_multisite_noparallel.toml"
RUNTEXT="running the AN Single protocols"
echo $RUNTEXT
for f in $CELLNAMES
do
    echo $f
    # python vcnmodel/model_run2.py VCN_c$f  -F -P initAN --configfile $CONFIG
    python vcnmodel/model_run2.py VCN_c$f  -F -P runANSingles -r 25 --configfile $CONFIG
done


wait
echo AN Singles complete
# for f in $FILES
# do
#     echo "Cell: <$f>"
#     ls -lat VCN_Cells/$f/Simulations/AN
#     echo " "
# done


