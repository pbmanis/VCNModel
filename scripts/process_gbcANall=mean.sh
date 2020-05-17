#!/bin/bash
CONFIG="autoscale_multisite_parallel.toml"
RUNTEXT="running the individual initialization and running AN PSTH protocols"
CELLNAMES="02 05 06 09 10 11 13 17 30"

for f in $CELLNAMES
do
    # python vcnmodel/model_run2.py $f -H -P initAN --configfile autoscale.toml
    python vcnmodel/model_run2.py VCN_c$f -H -P runANPSTH -r 10 --Spirou all=mean --configfile $CONFIG
    
done

wait

echo ANPSTH generators complete
# for f in $FILES
# do
#     echo "Cell: <$f>"
#     ls -lat VCN_Cells/$f/Simulations/AN
#     echo " "
# done


