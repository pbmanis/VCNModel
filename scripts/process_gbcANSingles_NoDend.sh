#!/bin/bash

CELLNAMES="02 05 06 09 10 11 13 17 30"
CELLNAMES="09 11 13 30"
#CONFIG="noscale.toml" #"autoscale.toml"
CONFIG="singles_autoscale_multisite_parallel_test.toml"
RUNTEXT="running the AN Single protocols NoDend"
TABLES="data_XM13A_nacncoop_normal" # " data_XM13A_nacncoop_pasdend data_XM13A_nacncoop_actdend"
DEND="NoDend"

echo $RUNTEXT
for t in $TABLES
    do
    for f in $CELLNAMES
        do
        echo $f
        python vcnmodel/model_run2.py VCN_c$f -P initAN -D $DEND --configfile $CONFIG  --datatable $t
        python vcnmodel/model_run2.py VCN_c$f -P runANSingles -r 5 -D $DEND --configfile $CONFIG  --datatable $t
     done
done

wait
echo AN Singles complete
# for f in $FILES
# do
#     echo "Cell: <$f>"
#     ls -lat VCN_Cells/$f/Simulations/AN
#     echo " "
# done


