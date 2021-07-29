#!/bin/bash

CELLNAMES="02 05 06 09 10 11 13 17 18 30"
#CONFIG="noscale.toml" #"autoscale.toml"
CONFIG="xm13a_multisite_parallel.toml"
RUNTEXT="running the AN Single protocols"
TABLES="data_XM13A_nacncoop_normal data_XM13A_nacncoop_pasdend data_XM13A_nacncoop_actdend"

echo $RUNTEXT
for t in $TABLES
    do
    for f in $CELLNAMES
        do
        echo $f
        python src/vcnmodel/model_run2.py VCN_c$f  -D Full -P initAN --configfile $CONFIG  --datatable $t
        python src/vcnmodel/model_run2.py VCN_c$f  -D Full -P runANSingles -r 5 --configfile $CONFIG  --datatable $t
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


