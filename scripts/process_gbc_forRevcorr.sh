#!/bin/bash

CELLNAMES="02" # " 05 06 09 10 11 13 17 18 30"
#CONFIG="noscale.toml" #"autoscale.toml"
CONFIG="toml/autoscale_multisite_10dB_parallel.toml"
RUNTEXT="running the AN PSTH for Revcorr"
TABLES="data_XM13A_nacncoop_normal data_XM13A_nacncoop_pasdend data_XM13A_nacncoop_actdend"

echo $RUNTEXT
# for t in $TABLES
#     do
#     echo $t
#     for f in $CELLNAMES
#     do
#         echo $f_
#         python vcnmodel/model_run2.py VCN_c$f  -D Full-P initAN --dendriteMode --configfile $CONFIG --datatable $t
#         python vcnmodel/model_run2.py VCN_c$f -D Full -P runANPSTH -r 25 --configfile $CONFIG --datatable $t
#     done
# done
#
# wait
# echo AN PSTH for Revcorr complete
# for f in $FILES
# do
#     echo "Cell: <$f>"
#     ls -lat VCN_Cells/$f/Simulations/AN
#     echo " "
# done

echo "running the individual initialization and/or running of IV protocols"
for f in $CELLNAMES
do
    echo $f
    python vcnmodel/model_run2.py VCN_c$f -P initIV --configfile $CONFIG -D Full --dendritemode normal --datatable data_XM13A_nacncoop
    python vcnmodel/model_run2.py VCN_c$f -P runANPSTH  -r 25 --configfile $CONFIG -D Full --dendritemode normal --datatable data_XM13A_nacncoop
done

wait

for f in $CELLNAMES
do
    echo $f
    python vcnmodel/model_run2.py VCN_c$f -P initIV --configfile $CONFIG -D Full --dendritemode passive --datatable data_XM13A_nacncoop_pasdend
    python vcnmodel/model_run2.py VCN_c$f -P runANPSTH  -r 25 --configfile $CONFIG -D Full --dendritemode passive --datatable data_XM13A_nacncoop_pasdend
done
wait

for f in $CELLNAMES
do
    echo $f
    python vcnmodel/model_run2.py VCN_c$f -P initIV --configfile $CONFIG -D Full --dendritemode active --datatable data_XM13A_nacncoop_actdend
    python vcnmodel/model_run2.py VCN_c$f -P runANPSTH -r 25 --configfile $CONFIG -D Full --dendritemode active --datatable data_XM13A_nacncoop_actdend
done
wait
