#!/bin/bash

#CONFIG="noscale.toml" #"autoscale.toml"
CONFIG="--configfile singles_multisite_parallel.toml"
RUNTEXT="running the AN Single protocols Full"

# sequentials:
CELLNAMES="02 05 06 09 10 11 13 17 18 30"
TABLES="data_XM13A_nacncoop_normal" # " data_XM13A_nacncoop_pasdend data_XM13A_nacncoop_actdend"
DENDRITES="Full" # "NoDend"
NREPS="25"
INPUTSIZE="150"  # ASA, um2
echo $RUNTEXT

for t in $TABLES
    do
    for f in $CELLNAMES
        do
        for d in $DENDRITES
            do
            AXON=""
            case $f in
                02 | 05)
                    AXON="-A standardized"
                    ;;
            esac
            echo $f
            # python vcnmodel/model_run2.py VCN_c$f -P initAN -D $d $AXON $CONFIG --test_input $INPUTSIZE--datatable $t
            python vcnmodel/model_run2.py VCN_c$f -P runANSingles --test_input $INPUTSIZE -r $NREPS -D $d $AXON $CONFIG  --datatable $t
         done
    done
done

wait
echo AN Singles oneinput testing complete
# for f in $FILES
# do
#     echo "Cell: <$f>"
#     ls -lat VCN_Cells/$f/Simulations/AN
#     echo " "
# done

