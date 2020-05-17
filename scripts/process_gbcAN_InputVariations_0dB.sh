#!/bin/bash
CONFIG="autoscale_multisite_0dB_parallel.toml"
RUNTEXT="running the individual initialization and running AN PSTH protocols"
CELLNAMES="02 05 06 09 10 11 13 17 30"
REPS=50
for f in $CELLNAMES
do
    # python vcnmodel/model_run2.py VCN_c$f -H -P runANPSTH -r $REPS --configfile $CONFIG
    # python vcnmodel/model_run2.py VCN_c$f -H -P runANPSTH -r $REPS --Spirou all=mean --configfile $CONFIG
    # python vcnmodel/model_run2.py VCN_c$f -H -P runANPSTH -r $REPS --Spirou max=mean --configfile $CONFIG
    python vcnmodel/model_run2.py VCN_c$f -H -P runANPSTH -r $REPS --Spirou removelargest --configfile $CONFIG
    # python vcnmodel/model_run2.py VCN_c$f -H -P runANPSTH -r $REPS --Spirou largestonly --configfile $CONFIG
    # python vcnmodel/model_run2.py VCN_c$f -H -P runANPSTH -r $REPS --Spirou twolargest --configfile $CONFIG
done

wait

echo All Variations Complete



