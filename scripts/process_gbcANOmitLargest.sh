#!/bin/bash

FILES="VCN_c09" # " VCN_c11 VCN_c17 VCN_c18"

for f in $FILES
do
    python vcnmodel/model_run.py $f -H -P initAN --configfile autoscale.toml 
    python vcnmodel/model_run.py $f -H -P runANPSTH -r 1 --Spirou removelargest --configfile autoscale.toml &
done

echo ANPSTH generators complete
# for f in $FILES
# do
#     echo "Cell: <$f>"
#     ls -lat VCN_Cells/$f/Simulations/AN
#     echo " "
# done


