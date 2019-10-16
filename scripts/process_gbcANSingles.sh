#!/bin/bash


FILES="VCN_c09 VCN_c11 VCN_c17 VCN_c18"
for f in $FILES
do
    echo "Cell: $f"
    python vcnmodel/model_run.py $f --protocol initAN -H \
         -r 1 -d 40 -f 4000. --sgcmodel cochlea -S MS --configfile autoscale.toml\
             --saveall

    python vcnmodel/model_run.py $f --protocol runANSingles -H \
         -r 1 -d 40 -f 4000. --sgcmodel cochlea -S MS --configfile autoscale.toml\
             --saveall
done
wait
echo ANPSTH generators complete
# for f in $FILES
# do
#     echo "Cell: <$f>"
#     ls -lat VCN_Cells/$f/Simulations/AN
#     echo " "
# done


