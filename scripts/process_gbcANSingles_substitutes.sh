#!/bin/bash

mod_params="-r 1 --stimulus tonepip"
stim_params="-f 16000. --dB 30."
AN_params=" --sgcmodel cochlea --SR MS"

targetcell="VCN_c09"
FILES="VCN_c02 VCN_c05 VCN_c10 VCN_c13"
for f in $FILES
do
    echo "Cell: $f"
    python vcnmodel/model_run.py $targetcell -i $f --protocol initAN -H \
          ${mod_params} ${stim_params} ${AN_params} --configfile autoscale.toml\
             --saveall

    python vcnmodel/model_run.py $targetcell -i $f --protocol runANSingles -H \
          ${mod_params} ${stim_params} ${AN_params} --configfile autoscale.toml\
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


