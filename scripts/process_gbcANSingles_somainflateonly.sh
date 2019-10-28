#!/bin/bash
mod_params="-r 1 --stimulus tonepip"
stim_params="-f 16000. --dB 30."
AN_params=" --sgcmodel cochlea --SR MS"

FILES="VCN_c09 VCN_c11 VCN_c17 VCN_c18"
for f in $FILES
do
    echo "Cell: $f, self inputs"
    python vcnmodel/model_run.py $f --protocol initAN -H \
         ${mod_params} ${stim_params} ${AN_params} --configfile autoscale_somaonly.toml\
             --saveall

    python vcnmodel/model_run.py $f --protocol runANSingles -H \
         ${mod_params} ${stim_params} ${AN_params} --configfile autoscale_somaonly.toml\
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


