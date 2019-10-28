#!/bin/bash

mod_params="-r 50 --fmod 200. --dmod 100 --stimulus SAM"
stim_params="-f 16000. --dB 30"
AN_params=" --sgcmodel cochlea --SR MS"
# cellstuffs="--type Bushy --modeltype II -M XM13nacncoop -H"
# endbulbs="--inputpattern VCN_c10"
# endbulbs_self=" "
# inflation="--soma-autoinflate --dendrite-autoinflate"


# python model_run.py VCN_c09 ${cellstuffs} ${AN_params} ${stim_params} ${endbulbs} ${mod_params} ${inflation} -P runANPSTH --noparallel &
# python model_run.py VCN_c09 ${cellstuffs} ${AN_params} ${stim_params} ${endbulbs_self} ${mod_params} ${inflation} -P runANPSTH --noparallel &

((seed = 100))

FILES="VCN_c09 VCN_c11 VCN_c17 VCN_c18"
for f in $FILES
do
    echo "Cell: $f"
    python vcnmodel/model_run.py $f -H --protocol initAN --seed $seed -${mod_params} ${stim_params} ${AN_params} --configfile autoscale.toml &
    python vcnmodel/model_run.py $f -H --protocol runANPSTH --seed $seed -${mod_params} ${stim_params} ${AN_params} --configfile autoscale.toml &
    ((seed = seed + 12))
done
wait
echo ANPSTH generators complete

# for f in $FILES
# do
#     echo "Cell: <$f>"
#     ls -lat VCN_Cells/$f/Simulations/AN/AN_Result*
#     echo " "
# done


