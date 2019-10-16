#!/bin/bash

modstuffs="-r 50 --fmod 200. --dmod 100. --stimulus SAM"
stimstuff="-f 4000. --dB 40"
ANstuff=" --sgcmodel cochlea --SR MS"
# cellstuffs="--type Bushy --modeltype II -M XM13nacncoop -H"
# endbulbs="--inputpattern VCN_c10"
# endbulbs_self=" "
# inflation="--soma-autoinflate --dendrite-autoinflate"


# python model_run.py VCN_c09 ${cellstuffs} ${ANstuff} ${stimstuff} ${endbulbs} ${modstuffs} ${inflation} -P runANPSTH --noparallel &
# python model_run.py VCN_c09 ${cellstuffs} ${ANstuff} ${stimstuff} ${endbulbs_self} ${modstuffs} ${inflation} -P runANPSTH --noparallel &


FILES="VCN_c09 VCN_c11 VCN_c17 VCN_c18"
for f in $FILES
do
    echo "Cell: $f"
    python src/model_run.py $f --protocol runANPSTH -${modstuffs} ${stimstuff} ${ANstuff} --configfile autoscale.toml &
done
wait
echo ANPSTH generators complete

for f in $FILES
do
	echo "Cell: <$f>"
    ls -lat VCN_Cells/$f/Simulations/AN/AN_Result*
    echo " "
done


