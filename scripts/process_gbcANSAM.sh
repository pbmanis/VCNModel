#!/bin/bash

modstuffs="-r 50 --fmod 200. --dmod 100. --stimulus SAM"
stimstuff="-f 16000. --dB 30"
ANstuff=" --sgcmodel cochlea --SR MS"
cellstuffs="--type Bushy --modeltype II -M XM13nacncoop -H"
endbulbs="--inputpattern VCN_c10"
endbulbs_self=" "
inflation="--soma-autoinflate --dendrite-autoinflate"

python model_run.py VCN_c09 ${cellstuffs} ${ANstuff} ${stimstuff} ${endbulbs} ${modstuffs} ${inflation} -P runANPSTH --noparallel &
python model_run.py VCN_c09 ${cellstuffs} ${ANstuff} ${stimstuff} ${endbulbs_self} ${modstuffs} ${inflation} -P runANPSTH --noparallel &


wait
echo ANPSTH generators complete
FILES="VCN_c08 VCN_c09 VCN_c17 VCN_c18 VCN_c19 VCN_c20 VCN_c21 VCN_c22"
for f in $FILES
do
	echo "Cell: <$f>"
    ls -lat VCN_Cells/$f/Simulations/AN/AN_result*.p
    echo " "
done


