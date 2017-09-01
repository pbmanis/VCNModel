#!/bin/bash

#python model_run.py VCN_c08 -P initAN --model XM13 -r 1 --sgcmodel cochlea -S MS -a 1.0 --noparallel
#python model_run.py VCN_c09 -P initAN -M mGBC -r 1 --sgcmodel cochlea -S MS -a 1.5 --noparallel
#python model_run.py VCN_c09h -P initAN -M mGBC -r 1 --sgcmodel cochlea -S MS -a 1.5
#python model_run.py VCN_c09nd -P initAN -M mGBC -r 1 --sgcmodel cochlea -S MS -a 1.5 --noparallel
#python model_run.py VCN_c17 -P initAN -M mGBC -r 1 --sgcmodel cochlea -S MS -a 3.0 --noparallel
#python model_run.py VCN_c18 -P initAN -M mGBC -r 1 --sgcmodel cochlea -S MS -a 3.0 --noparallel
#python model_run.py VCN_c19 -P initAN -M mGBC -r 1 --sgcmodel cochlea -S MS -a 1.5 --noparallel
#python model_run.py VCN_c20 -P initAN -M mGBC -r 1 --sgcmodel cochlea -S MS -a 1.5 --noparallel
#python model_run.py VCN_c21 -P initAN -M mGBC -r 1 --sgcmodel cochlea -S MS -a 3.0 --noparallel
#python model_run.py VCN_c22 -P initAN -M mGBC -r 1 --sgcmodel cochlea -S MS -a 3.0 --noparallel

modstuffs="-r 20 --fmod 50. --depth 100. --stimulus SAM"

python model_run.py VCN_c08 -T Bushy -M mGBC --sgcmodel cochlea --SR MS ${modstuffs}  -P runANPSTH --noparallel &
python model_run.py VCN_c09 -T Bushy -M mGBC --sgcmodel cochlea --SR MS ${modstuffs}  -P runANPSTH -a 1.0 --noparallel &
python model_run.py VCN_c09h -T Bushy -M mGBC --sgcmodel cochlea --SR MS ${modstuffs}  -P runANPSTH -a 1.5 --noparallel  &
python model_run.py VCN_c09nd -T Bushy -M mGBC --sgcmodel cochlea --SR MS ${modstuffs}  -P runANPSTH -a 1.5 --noparallel  &
python model_run.py VCN_c17 -T Bushy -M mGBC --sgcmodel cochlea --SR MS ${modstuffs}  -P runANPSTH -a 1.5 --noparallel  &
python model_run.py VCN_c18 -T Bushy -M mGBC --sgcmodel cochlea --SR MS ${modstuffs}  -P runANPSTH -a 3.0 --noparallel  &
python model_run.py VCN_c19 -T Bushy -M mGBC --sgcmodel cochlea --SR MS ${modstuffs}  -P runANPSTH -a 3.0 --noparallel  &
python model_run.py VCN_c20 -T Bushy -M mGBC --sgcmodel cochlea --SR MS ${modstuffs}  -P runANPSTH -a 1.5 --noparallel  &
python model_run.py VCN_c21 -T Bushy -M mGBC --sgcmodel cochlea --SR MS ${modstuffs}  -P runANPSTH -a 3.0 --noparallel  &
python model_run.py VCN_c22 -T Bushy -M mGBC --sgcmodel cochlea --SR MS ${modstuffs}  -P runANPSTH -a 3.0 --noparallel

wait
echo ANPSTH generators complete
FILES="VCN_c08 VCN_c09 VCN_c17 VCN_c18 VCN_c19 VCN_c20 VCN_c21 VCN_c22"
for f in $FILES
do
	echo "Cell: <$f>"
    ls -lat VCN_Cells/$f/Simulations/AN/AN_result*.p
    echo " "
done


