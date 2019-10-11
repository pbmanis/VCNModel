#!/bin/bash

#python model_run.py VCN_c08 --protocol initAN --model XM13 -r 1 --sgcmodel cochlea -S MS -a 1.0 --noparallel
#python model_run.py VCN_c09 --protocol initAN --model mGBC -r 1 --sgcmodel cochlea -S MS -a 1.5 --noparallel
#python model_run.py VCN_c09h --protocol initAN --model mGBC -r 1 --sgcmodel cochlea -S MS -a 1.5
#python model_run.py VCN_c09nd --protocol initAN --model mGBC -r 1 --sgcmodel cochlea -S MS -a 1.5 --noparallel
#python model_run.py VCN_c17 --protocol initAN --model mGBC -r 1 --sgcmodel cochlea -S MS -a 3.0 --noparallel
#python model_run.py VCN_c18 --protocol initAN --model mGBC -r 1 --sgcmodel cochlea -S MS -a 3.0 --noparallel
#python model_run.py VCN_c19 --protocol initAN --model mGBC -r 1 --sgcmodel cochlea -S MS -a 1.5 --noparallel
#python model_run.py VCN_c20 --protocol initAN --model mGBC -r 1 --sgcmodel cochlea -S MS -a 1.5 --noparallel
#python model_run.py VCN_c21 --protocol initAN --model mGBC -r 1 --sgcmodel cochlea -S MS -a 3.0 --noparallel
#python model_run.py VCN_c22 --protocol initAN --model mGBC -r 1 --sgcmodel cochlea -S MS -a 3.0 --noparallel


# python model_run.py VCN_c08 --protocol runANPSTH --model XM13 -r 50 --sgcmodel cochlea -S MS -a 1.0 --noparallel &
python src/model_run.py VCN_c09 --protocol runANPSTH --model mGBC -r 50 -d 40 -f 4000. --sgcmodel cochlea -S MS -a 1.5 --noparallel &
#python model_run.py VCN_c09h --protocol runANPSTH --model mGBC -r 50 --sgcmodel cochlea -S MS -a 1.5 &
#python model_run.py VCN_c09nd --protocol runANPSTH --model mGBC -r 50 --sgcmodel cochlea -S MS -a 1.5 --noparallel &
# python model_run.py VCN_c17 --protocol runANPSTH --model mGBC -r 50 --sgcmodel cochlea -S MS -a 3.0 --noparallel &
# python model_run.py VCN_c18 --protocol runANPSTH --model mGBC -r 50 --sgcmodel cochlea -S MS -a 3.0 --noparallel &
# python model_run.py VCN_c19 --protocol runANPSTH --model mGBC -r 50 --sgcmodel cochlea -S MS -a 1.5 --noparallel &
# python model_run.py VCN_c20 --protocol runANPSTH --model mGBC -r 50 --sgcmodel cochlea -S MS -a 1.5 --noparallel &
# python model_run.py VCN_c21 --protocol runANPSTH --model mGBC -r 50 --sgcmodel cochlea -S MS -a 3.0 --noparallel &
# python model_run.py VCN_c22 --protocol runANPSTH --model mGBC -r 50 --sgcmodel cochlea -S MS -a 3.0 --noparallel &

wait
echo ANPSTH generators complete
FILES="VCN_c08 VCN_c09 VCN_c17 VCN_c18 VCN_c19 VCN_c20 VCN_c21 VCN_c22"
for f in $FILES
do
	echo "Cell: <$f>"
    ls -lat VCN_Cells/$f/Simulations/AN/AN_result*.p
    echo " "
done


