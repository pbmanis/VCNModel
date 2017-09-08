#!/bin/bash

#python model_run.py VCN_c08 --protocol runANSingles --model XM13 -r 10 --sgcmodel cochlea -S MS -a 1.0 --noparallel &
# VCN_09 was -a 1.5
python model_run.py VCN_c09 --protocol runANSingles --model mGBC -r 10 --sgcmodel cochlea -S MS -a 1.0 --plot --noparallel &
# python model_run.py VCN_c17 --protocol runANSingles --model mGBC -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel &
# python model_run.py VCN_c18 --protocol runANSingles --model mGBC -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel &
# python model_run.py VCN_c19 --protocol runANSingles --model mGBC -r 10 --sgcmodel cochlea -S MS -a 1.5 --noparallel &
# python model_run.py VCN_c20 --protocol runANSingles --model mGBC -r 10 --sgcmodel cochlea -S MS -a 1.5 --noparallel &
# python model_run.py VCN_c21 --protocol runANSingles --model mGBC -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel &
# python model_run.py VCN_c22 --protocol runANSingles --model mGBC -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel &

wait
echo ANPSTH generators complete
FILES="VCN_c08 VCN_c09 VCN_c17 VCN_c18 VCN_c19 VCN_c20 VCN_c21 VCN_c22"
for f in $FILES
do
	echo "Cell: <$f>"
    ls -lat VCN_Cells/$f/Simulations/AN
    echo " "
done


