#!/bin/bash
echo "this is old, use process_gbcIV.sh instead"
exit
FILES="VCN_c08 VCN_c09 VCN_c17 VCN_c18 VCN_c19 VCN_c20 VCN_c21 VCN_c22"
python model_run.py VCN_c08 -H --protocol initIV --model XM13 -r 10 --sgcmodel cochlea -S MS -a 1.0 --noparallel &
python model_run.py VCN_c09 -H --protocol initIV --model mGBC -r 10 --sgcmodel cochlea -S MS -a 1.5 --noparallel &
python model_run.py VCN_c17 -H --protocol initIV --model mGBC -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel &
python model_run.py VCN_c18 -H --protocol initIV --model mGBC -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel &
python model_run.py VCN_c19 -H --protocol initIV --model mGBC -r 10 --sgcmodel cochlea -S MS -a 1.5 --noparallel &
python model_run.py VCN_c20 -H --protocol initIV --model mGBC -r 10 --sgcmodel cochlea -S MS -a 1.5 --noparallel &
python model_run.py VCN_c21 -H --protocol initIV --model mGBC -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel &
python model_run.py VCN_c22 -H --protocol initIV --model mGBC -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel &

wait
echo IV initialization complete

for f in $FILES
do
	echo "Cell: <$f>"
    ls -lat VCN_Cells/$f/Initialization/IV*.dat
    echo " "
done