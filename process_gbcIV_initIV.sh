#!/bin/bash
FILES="VCN_c08 VCN_c09 VCN_c11 VCN_c14 VCN_c17 VCN_c18 VCN_c19 VCN_c20 VCN_c21 VCN_c22"
#python model_run.py VCN_c08 --protocol runIV --model XM13 -r 10 --sgcmodel cochlea -S MS -a 1.0 --noparallel
python model_run.py VCN_c09 --protocol initIV --model mGBC -r 10 --sgcmodel cochlea -S MS -a 1.5 --noparallel
python model_run.py VCN_c17 --protocol initIV --model mGBC -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel
python model_run.py VCN_c11 --protocol initIV --model mGBC -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel
python model_run.py VCN_c14 --protocol initIV --model mGBC -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel
python model_run.py VCN_c18 --protocol initIV --model mGBC -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel
python model_run.py VCN_c19 --protocol initIV --model mGBC -r 10 --sgcmodel cochlea -S MS -a 1.5 --noparallel
python model_run.py VCN_c20 --protocol initIV --model mGBC -r 10 --sgcmodel cochlea -S MS -a 1.5 --noparallel
python model_run.py VCN_c21 --protocol initIV --model mGBC -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel
python model_run.py VCN_c22 --protocol initIV --model mGBC -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel

wait
echo ANPSTH generators complete

for f in $FILES
do
	echo "Cell: <$f>"
    ls -lat VCN_Cells/$f/Simulations/IV
    echo " "
done