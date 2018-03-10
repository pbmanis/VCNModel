#!/bin/bash
FILES="VCN_c08 VCN_c09 VCN_c09nd VCN_c11 VCN_c14 VCN_c17 VCN_c18 VCN_c19 VCN_c20 VCN_c21 VCN_c22"
#python model_run.py VCN_c08 --protocol runIV --model XM13 -r 10 --sgcmodel cochlea -S MS -a 1.0 --noparallel
python model_run.py VCN_c09nd --hoc VCN_c09nd.hoc --protocol initIV --model XM13 --noparallel &
python model_run.py VCN_c09 --hoc VCN_c09.hoc --protocol initIV --model XM13 --noparallel &
python model_run.py VCN_c11 --hoc VCN_c11.hoc --protocol initIV --model XM13 --noparallel &
python model_run.py VCN_c14 --hoc VCN_c14.hoc --protocol initIV --model XM13 --noparallel &
python model_run.py VCN_c17 --hoc VCN_c17.hoc --protocol initIV --model XM13 --noparallel &
python model_run.py VCN_c18 --hoc VCN_c18.hoc --protocol initIV --model XM13 --noparallel &
python model_run.py VCN_c19 --hoc VCN_c19.hoc --protocol initIV --model XM13 --noparallel &
python model_run.py VCN_c20 --hoc VCN_c20.hoc --protocol initIV --model XM13 --noparallel &
python model_run.py VCN_c21 --hoc VCN_c21.hoc --protocol initIV --model XM13 --noparallel &
python model_run.py VCN_c22 --hoc VCN_c22.hoc --protocol initIV --model XM13 --noparallel &

wait
echo initIV for all cells complete

for f in $FILES
do
	echo "Cell: <$f>"
    ls -lat VCN_Cells/$f/Initialization
    echo " "
done