#!/bin/bash

# do only cells with morphology
# ONLY "A" grade cells with complete dendrites:
FILES="VCN_c09 VCN_c11 VCN_c17 VCN_c18"
proto="initIV"

echo "running the individual initialization or IV protocols"
for f in $FILES
do
    echo $f
    python src/model_run.py ${f} -H -P ${proto} --configfile autoscale.toml &
done
# python model_run.py VCN_c08 -H --protocol initIV --model XM13 --noparallel &
# python model_run.py VCN_c09 -H --protocol initIV --model mGBC --noparallel &
# python model_run.py VCN_c17 -H --protocol initIV --model mGBC --noparallel &
# python model_run.py VCN_c18 -H --protocol initIV --model mGBC --noparallel &
# python model_run.py VCN_c19 -H --protocol initIV --model mGBC --noparallel &
# python model_run.py VCN_c20 -H --protocol initIV --model mGBC --noparallel &
# python model_run.py VCN_c21 -H --protocol initIV --model mGBC --noparallel &
# python model_run.py VCN_c22 -H --protocol initIV --model mGBC --noparallel &

wait
echo "IV initialization complete"

for f in $FILES
do
	echo "Cell: <$f>"
    ls -lat VCN_Cells/$f/Initialization/IV*.dat
    echo " "
done
