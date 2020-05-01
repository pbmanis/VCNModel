#!/bin/bash

# do only cells with morphology
# ONLY "A" grade cells with complete dendrites:
#    gradeA = [2, 5, 6, 9, 10, 11, 13, 17, 24, 29, 30]
FILES="02 05 06 09 10 11 13 17 24 29 30"
proto="initIV"

echo "running the individual initialization or IV protocols"
for f in $FILES
do
    echo $f
    python vcnmodel/model_run2.py VCN_c${f} -F -P ${proto} --configfile autoscale.toml &
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
    ls -lat ../VCN-SBEM-Data/VCN_Cells/VCN_c${f}/Initialization/IV*.dat
    echo " "
done
