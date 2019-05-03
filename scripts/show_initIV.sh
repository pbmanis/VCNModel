#!/bin/bash
FILES="VCN_c08 VCN_c09 VCN_c09nd VCN_c11 VCN_c14 VCN_c17 VCN_c18 VCN_c19 VCN_c20 VCN_c21 VCN_c22"

for f in $FILES
do
	echo "Cell: <$f>"
    # rm VCN_Cells/$f/Initialization/IVneuronState.dat
    ls -lat VCN_Cells/$f/Initialization/IVneuronState_mGBC.dat
    echo " "
done