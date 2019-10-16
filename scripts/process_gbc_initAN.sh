#!/bin/bash

# python model_run.py VCN_c08 -M ${model} -T Bushy --sgcmodel cochlea -S MS -r 20 -P initAN &
# python model_run.py VCN_c09 -M ${model} -T Bushy --sgcmodel cochlea -S MS -r 20 -P initAN &
# python model_run.py VCN_c09h -M ${model} -T Bushy --sgcmodel cochlea -S MS -r 20 -P initAN &
#python model_run.py VCN_c09nd -M ${model} -T Bushy --sgcmodel cochlea -S MS -r 20 -P initAN &
# python model_run.py VCN_c17  -M ${model} -T Bushy --sgcmodel cochlea -S MS -r 20 -P initAN &
# python model_run.py VCN_c18  -M ${model} -T Bushy --sgcmodel cochlea -S MS -r 20 -P initAN &
# python model_run.py VCN_c19  -M ${model} -T Bushy --sgcmodel cochlea -S MS -r 20 -P initAN &
# python model_run.py VCN_c20 -M ${model} -T Bushy --sgcmodel cochlea -S MS -r 20 -P initAN &
# python model_run.py VCN_c21  -M ${model} -T Bushy --sgcmodel cochlea -S MS -r 20 -P initAN &
# python model_run.py VCN_c22  -M ${model} -T Bushy --sgcmodel cochlea -S MS -r 20 -P initAN &
FILES="VCN_c09 VCN_c11 VCN_c17 VCN_c18"
for f in $FILES
do
	echo "Cell: <$f>"
    python vcnmodel/model_run.py $f --protocol initAN -r 50 -d 40 -f 4000. --sgcmodel cochlea -S MS --configfile autoscale.toml &
done
wait
echo "AN initiation complete"

# for f in $FILES
# do
#     echo "Cell: <$f>"
#     ls -lat VCN_Cells/$f/Simulations/AN/AN_result*.p
#     echo " "
# done
