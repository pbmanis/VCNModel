#!/bin/bash

#python model_run.py VCN_c08 --protocol runANSingles --model ${model} -r 10 --sgcmodel cochlea -S MS -a 1.0 --noparallel &
# VCN_09 was -a 1.5
#python model_run.py VCN_c09 --protocol runANSingles -M ${model} -r 10 --sgcmodel cochlea -S MS -a 1.5 --plot 
#python model_run.py VCN_c09 --type Bushy --modeltype II --protocol runANSingles -H --model XM13nacncoop  -r 5 --sgcmodel cochlea -S MS  -f 16000. --soma-autoinflate --dendrite-autoinflate --noparallel
# python model_run.py VCN_c09 --type Bushy --modeltype II --protocol runANSingles -H --model XM13nacncoop --inputpattern VCN_c10 -r 5 --sgcmodel cochlea -S MS  -f 16000. --soma-autoinflate --dendrite-autoinflate --noparallel
# python model_run.py VCN_c17 --protocol runANSingles --model ${model} -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel &
# python model_run.py VCN_c18 --protocol runANSingles --model ${model} -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel &
# python model_run.py VCN_c19 --protocol runANSingles --model ${model} -r 10 --sgcmodel cochlea -S MS -a 1.5 --noparallel &
# python model_run.py VCN_c20 --protocol runANSingles --model ${model} -r 10 --sgcmodel cochlea -S MS -a 1.5 --noparallel &
# python model_run.py VCN_c21 --protocol runANSingles --model ${model} -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel &
# python model_run.py VCN_c22 --protocol runANSingles --model ${model} -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel &

FILES="VCN_c09 VCN_c11 VCN_c17 VCN_c18"
for f in $FILES
do
    echo "Cell: $f"

    python src/model_run.py $f --protocol runANSingles -H \
         -r 1 -d 40 -f 4000. --sgcmodel cochlea -S MS --configfile autoscale.toml\
             --saveall
done
wait
echo ANPSTH generators complete
# FILES="VCN_c08 VCN_c09 VCN_c09nd VCN_c17 VCN_c18 VCN_c19 VCN_c20 VCN_c21 VCN_c22"
# for f in $FILES
# do
#     echo "Cell: <$f>"
#     ls -lat VCN_Cells/$f/Simulations/AN
#     echo " "
# done


