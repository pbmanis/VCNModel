

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
# python src/model_run.py VCN_c09 --protocol runANPSTH --model XM13_nacncoop -r 50 -d 40 -f 4000. --sgcmodel cochlea -S MS -a 1.5 --configfile autoscale.toml &
#python model_run.py VCN_c09h --protocol runANPSTH --model mGBC -r 50 --sgcmodel cochlea -S MS -a 1.5 &
#python model_run.py VCN_c09nd --protocol runANPSTH --model mGBC -r 50 --sgcmodel cochlea -S MS -a 1.5 --noparallel &
# python model_run.py VCN_c17 --protocol runANPSTH --model mGBC -r 50 --sgcmodel cochlea -S MS -a 3.0 --noparallel &
# python model_run.py VCN_c18 --protocol runANPSTH --model mGBC -r 50 --sgcmodel cochlea -S MS -a 3.0 --noparallel &
# python model_run.py VCN_c19 --protocol runANPSTH --model mGBC -r 50 --sgcmodel cochlea -S MS -a 1.5 --noparallel &
# python model_run.py VCN_c20 --protocol runANPSTH --model mGBC -r 50 --sgcmodel cochlea -S MS -a 1.5 --noparallel &
# python model_run.py VCN_c21 --protocol runANPSTH --model mGBC -r 50 --sgcmodel cochlea -S MS -a 3.0 --noparallel &
# python model_run.py VCN_c22 --protocol runANPSTH --model mGBC -r 50 --sgcmodel cochlea -S MS -a 3.0 --noparallel &
FILES="VCN_c09 VCN_c11 VCN_c17 VCN_c18"
((seed = 100))
((reps = 50))
for f in $FILES
do
    echo "Cell: $f"
    echo "seed: $seed"
    python vcnmodel/model_run.py $f -H --protocol initAN --seed $seed -r $reps -d 30 -f 16000. --sgcmodel cochlea -S MS --configfile autoscale.toml &
    python vcnmodel/model_run.py $f -H --protocol runANPSTH --seed $seed -r $reps -d 30 -f 16000. --sgcmodel cochlea -S MS --configfile autoscale.toml &
    ((seed = seed + 1))
done
wait
echo ANPSTH generators complete for $f

# for f in $FILES
# do
#     echo "Cell: <$f>"
#     ls -lat VCN_Cells/$f/Simulations/AN/AN_Result*
#     echo " "
# done


