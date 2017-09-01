#!/bin/bash
FILES="VCN_c08 VCN_c09 VCN_c11 VCN_c14 VCN_c17 VCN_c18 VCN_c19 VCN_c20 VCN_c21 VCN_c22"
proto="runIV"
case $1 in
 init)
    proto="initIV"
    ;;
 justplot)
    proto="justplot"
    ;;
esac

model="XM13"
echo $proto

if ! [[ "$proto" = "justplot" ]] ; then
    #python model_run.py VCN_c08 -P  --model XM13 -r 10 --sgcmodel cochlea -S MS -a 1.0 --noparallel
    echo "what are we doing here?"
    python model_run.py VCN_c09nd -P ${proto} -M ${model} -r 10 --sgcmodel cochlea -S MS -a 1.5 --noparallel &
    python model_run.py VCN_c09 -P ${proto} -M ${model} -r 10 --sgcmodel cochlea -S MS -a 1.5 --noparallel &
    python model_run.py VCN_c11 -P ${proto} -M ${model} -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel &
    python model_run.py VCN_c14 -P ${proto} -M ${model} -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel &
    python model_run.py VCN_c17 -P ${proto} -M ${model} -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel &
    python model_run.py VCN_c18 -P ${proto} -M ${model} -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel &
    python model_run.py VCN_c19 -P ${proto} -M ${model} -r 10 --sgcmodel cochlea -S MS -a 1.5 --noparallel &
    python model_run.py VCN_c20 -P ${proto} -M ${model} -r 10 --sgcmodel cochlea -S MS -a 1.5 --noparallel &
    python model_run.py VCN_c21 -P ${proto} -M ${model} -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel &
    python model_run.py VCN_c22 -P ${proto} -M ${model} -r 10 --sgcmodel cochlea -S MS -a 3.0 --noparallel &

    wait
    echo ANPSTH generators complete

    for f in $FILES
    do
    	echo "Cell: <$f>"
        ls -lat VCN_Cells/$f/Simulations/IV
        echo " "
    done
fi

if [[ "$proto" = "runIV" ]] || [[ "$proto" = "justplot" ]] ; then
    python all_gbc_iv.py ${model}
fi
