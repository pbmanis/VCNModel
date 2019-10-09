#!/bin/bash
#FILES="VCN_c08 VCN_c09 VCN_c09nd VCN_c11 VCN_c14 VCN_c16 VCN_c17 VCN_c18 VCN_c19 VCN_c20 VCN_c21 VCN_c22"
model="XM13_nabu"
proto="runIV"
case $1 in
 run)
    proto="runIV"
    ;;
 init)
    proto="initIV"
    ;;
 plot)
    proto="plot"
    ;;
 list)
    proto="list"
    ;;
esac

case $2 in
 all)
    # FILES="VCN_c08 VCN_c09 VCN_c11 VCN_c14 VCN_c16 VCN_c17 VCN_c18 VCN_c19 VCN_c20 VCN_c21 VCN_c22"
    FILES="VCN_c18 VCN_c19 VCN_c20 VCN_c21 VCN_c22"
    ;;
 *)
    FILES=$2
    ;;
esac

case $3 in
  plot)
    PLOT="--plot"
    ;;
  *)
    PLOT=""
    ;;
esac

echo $proto

if [[ "$proto" = "initIV" ]] || [[ "$proto" = "runIV" ]] ; then
    #MS -a 1.0 --noparallel
    echo "running the individual initialization and/or running of IV protocols"
    for f in $FILES
    do
        echo $f
        python src/model_run.py ${f} -H -P ${proto} -M ${model} ${PLOT} --soma-autoinflate --noparallel &
    # python model_run.py VCN_c09 -P ${proto} -M ${model} --noparallel &
    # python model_run.py VCN_c11 -P ${proto} -M ${model} --noparallel &
    # python model_run.py VCN_c14 -P ${proto} -M ${model} --noparallel &
    # python model_run.py VCN_c16 -P ${proto} -M ${model} --noparallel &
    # python model_run.py VCN_c17 -P ${proto} -M ${model} --noparallel &
    # python model_run.py VCN_c18 -P ${proto} -M ${model} --noparallel &
    # python model_run.py VCN_c19 -P ${proto} -M ${model} --noparallel &
    # python model_run.py VCN_c20 -P ${proto} -M ${model} --noparallel &
    # python model_run.py VCN_c21 -P ${proto} -M ${model} --noparallel &
    # python model_run.py VCN_c22 -P ${proto} -M ${model} --noparallel &
    done

    wait
    echo IV runs complete

    for f in $FILES
    do
    	echo "Cell: <$f>"
        echo "VCN_Cells/$f/Simulations/IV"
        ls -lat VCN_Cells/$f/Simulations/IV
        echo " "
    done
fi

if [[ "$proto" = "runIV" ]] || [[ "$proto" = "plot" ]] ; then
    python src/all_gbc_iv.py ${model}
fi
#
if [[ "$proto" = 'list' ]] ; then
    for f in $FILES
    do
        echo "Cell: <$f>"
        ls -lat VCN_Cells/$f/Simulations/IV
        echo " "
    done
fi