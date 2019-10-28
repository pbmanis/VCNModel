# Example:
# scripts/process_gbcIV.sh run all 
#
proto="runIV"
FILES="VCN_c09 VCN_c11 VCN_c17 VCN_c18"
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
    FILES="VCN_c09 VCN_c11 VCN_c17 VCN_c18"
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
        python src/model_run.py $f -H -P initIV ${PLOT} --noparallel --configfile autoscale.toml &
        python src/model_run.py $f -H -P ${proto} ${PLOT} --noparallel --configfile autoscale.toml &
    done

    wait
    echo IV runs complete

    # for f in $FILES
    # do
    #     echo "Cell: <$f>"
    #     echo "VCN_Cells/$f/Simulations/IV"
    #     ls -lat VCN_Cells/$f/Simulations/IV
    #     echo " "
    # done
fi

echo $proto
# if [[ "$proto" = "runIV" ]] || [[ "$proto" = "plot" ]] ; then
#     python src/all_gbc_iv.py ${model}
# fi

if [[ "$proto" = 'list' ]]
    then
        echo $FILES
        for f in $FILES
        do
            echo "Cell: <$f>"
            ls -lat VCN_Cells/$f/Simulations/IV
            echo " "
        done
fi