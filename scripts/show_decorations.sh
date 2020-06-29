CH="klt kht ihvcn nacncoop"
CELLNAMES="02 05 06 09 10 11 13 17 24 29 30"
DM="normal active passive"
for f in $CELLNAMES
do
    for C in $CH
    do
       for D in $DM
       do
           echo "Plotting VCN_c$f $C $D"
            python vcnmodel/model_run2.py VCN_c$f -F -v mechanism --configfile autoscale_multisite_parallel.toml --dendritemode $D -c $C
       done
    done
done

wait