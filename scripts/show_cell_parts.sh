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
           python vcnmodel/model_run2.py VCN_c$f -F --configfile autoscale_multisite_parallel.toml  --dendritemode $D --displaymode mechanism --mechanism $C --style cylinders --displayscale
           # hocRender ../VCN-SBEM-Data/VCN_Cells/VCN_c$f/Morphology/VCN_c$f\_Full.hoc -m cylinders -r pyqtgraph -d sec-type -M $C
       done
    done
done

wait