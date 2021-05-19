# Example:
# scripts/process_gbcIV.sh run all 
#
proto="runIV"
#######################################################
# Full models are from data/reconstuctions Matthew Kersting sent on
# March 6, 2020. 
# Note we do not have a full reconstruction for cell 18
# in that dataset.
#######################################################
CELLNAMES="18" # "02 05 06 09 10 11 13 17 30"
CONFIG="autoscale_xm13a_multisite_parallel.toml"
echo "running the individual initialization and/or running of IV protocols"
for f in $CELLNAMES
do
    echo $f
    python src/vcnmodel/model_run2.py VCN_c$f -P initIV --configfile $CONFIG -D Full --dendritemode normal --datatable data_XM13A_nacncoop
    python src/vcnmodel/model_run2.py VCN_c$f -P runIV  --configfile $CONFIG -D Full --dendritemode normal --datatable data_XM13A_nacncoop
done

wait

for f in $CELLNAMES
do
    echo $f
    python src/vcnmodel/model_run2.py VCN_c$f -P initIV --configfile $CONFIG -D Full --dendritemode passive --datatable data_XM13A_nacncoop_pasdend
    python src/vcnmodel/model_run2.py VCN_c$f -P runIV  --configfile $CONFIG -D Full --dendritemode passive --datatable data_XM13A_nacncoop_pasdend
done

for f in $CELLNAMES
do
    echo $f
    python src/vcnmodel/model_run2.py VCN_c$f -P initIV --configfile $CONFIG -D Full --dendritemode active --datatable data_XM13A_nacncoop_actdend
    python src/vcnmodel/model_run2.py VCN_c$f -P runIV  --configfile $CONFIG -D Full --dendritemode active --datatable data_XM13A_nacncoop_actdend
done

echo IV runs complete
# with "A", we do all cells in grade A
#python vcnmodel/plotters/plot_gbc_ivresults_2.py "A" -s