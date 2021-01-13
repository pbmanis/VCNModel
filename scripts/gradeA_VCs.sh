# Example:
# scripts/process_gbcIV.sh run all 
#
proto="runVC"
#######################################################
# Full models are from data/reconstuctions Matthew Kersting sent on
# March 6, 2020. 
# Note we do not have a full reconstruction for cell 18
# in that dataset.
#######################################################
CELLNAMES="17" # "02 05 06 09 10 11 13 17 24 29 30"
CONFIG="autoscale_xm13a_multisite_parallel.toml"
echo "running the individual initialization and/or running of VC protocols"
for f in $CELLNAMES
do
    echo $f
    python vcnmodel/model_run2.py VCN_c$f -D Full -P initVC --configfile $CONFIG --dendritemode normal --datatable data_XM13A_nacncoop
    python vcnmodel/model_run2.py VCN_c$f -D Full -P runVC  --configfile $CONFIG --dendritemode normal --datatable data_XM13A_nacncoop
done

wait

for f in $CELLNAMES
do
    echo $f
    python vcnmodel/model_run2.py VCN_c$f -D Full -P initVC --configfile $CONFIG --dendritemode passive --datatable data_XM13A_nacncoop_pasdend
    python vcnmodel/model_run2.py VCN_c$f -D Full -P runVC  --configfile $CONFIG --dendritemode passive --datatable data_XM13A_nacncoop_pasdend
done

wait

for f in $CELLNAMES
do
    echo $f
    python vcnmodel/model_run2.py VCN_c$f -D Full -P initVC --configfile $CONFIG --dendritemode active --datatable data_XM13A_nacncoop_actdend
    python vcnmodel/model_run2.py VCN_c$f -D Full -P runVC  --configfile $CONFIG --dendritemode active --datatable data_XM13A_nacncoop_actdend
done

echo VC runs complete
# with "A", we do all cells in grade A
#python vcnmodel/plotters/plot_gbc_ivresults_2.py "A" -s