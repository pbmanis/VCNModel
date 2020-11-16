# Example:
# scripts/process_gbcIV.sh run all 
#
proto="testIV"
#######################################################
# Full models are from data/reconstuctions Matthew Kersting sent on
# March 6, 2020. 
# Note we do not have a full reconstruction for cell 18
# in that dataset.
#######################################################
CELLNAMES="18" # "02 05 06 09 10 11 13 17 30"
CONFIG="toml/autoscale_multisite_parallel.toml"
# CONFIG="autoscale_xm13a_multisite_parallel.toml"
echo "computing Zin for each gradeA Cell"


# for f in $CELLNAMES
# do
#     echo $f
#     python vcnmodel/model_run2.py VCN_c$f -P initIV --configfile $CONFIG -D Full  --dendritemode normal --datatable data_XM13A_nacncoop
#     python vcnmodel/model_run2.py VCN_c$f -P Zin --configfile $CONFIG -D Full  --dendritemode normal --datatable data_XM13A_nacncoop
# done
# wait
#
# for f in $CELLNAMES
# do
#     echo $f
#     python vcnmodel/model_run2.py VCN_c$f -P initIV --configfile $CONFIG -D Full --dendritemode passive --datatable data_XM13A_nacncoop_pasdend
#     python vcnmodel/model_run2.py VCN_c$f -P Zin --configfile $CONFIG -D Full  --dendritemode passive --datatable data_XM13A_nacncoop_pasdend
# done
# wait
# for f in $CELLNAMES
# do
#     echo $f
#     python vcnmodel/model_run2.py VCN_c$f -P initIV --configfile $CONFIG -D Full  --dendritemode active --datatable data_XM13A_nacncoop_actdend
#     python vcnmodel/model_run2.py VCN_c$f -P Zin --configfile $CONFIG -D Full  --dendritemode active --datatable data_XM13A_nacncoop_actdend
# done
# wait

for f in $CELLNAMES
do
    echo $f
    python vcnmodel/model_run2.py VCN_c$f -P initIV --configfile $CONFIG -D NoDend  --dendritemode normal --datatable data_XM13A_nacncoop
    python vcnmodel/model_run2.py VCN_c$f -P Zin --configfile $CONFIG -D NoDend  --dendritemode normal --datatable data_XM13A_nacncoop
done
wait

echo Zin runs complete