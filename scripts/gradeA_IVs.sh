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
CELLNAMES="02 05 06 09 10 11 13 17 24 29 30"
CONFIG="autoscale.toml"
echo "running the individual initialization and/or running of IV protocols"
for f in $CELLNAMES
do
    echo $f
    #python vcnmodel/model_run2.py VCN_c$f -F -P initIV --configfile $CONFIG
    python vcnmodel/model_run2.py VCN_c$f -F -P runIV  --configfile $CONFIG &
done

wait

echo IV runs complete
# with "A", we do all cells in grade A
#python vcnmodel/plotters/plot_gbc_ivresults_2.py "A" -s