# Example:
# scripts/process_gbcIV.sh run all 
#

#######################################################
# Full models are from data/reconstuctions Matthew Kersting sent on
# March 6, 2020. 
# Note we do not have a full reconstruction for cell 18
# in that dataset.
#######################################################
CELLNAMES="13" # " 13 24" # "02" #" 05 06 09 10 11 13 17 24 29 30"
CONFIG="noscale.toml" #"autoscale.toml"
#CONFIG="autoscale_simple.toml"
RUNTEXT="running the individual initialization and running AN PSTH protocols"
echo $RUNTEXT
for f in $CELLNAMES
do
    echo $f
    python vcnmodel/model_run2.py VCN_c$f  -F -P initAN --configfile $CONFIG
    python vcnmodel/model_run2.py VCN_c$f  -F -P runANPSTH -r 10 --configfile $CONFIG
done

wait

echo AN runs complete
# with "A", we do all cells in grade A
# python vcnmodel/plotters/plot_gbc_ivresults_2.py "A" -s