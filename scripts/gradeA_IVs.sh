# Example:
# scripts/process_gbcIV.sh run all 
#
proto="runIV"
#######################################################
# Full models are from data Matthew Kersting sent on
# March 6, 2020. 
# Note we do not have a full reconstruction for cell 18
# in that dataset.
#######################################################
CELLNAMES="VCN_c02 VCN_c05 VCN_c06 VCN_c09 VCN_c10 VCN_c11 VCN_c13 VCN_c17 VCN_c24 VCN_c30 VCN_c31"
CELLNAMES="VCN_c02"
CELLNO="02 05 06 09 10 11 13 17 24 30 31"
echo "running the individual initialization and/or running of IV protocols"
for f in $CELLNAMES
do
    echo $f
    python vcnmodel/model_run.py $f -F -P initIV --configfile autoscale.toml
    python vcnmodel/model_run.py $f -F -P runIV  --configfile autoscale.toml
done

wait

echo IV runs complete

python vcnmodel/plotters/plot_gbc_ivresults_2.py $CELLNO