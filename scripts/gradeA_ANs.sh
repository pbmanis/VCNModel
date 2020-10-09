# Example:
# scripts/process_gbcIV.sh run all 
#

#######################################################
# Full models are from data/reconstuctions Matthew Kersting sent on
# March 6, 2020. 
# Note we do not have a full reconstruction for cell 18
# in that dataset.
#######################################################
CELLNAMES="30" # "02 05 06 09 10 13 17"
#CONFIG="noscale.toml" #"autoscale.toml"
CONFIG="autoscale_xm13a_multisite_parallel.toml"
RUNTEXT="running the individual initialization and running AN PSTH protocols"
WORKERS="16"
REPS="100"
echo $RUNTEXT
for f in $CELLNAMES
do
    echo $f
    python vcnmodel/model_run2.py VCN_c$f  -F -P initAN --configfile $CONFIG  --datatable data_XM13A_nacncoop
    python vcnmodel/model_run2.py VCN_c$f  -F -P runANPSTH -r $REPS --dB 10 --Spirou all --configfile $CONFIG --datatable data_XM13A_nacncoop
    if [ $? -ne 0 ]; then
        exit 1
    fi
    python vcnmodel/model_run2.py VCN_c$f  -F -P runANPSTH -r $REPS --dB 10 --Spirou largestonly --configfile $CONFIG --datatable data_XM13A_nacncoop
    if [ $? -ne 0 ]; then
        exit 1
    fi
    python vcnmodel/model_run2.py VCN_c$f  -F -P runANPSTH -r $REPS --dB 10 --Spirou removelargest --configfile $CONFIG --datatable data_XM13A_nacncoop
   if [ $? -ne 0 ]; then
        exit 1
    fi
    # python vcnmodel/model_run2.py VCN_c$f  -F -P runANPSTH -r $REPS --dB 20 --Spirou largestonly --workers $WORKERS --configfile $CONFIG --datatable data_XM13A_nacncoop
    # if [ $? -ne 0 ]; then
    #     exit 1
    # fi
    # python vcnmodel/model_run2.py VCN_c$f  -F -P runANPSTH -r $REPS --dB 20 --Spirou removelargest --workers $WORKERS --configfile $CONFIG --datatable data_XM13A_nacncoop
    # if [ $? -ne 0 ]; then
    #     exit 1
    # fi
done

# wait
# echo $RUNTEXT
# for f in $CELLNAMES
# do
#     echo $f
#     python vcnmodel/model_run2.py VCN_c$f  -F -P initAN --configfile $CONFIG  --datatable data_XM13A_nacncoop_pasdend
#     python vcnmodel/model_run2.py VCN_c$f  -F -P runANPSTH -r 50 --configfile $CONFIG  --datatable data_XM13A_nacncoop_pasdend
# done
# echo $RUNTEXT
# for f in $CELLNAMES
# do
#     echo $f
#     python vcnmodel/model_run2.py VCN_c$f  -F -P initAN --configfile $CONFIG  --datatable data_XM13A_nacncoop_actdend
#     python vcnmodel/model_run2.py VCN_c$f  -F -P runANPSTH -r 50 --configfile $CONFIG  --datatable data_XM13A_nacncoop_actdend
# done
echo AN runs complete
# with "A", we do all cells in grade A
# python vcnmodel/plotters/plot_gbc_ivresults_2.py "A" -s
