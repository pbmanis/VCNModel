# Example:
# scripts/process_gbcIV.sh run all 
#

#######################################################
# Full models are from data/reconstuctions Matthew Kersting sent on
# March 6, 2020. 
# Note we do not have a full reconstruction for cell 18
# in that dataset.
#######################################################
CELLNAMES="17" # "02 05 06 09 " # 
#CONFIG="noscale.toml" #"autoscale.toml"
CONFIG1="autoscale_xm13a_multisite_parallel_SAM_mixed.toml"
RUNTEXT1="running the individual initialization and running AN PSTH protocols"
CONFIG2="autoscale_xm13a_multisite_parallel_SAM_allHS.toml"
RUNTEXT2="running the individual initialization and running AN PSTH protocols"
NREPS="20"
FMOD="400"
DEP=1

echo $RUNTEXT1
for f in $CELLNAMES
do
    echo $f
    # python vcnmodel/model_run2.py VCN_c$f  -F -P initAN --configfile $CONFIG
    # if [ $? -ne 0 ]; then
    #     exit 1
    # fi
    # python vcnmodel/model_run2.py VCN_c$f  -D Full -P runANPSTH -r 1 --dB 10 --Spirou twolargest --configfile $CONFIG --datatable data_XM13A_nacncoop
    # if [ $? -ne 0 ]; then
    #     exit 1
    # fi
    python vcnmodel/model_run2.py VCN_c$f  -D Full -P runANPSTH -r $NREPS --dB 20 --fmod $FMOD --depression $DEP --Spirou all --configfile $CONFIG1 --datatable data_XM13A_nacncoop
    if [ $? -ne 0 ]; then
        exit 1
    fi
    python vcnmodel/model_run2.py VCN_c$f  -D Full -P runANPSTH -r $NREPS --dB 50 --fmod $FMOD --depression $DEP --Spirou all --configfile $CONFIG1 --datatable data_XM13A_nacncoop
    if [ $? -ne 0 ]; then
        exit 1
    fi
    python vcnmodel/model_run2.py VCN_c$f  -D Full -P runANPSTH -r $NREPS --dB 80 --fmod $FMOD --depression $DEP --Spirou all --configfile $CONFIG1 --datatable data_XM13A_nacncoop
    if [ $? -ne 0 ]; then
        exit 1
    fi
    python vcnmodel/model_run2.py VCN_c$f  -D Full -P runANPSTH -r $NREPS --dB 20 --fmod $FMOD --depression $DEP --Spirou all --configfile $CONFIG2 --datatable data_XM13A_nacncoop
    if [ $? -ne 0 ]; then
        exit 1
    fi
    python vcnmodel/model_run2.py VCN_c$f  -D Full -P runANPSTH -r $NREPS --dB 50 --fmod $FMOD --depression $DEP --Spirou all --configfile $CONFIG2 --datatable data_XM13A_nacncoop
    if [ $? -ne 0 ]; then
        exit 1
    fi
    python vcnmodel/model_run2.py VCN_c$f  -D Full -P runANPSTH -r $NREPS --dB 80 --fmod $FMOD --depression $DEP --Spirou all --configfile $CONFIG2 --datatable data_XM13A_nacncoop
    if [ $? -ne 0 ]; then
        exit 1
    fi
done

wait

echo AN runs complete
# with "A", we do all cells in grade A
# python vcnmodel/plotters/plot_gbc_ivresults_2.py "A" -s