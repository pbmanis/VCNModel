# Example:
# scripts/process_gbcIV.sh run all 
#

#######################################################
# Full models are from data/reconstuctions Matthew Kersting sent on
# March 6, 2020. 
# Note we do not have a full reconstruction for cell 18
# in that dataset.
#######################################################
CELLNAMES="11" # "02 05 06 09 10 11 13 17 18 30"
CONFIG="--configfile singles_multisite_parallel_regularclicktrain_40dB.toml"
RUNTEXT="running the individual initialization and running AN SAM protocols"
TABLE="--datatable data_XM13A_nacncoop_normal" # " data_XM13A_nacncoop_pasdend data_XM13A_nacncoop_actdend"
REPS="-r 2"
DMOD="--dmod 100."
WKN="--workers 16"
TEST="" # "--testsetup"
EXP1="--Spirou all" # all=mean max=mean removelargest largestonly twolargest"
EXP2="--Spirou largestonly" # all=mean max=mean removelargest largestonly twolargest"
EXP3="--Spirou removelargest" # all=mean max=mean removelargest largestonly twolargest"
TESTSCRIPT="False"  # or "False"
echo $RUNTEXT
for cell in $CELLNAMES
do
    # Handle case when simulations failed because of disk space  
    # These statements test cell ID and frequency.
    # you can verify what will be run by seting TESTSCRIPT to "True"
    # if [[ ($cell = "02" || $cell = "05") && $freq = "500" ]]
    # then
    #     echo Not running VCN_c$cell $freq
    #     continue
    # fi
    # if [[ ($freq -le "400") ]]  # account for runs all completed for lower frqeuencies
    # then
    #     echo Skipping VCN_c$cell $freq
    #     continue
    # fi
    #
    # echo Running Cell VCN_c$cell $freq
    #
    # AXON=""
    # case $cell in
    #     02 | 05)
    #         AXON="-A standardized"
    #         echo "   " Using standardized axon on cell
    #         ;;
    #
    # esac
    if [[ ($TESTSCRIPT = "True") ]]
    then
        echo "       (testing)"
        continue
    fi
    python vcnmodel/model_run2.py VCN_c$cell -D Full $AXON -P initAN $CONFIG $TABLE $TEST
    python vcnmodel/model_run2.py VCN_c$cell -D Full $AXON -P runANPSTH $REPS $EXP1 $CONFIG $TABLE $TEST $WKN
    if [ $? -ne 0 ]; then
        exit 1
    fi
    # python vcnmodel/model_run2.py VCN_c$cell -D Full $AXON -P runANPSTH $REPS $EXP2 $CONFIG $TABLE $TEST $WKN
    # if [ $? -ne 0 ]; then
    #     exit 1
    # fi
    # python vcnmodel/model_run2.py VCN_c$cell -D Full $AXON -P runANPSTH $REPS $EXP3 $CONFIG $TABLE $TEST $WKN
    # if [ $? -ne 0 ]; then
    #     exit 1
    # fi

done


wait

echo AN frozen noise runs complete
# with "A", we do all cells in grade A
# python vcnmodel/plotters/plot_gbc_ivresults_2.py "A" -s
