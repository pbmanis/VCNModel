# Example:
# scripts/process_gbcIV.sh run all 
#

#######################################################
# Full models are from data/reconstuctions Matthew Kersting sent on
# March 6, 2020. 
# Note we do not have a full reconstruction for cell 18
# in that dataset.
#######################################################
CELLNAMES="02 05 06 09 10 11 13 17 18 30"
# CELLNAMES="09 11 17"  # just for twolargest - only those for which second input is also suprathreshold
#CONFIG="noscale.toml" #"autoscale.toml"
CONFIG="--configfile singles_multisite_parallel_SAM30dB.toml"
RUNTEXT="running the individual initialization and running AN SAM protocols"
TABLE="--datatable data_XM13A_nacncoop_normal" # " data_XM13A_nacncoop_pasdend data_XM13A_nacncoop_actdend"
REPS="-r 100"
FREQS="50" # 100 200 300 400 500 750 1000"
DMOD="--dmod 100."
WKN="--workers 16"
TEST="" #"--testsetup"
EXP1="--Spirou all" # all=mean max=mean removelargest largestonly twolargest"
EXP2="--Spirou largestonly" # all=mean max=mean removelargest largestonly twolargest"
EXP3="--Spirou removelargest" # all=mean max=mean removelargest largestonly twolargest"
echo $RUNTEXT
for freq in $FREQS
do
    echo "Freq: " $freq
    for cell in $CELLNAMES
    do
        # Handle case when simulations failed because of disk space
        # Generally, this will not be needed, but just in case... 
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
        # echo Running Cell VCN_c$cell $freq
        #
        
        ################################################################
        # DO NOT DELETE THE FOLLOWING SECTION                          #
        # This makes sure we run cells 2 and 5 in the appropriate mode #
        ################################################################
        AXON=""
        case $cell in
            02 | 05)
                AXON="-A standardized"
                echo Using standardized axon on cell
                ;;
        esac
        ################################################################
        # DO NOT DELETE THE ABOVE SECTION                              #
        # This makes sure we run cells 2 and 5 in the appropriate mode #
        ################################################################
        if [[ ($TESTSCRIPT = "True") ]]
        then
            echo "       (testing)"
            continue
        fi

        python vcnmodel/model_run2.py VCN_c$cell -D Full $AXON -P initAN $CONFIG $TABLE $TEST
        python vcnmodel/model_run2.py VCN_c$cell -D Full $AXON -P runANPSTH --fmod $freq $DMOD $REPS $EXP1 $CONFIG $TABLE $TEST $WKN
        if [ $? -ne 0 ]; then
            exit 1
        fi
        python vcnmodel/model_run2.py VCN_c$cell -D Full $AXON -P runANPSTH --fmod $freq $DMOD $REPS $EXP2 $CONFIG $TABLE $TEST $WKN
        if [ $? -ne 0 ]; then
            exit 1
        fi
        python vcnmodel/model_run2.py VCN_c$cell -D Full $AXON -P runANPSTH --fmod $freq $DMOD $REPS $EXP3 $CONFIG $TABLE $TEST $WKN
        if [ $? -ne 0 ]; then
            exit 1
        fi

    done
done

wait

echo AN SAM runs complete
# with "A", we do all cells in grade A
# python vcnmodel/plotters/plot_gbc_ivresults_2.py "A" -s
