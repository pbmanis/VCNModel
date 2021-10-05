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
FREQS="100" # 200 300 400 500 750 1000"
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
        AXON=""
        case $cell in
            02 | 05)
                AXON="-A standardized"
                echo Using standardized axon on cell
                ;;

        esac
        echo Cell VCN_c$cell
        python src/vcnmodel/model_run2.py VCN_c$cell -D Full $AXON -P initAN $CONFIG $TABLE $TEST
        python src/vcnmodel/model_run2.py VCN_c$cell -D Full $AXON -P runANPSTH --fmod $freq $DMOD $REPS $EXP1 $CONFIG $TABLE $TEST $WKN
        if [ $? -ne 0 ]; then
            exit 1
        fi
        python src/vcnmodel/model_run2.py VCN_c$cell -D Full $AXON -P runANPSTH --fmod $freq $DMOD $REPS $EXP2 $CONFIG $TABLE $TEST $WKN
        if [ $? -ne 0 ]; then
            exit 1
        fi
        python src/vcnmodel/model_run2.py VCN_c$cell -D Full $AXON -P runANPSTH --fmod $freq $DMOD $REPS $EXP3 $CONFIG $TABLE $TEST $WKN
        if [ $? -ne 0 ]; then
            exit 1
        fi

    done
done

wait

echo AN SAM runs complete
# with "A", we do all cells in grade A
# python vcnmodel/plotters/plot_gbc_ivresults_2.py "A" -s
