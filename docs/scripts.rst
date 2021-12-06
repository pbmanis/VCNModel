Scripts
=======

Many simulations with parametric variation are controlled by bash (or zsh) scripts, that in turn set the parameters for a given simulation.
The scripts can have nested loops for exploration of parameter spaces, and provide for the generation of initialization files specific
to the simulations that are being run. In general, the scripts should not be "shelled" out (with an '&') as subsequent steps may depend
on completion of prior steps such as initialization. Parallelization is handled within model_run2.py, where repeated simulations varying
only by a single parameter or selection of a random number generator (for auditory nerve bases simulations) use the parallel processing
capabilities of a given computer. 

Scripts should generally take no arguments, and should be saved in the "scripts" directory.

A simple example script that runs the IV protocol for all cells, and for 3 different decorations (gradeA_IVs.sh) looks like this::

    # Example:
    # scripts/process_gbcIV.sh
    #
    PROTO="runIV"
    #######################################################
    # Full models are from data/reconstuctions Matthew Kersting sent on
    # March 6, 2020. 
    # 
    #######################################################
    CELLNAMES="02 05 06 09 10 11 13 17 18 30"
    CONFIG="--configfile xm13a_multisite_parallel.toml"
    DATATABLE="--datatable data_XM13A_nacncoop_normal"
    DATATABLEA="--datatable data_XM13A_nacncoop_actdend"
    DATATABLEP="--datatable data_XM13A_nacncoop_pasdend"

    echo "running the individual initialization and/or running of IV protocols"
    for f in $CELLNAMES
    do
        AXON=""
        case $f in
            02 | 05)
                AXON="-A standardized"
                ;;
        esac
        echo $f
        echo $AXON
        python src/vcnmodel/model_run2.py VCN_c$f -P initIV $CONFIG -D Full $AXON $DATATABLE --dendritemode normal
        python src/vcnmodel/model_run2.py VCN_c$f -P runIV  $CONFIG -D Full $AXON $DATATABLE --dendritemode normal
    done
    wait

    for f in $CELLNAMES
    do
        AXON=""
        case $f in
            02 | 05)
                AXON="-A standardized"
                ;;
        esac
        echo $f
        echo $AXON
        python src/vcnmodel/model_run2.py VCN_c$f -P initIV $CONFIG -D Full $AXON $DATATABLEP --dendritemode passive
        python src/vcnmodel/model_run2.py VCN_c$f -P runIV  $CONFIG -D Full $AXON $DATATABLEP --dendritemode passive
    done
    wait

    for f in $CELLNAMES
    do
        AXON=""
        case $f in
            02 | 05)
                AXON="-A standardized"
                ;;
        esac
        echo $f
        echo $AXON
        python src/vcnmodel/model_run2.py VCN_c$f -P initIV $CONFIG -D Full $AXON $DATATABLEA --dendritemode active
        python src/vcnmodel/model_run2.py VCN_c$f -P runIV  $CONFIG -D Full $AXON $DATATABLEA --dendritemode active
    done
    wait

    echo IV runs complete


6 December 2021 pbm