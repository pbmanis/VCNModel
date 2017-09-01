#!/bin/bash

python model_run.py VCN_c08 -M mGBC -S MS -r 20 -P runANIO --noparallel &
python model_run.py VCN_c09 -M mGBC -S MS -r 20 -P runANIO --noparallel & 
python model_run.py VCN_c09h -M mGBC -S MS -r 20 -P runANIO --noparallel &
python model_run.py VCN_c09nd -M mGBC -S MS -r 20 -P runANIO --noparallel &
python model_run.py VCN_c17  -M mGBC -S MS -r 20 -P runANIO --noparallel &
python model_run.py VCN_c18  -M mGBC -S MS -r 20 -P runANIO --noparallel &
python model_run.py VCN_c19  -M mGBC -S MS -r 20 -P runANIO --noparallel &
python model_run.py VCN_c20 -M mGBC -S MS -r 20 -P runANIO --noparallel &
python model_run.py VCN_c21  -M mGBC -S MS -r 20 -P runANIO --noparallel &
python model_run.py VCN_c22  -M mGBC -S MS -r 20 -P runANIO --noparallel &
