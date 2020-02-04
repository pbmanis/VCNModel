#!/bin/bash
python vcnmodel/model_run.py VCN_c09 --type Bushy --modeltype II --protocol initAN --hocfile VCN_c09_FullCell_edited.hoc --model XM13_nacncoop \
         -r 1 --sgcmodel cochlea -S MS -f 16000. --soma_autoinflate --dendrite_autoinflate\
         --tagstring noDendNa\
         --saveall --noparallel


python vcnmodel/model_run.py VCN_c09 --type Bushy --modeltype II --protocol runANPSTH --hocfile VCN_c09_FullCell_edited.hoc --model XM13_nacncoop \
         -r 1 --sgcmodel cochlea -S MS -f 16000. --soma_autoinflate --dendrite_autoinflate\
         --tagstring noDendNa\
         --saveall --noparallel

wait
echo AN saveallrun complete


