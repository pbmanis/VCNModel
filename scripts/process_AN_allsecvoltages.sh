#!/bin/bash
# python vcnmodel/model_run.py VCN_c06 --type Bushy --modeltype II --protocol initIV --hocfile VCN_c06_Full.hocx --model XM13_nacncoop \
#          -r 1 \
#          --tagstring noDendNa\
#          --saveall --noparallel\
#          #--soma_autoinflate --dendrite_autoinflate
#
# python vcnmodel/model_run.py VCN_c06 --type Bushy --modeltype II --protocol runIV --hocfile VCN_c06_Full.hocx --model XM13_nacncoop \
#          -r 1 \
#          --tagstring noDendNa\
#          --saveall --noparallel\
#         #--soma_autoinflate --dendrite_autoinflate
#
# echo IV saveallrun complete

python vcnmodel/model_run.py VCN_c06 --type Bushy --modeltype II --protocol initIV --hocfile VCN_c06_Full.hocx --model XM13_nacncoop \
         -r 1 \
         --tagstring noDendNa_Full
         --saveall --noparallel\
         #--soma_autoinflate --dendrite_autoinflate

python vcnmodel/model_run.py VCN_c06 --type Bushy --modeltype II --protocol runIV --hocfile VCN_c06_Full.hocx --model XM13_nacncoop \
         -r 1 \
         --tagstring noDendNa_Full\
         --saveall --noparallel\
        #--soma_autoinflate --dendrite_autoinflate

echo IV saveallrun complete

# python vcnmodel/model_run.py VCN_c06 --type Bushy --modeltype II --protocol initAN --hocfile VCN_c06_Full.hocx --model XM13_nacncoop \
#          -r 1 --sgcmodel cochlea -S MS -f 16000. \
#          --saveall --noparallel \
#          --tagstring noDendNa \
#              # --soma_autoinflate --dendrite_autoinflate\
#     #
#
#
# python vcnmodel/model_run.py VCN_c06 --type Bushy --modeltype II --protocol runANPSTH --hocfile VCN_c06_Full.hocx --model XM13_nacncoop \
#          -r 1 --sgcmodel cochlea -S MS -f 16000. \
#          --tagstring noDendNa \
#          --saveall --noparallel \
#         #     --soma_autoinflate --dendrite_autoinflate\
#
# wait
# echo AN saveallrun complete


