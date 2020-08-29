#!/bin/bash
# Use this to kill stopped zombies processes:
# ps -x | grep -0 "/usr/local/bin/rsync"
# then kill the process
TIME=$1
DRYRUN=$2
SIM="AN"
datapath="/Users/pbmanis/Desktop/Python/VCN-SBEM-Data/VCN_Cells/VCN_c$cn/Simulations/$SIM/"
if [[ -z $TIME ]]; then
  echo "Error: no time argument."
  echo "Please enter the number of days to sync."
  exit 1
fi

if [[ $DRYRUN = "dry" ]]; then
  DRYRUNCMD="--dry-run"
  echo "Dry run initiated..."
fi
CELLNO="30"

# do a retrieval.
for cn in $CELLNO
do
    datapath="/Users/pbmanis/Desktop/Python/VCN-SBEM-Data/VCN_Cells/VCN_c$cn/Simulations/$SIM/"
    echo From Cell VCN_$cn
    /usr/local/bin/rsync -avrz -v \
        $DRYRUNCMD --files-from=<(ssh \
        pbmanis@152.19.86.111 "find $datapath \
        -mtime -$TIME -type d -exec ls $(basename {}) \;") \
        pbmanis@152.19.86.111:$datapath/ $datapath
        echo " "
done
echo Finished sync

# just check script against local files
# for cn in $CELLNO
# do
#     datapath="/Users/pbmanis/Desktop/Python/VCN-SBEM-Data/VCN_Cells/VCN_c$cn/Simulations/$SIM"
#     echo VCN_$cn
#     find $datapath \
#         -mtime -$TIME -type f \
#         -exec ls $(basename {}) \; | cut -d '/' -f 11-
#         #| sed -e 's|^$datapath/*||p' \;
#         # | cut -d '/' -f 2
#     echo " "
# done

# check script against remote files
# for cn in $CELLNO
# do
#         datapath="/Users/pbmanis/Desktop/Python/VCN-SBEM-Data/VCN_Cells/VCN_c$cn/Simulations/$SIM"
#         ssh \
#         pbmanis@152.19.86.111 "find $datapath \
#         -mtime -$TIME -type f \
#         | cut -d '\' -f 2- \
#         -exec ls -lat $(basename {}) \;"
#         echo " "
# done
