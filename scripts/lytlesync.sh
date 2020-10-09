#!/bin/bash
# Use this to kill stopped zombies processes:
# ps -x | grep -0 "/usr/local/bin/rsync"
# then kill the process
TIME=$1
DRYRUN=$2
SIM="AN"
if [[ -z $TIME ]]; then
  echo "Error: no time argument."
  echo "Please enter the number of days to sync."
  exit 1
fi

if [[ $DRYRUN = "dry" ]]; then
  DRYRUNCMD="--dry-run"
  echo "Dry run initiated..."
fi

CELLNO="06 09 10 11 13 17 30"
# -exec ls $(basename {}) \;") \
#-exec ls $(basename {}) \;
# do a retrieval.
for cn in $CELLNO
do
    datapath="/Users/pbmanis/Desktop/Python/VCN-SBEM-Data/VCN_Cells/VCN_c$cn/Simulations/$SIM/"
    echo From Cell VCN_$cn
    echo $FROMFI
    /usr/local/bin/rsync -avrz -R \
        --exclude 'archived simulations' \
        --exclude 'archived simulations/**' \
        $DRYRUNCMD --files-from=<(ssh \
        pbmanis@152.19.86.111 "find $datapath -type d -name '*2020-09-02*'\
        ") \
        pbmanis@152.19.86.111:/ /
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
