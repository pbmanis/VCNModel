CELLS="02 05 06 09 10 11 13 17 18 30"
DIR="/Users/pbmanis/Desktop/Python/VCN-SBEM-Data" # "/VCN_Cells/"
REMOTESYSTEM="pbmanis@lytle.med.unc.edu"
REMOTEDIR="pbmanis@lytle.med.unc.edu:"
FLAGS="-RDa0P -r -v" #" --dry-run"
EXC="--exclude=Z*"
# RS1="--files-from=<(ssh pbmanis@lytle.med.unc.edu find /Users/pbmanis/Desktop/Python/VCN-SBEM-DATA/VCN_Cells/ -mtime -1 -print0)"
SIM="VC" # "IV"
#
# The goal here is to only check the relevant directories,
# and only for data from the last day

cd $DIR
touch files.txt
"" > files.txt

for cell in $CELLS
do
	echo "getting remote file list: " VCN_c$cell $FLAGS
	# RCMD=${RS1}VCN_c$cell/Simulations/$SIM
	# echo $RCMD
	ssh $REMOTESYSTEM "cd /Users/pbmanis/Desktop/Python/VCN-SBEM-Data; find -f VCN_Cells/VCN_c$cell/Simulations/$SIM/* -mtime -400 -print0" > files.txt
#	ssh $REMOTESYSTEM "cd /Users/pbmanis/Desktop/Python/VCN-SBEM-Data; find -f VCN_Cells/VCN_c$cell/Morphology/* -mtime -800 -print0" > files.txt
	rsync $FLAGS  --files-from=files.txt $REMOTESYSTEM:$DIR .
done

# less files_VCN_c02.txt
# # exit
# # # wait
# # #
# # for cell in $CELLS
# # do
# # 	echo "retrieving for " $cell $FLAGS
# # 	# RCMD=${RS1}VCN_c$cell/Simulations/$SIM
# # 	# echo $RCMD
# # 	# ssh $REMOTESYSTEM "cd /Users/pbmanis/Desktop/Python/VCN-SBEM-Data; find -f VCN_Cells/VCN_c$cell/Simulations/$SIM/* -mtime -2 -print0"> files.txt
# # 	rsync $FLAGS  --files-from=files.txt $REMOTESYSTEM:$DIR . &
# # 	# rsync $FLAGS --files-from=file.tt -mtime -1 -print0") $DIR/VCN_c$cell/Simulations/IV
# # done
# # wait