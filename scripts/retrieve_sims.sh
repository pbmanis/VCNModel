CELLS="02" #" 05 06 09 10 11 13 17 18 30"
DIR="/Users/pbmanis/Desktop/Python/VCN-SBEM-DATA/VCN_Cells/"
REMOTESYSTEM="ssh pbmanis@lytle.med.unc.edu"
REMOTEDIR="pbmanis@lytle.med.unc.edu:"
FLAGS="-RDa0P -rr --dry-run"
EXC="--exclude=Z*"
RS1="--files-from=<(ssh pbmanis@lytle.med.unc.edu find /Users/pbmanis/Desktop/Python/VCN-SBEM-DATA/VCN_Cells/ -mtime -1 -print0)"

#
# The goal here is to only check the relevant directories,
# and only for data from the last day

for cell in $CELLS
do
	echo $FLAGS
	RCMD=${RS1}VCN_c$cell/Simulations/IV
	# echo $RCMD
	rsync $FLAGS --files-from=<(ssh pbmanis@lytle.med.unc.edu "find /Users/pbmanis/Desktop/Python/VCN-SBEM-DATA/VCN_Cells/ -mtime -1 -print0") $DIR/VCN_c$cell/Simulations/IV
done
wait
