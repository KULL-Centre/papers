#!/bin/sh
#
if [ $# -lt 2 ]; then
	echo "$0 pdb1 model_1.pdb"
	exit
else
	dir=$1
	output=$2
fi
rm $output 2>&1 >/dev/null
for i in `seq 1 1 31`
do
j=`printf %03d $i`
echo "MODEL $i" >>$output
cat $dir/${j}.pdb >> $output
echo "ENDMDL" >>$output
done
