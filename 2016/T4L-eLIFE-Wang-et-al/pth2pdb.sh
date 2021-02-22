#!/bin/sh
# wyong 2009.08.19
#
if [ $# -lt 1 ]; then
	echo "$0 sdp_bind_ENM.pth path1"
	exit
else
	pth=$1
	dir=$2
fi
bindir="/nike-storage/wyong/software/moil11/moil.source/exe"
wcon="T4L_G.wcon"
all=31

mkdir $dir 2>/dev/null
i=0
while [ $i -lt $all ]
do
	let i=$i+1
	echo $i
cat << EOF > .inp
~  basic ccrd input file
~
file conn name=($wcon) unit=10 read
file rcrd name=($pth) bina unit=11 read
file wcrd name=(foo.crd) unit=12 wovr
fpth tchr lpst=$i
action
*EOD
EOF
	$bindir/ccrd < .inp > .log
	$bindir/crd2pdb < foo.crd >foo.pdb
	ipdb=`printf %03d $i`
	awk '$0 !~ /[NC]TER/' foo.pdb > $dir/$ipdb.pdb
done

