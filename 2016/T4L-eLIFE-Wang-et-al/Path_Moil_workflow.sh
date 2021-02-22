#!/bin/sh
# wyong 2007.1.30
log=0
grid=31
gamalist="1 2 10 100 1000 10000"
repllist="1 10 100 1000 4000 40000"
date |awk '{print $4}' |awk -F: '{print $1,$2,$3}' >time1
for gama in $gamalist
do
for repl in $repllist 
do
log=`expr $log + 1`
echo '####################################################'
echo "gama= $gama || repl= $repl || logfile=chmin.log.$log" 
cat <<EOF >chmin${log}.inp
file conn name=(T4L_G.wcon) read
file rcrd unit=5 read
file wcrd name=(chmin_${log}.pth) bina wovr
#ste=400 #pri=100 #wcr=100 list=100
gama=$gama repl=$repl
lmbd=2.
grid=$grid
relx=12. rvmx=9. epsi=1. cdie cini
hvdw
~gbsa
action
file name=(T4L_G_min.crd) read
file name=(T4L_E_min.crd) read
EOF
echo "chmin <chmin${log}.inp >chmin.log.$log ......"
/nike-storage/wyong/software/moil11/moil.source/exe/chmin <chmin${log}.inp >chmin.log.$log

echo '####################################################'

#==========================================
# Step2. do SDP optimization
cat <<EOF >sdp${log}.inp
file conn name=(T4L_G.wcon) read
file rcrd name=(chmin_${log}.pth) bina unit=14 read
file wcrd name=(sdp_${log}.pth) bina unit=12 wovr
#ste=1000 #pri=100 #wcr=100 list=200 tolg=0.001 estred=0.1
gama=200 grid=31 hami=1.d-5 
pdqe=-42.2
rmax=9999. epsi=1. cdie v14f=8. el14=2. cpth
tmpr=300.0 proc=1 nupd=1
~gbsa
hvdw
action
EOF
/nike-storage/wyong/software/moil11/moil.source/exe/sdpS <sdp${log}.inp >sdp.log.${log}

#=========================================
# Step3. convert pth to pdb structures
bash pth2pdb.sh sdp_${log}.pth pdb${log}


#=========================================
# Step4. combine pdb structures into one pdb file 
sh pdb2multimodel.sh pdb${log}

#tail -40 chmin.log.$log |head -33 >temp
#awk 'NR>1 && NR<33 {print $1,$4}' temp >$gama.$repl.dat
#echo "$gama $repl" >>$gama.$repl.dat
done
#paste $gama.*.dat |awk '{printf("%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f \n",$2,$4,$6,$8,$10,$12)}' >$gama.dat
done
date |awk '{print $4}' |awk -F: '{print $1,$2,$3}' >time2
paste time1 time2 |awk '{printf("The total time is :\n %3d hours,%3d minutes,%3d seconds\n",$4-$1,$5-$2,$6-$3)}'
