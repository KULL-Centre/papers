pepsi_path=/groups/sbinlab/fpesce/software/Pepsi-SAXS
structures=$1
exp_path=$2

theta=$3

ens_size=`find $structures/ -name "frame*pdb" | wc -l`

gl=$4

grid_point=`grep -v "#" GRID | sed -n $gl'p' | awk '{print $1}'`
echo "Pepsi-SAXS in GP"$grid_point
mkdir GP$grid_point
mkdir GP$grid_point/pepsi
dro=`grep -v "#" GRID | sed -n $gl'p' | awk '{print $2}'`
r0=`grep -v "#" GRID | sed -n $gl'p' | awk '{print $3}'`
for i in `seq 0 $(($ens_size-1))`
do
	$pepsi_path $structures/frame$i.pdb $exp_path -o GP$grid_point/pepsi/saxs$i.dat -cst --cstFactor 0 --I0 1.0 --dro $dro --r0_min_factor $r0 --r0_max_factor $r0 --r0_N 1
	#echo `cat GP$grid_point/pepsi/saxs$i.dat | awk '{if (NF==4) printf $4" "}' ; echo " " | awk '{printf "\n"}'` >> GP$grid_point/calc_saxs.txt
	#rm GP$grid_point/pepsi/saxs$i.dat
	#mv GP$grid_point/pepsi/saxs$i.log GP$grid_point/pepsi.log	
done > GP$grid_point/pepsi/logPEPSI

echo "Packing Pepsi-SAXS output as BME input in GP"$grid_point
python3.8 /groups/sbinlab/fpesce/scripts/BME/bme_prep.py GP$grid_point/pepsi $ens_size GP$grid_point
echo "Cleaning Pepsi-SAXS data in GP"$grid_point
mv GP$grid_point/pepsi/saxs0.log GP$grid_point/pepsi.log
rm -r GP$grid_point/pepsi &

echo "Starting iBME optimization in GP"$grid_point
python3.8 /groups/sbinlab/fpesce/scripts/BME2/iBME.py $exp_path GP$grid_point/calc_saxs.txt $theta gp$grid_point
mv gp$(($grid_point))_* GP$grid_point/
mv gp$(($grid_point)).log GP$grid_point/
echo "DONE GP"$grid_point

