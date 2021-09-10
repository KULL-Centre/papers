pepsi_path=/groups/sbinlab/fpesce/software/Pepsi-SAXS
structures=$1
exp_path=$2

theta=$3

ens_size=`find $structures/ -name "frame*pdb" | wc -l`

gl=$4

grid_point=`grep -v "#" GRID_Pepsi-SAXS | sed -n $gl'p' | awk '{print $1}'`
echo "Pepsi-SAXS in GP"$grid_point
mkdir GP$grid_point
dro=`grep -v "#" GRID | sed -n $gl'p' | awk '{print $2}'`
r0=`grep -v "#" GRID | sed -n $gl'p' | awk '{print $3}'`
for i in `seq 0 $(($ens_size-1))`
do
    $pepsi_path $structures/frame$i.pdb $exp_path -o GP$grid_point/saxs$i.dat -cst --cstFactor 0 --I0 1.0 --dro $dro --r0_min_factor $r0 --r0_max_factor $r0 --r0_N 1
    grep -v "#" GP$grid_point/saxs$i.dat | awk -v i=$i 'BEGIN{printf i" "}{printf $4" "}END{printf "\n"}' >> GP$grid_point/calc_saxs.txt
    rm GP$grid_point/saxs$i.dat
    if [ $grid_point -eq 0 ]
    then
        more GP$grid_point/saxs$i.log | grep "Radius of gyration of the envelope" | awk '{print $8}' >> GP$grid_point/Rg_env.dat
    fi
    rm GP$grid_point/saxs$i.log
done > GP$grid_point/logPEPSI

echo "Starting iBME optimization in GP"$grid_point
python3.8 /groups/sbinlab/fpesce/scripts/BME2/iBME.py $exp_path GP$grid_point/calc_saxs.txt $theta gp$grid_point
mv gp$(($grid_point))_* GP$grid_point/
mv gp$(($grid_point)).log GP$grid_point/
echo "DONE GP"$grid_point

