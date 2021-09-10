structures=$1
exp_path=$2

theta=$3

gl=$4

ens_size=`find $structures/ -name "frame*pdb" | wc -l`

grid_point=`grep -v "#" GRID_FoXS | sed -n $gl'p' | awk '{print $1}'`
echo "FoXS in GP"$grid_point
mkdir GP$grid_point
c2=`grep -v "#" GRID | sed -n $gl'p' | awk '{print $2}'`
c1=`grep -v "#" GRID | sed -n $gl'p' | awk '{print $3}'`
for i in `seq 0 $(($ens_size-1))`
do
    cp $structures/frame$i.pdb GP$grid_point/frame.pdb
    foxs GP$grid_point/frame.pdb $exp_path --min_c1 $c1 --max_c1 $c1 --min_c2 $c2 --max_c2 $c2
    grep -v "#" GP$grid_point/frame*.fit | awk -v i=$i 'BEGIN{printf i" "}{printf $4" "}END{printf "\n"}' >> GP$grid_point/calc_saxs.txt
done > GP$grid_point/logFoXS

python /groups/sbinlab/fpesce/scripts/HYDRA_FoXS/rescale.py $grid_point

echo "Starting iBME optimization in GP"$grid_point
python /groups/sbinlab/fpesce/scripts/BME2/iBME.py $exp_path GP$grid_point/calc_saxs_res.txt $theta gp$grid_point
mv gp$(($grid_point))_* GP$grid_point/
mv gp$(($grid_point)).log GP$grid_point/
echo "DONE GP"$grid_point

echo
echo "ENDING TIME: " `date`
