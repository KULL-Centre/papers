#!/bin/bash


cd two_villin_h36_init

for j in 1.00
do

cd lambda_${j}
cp ../../relax_mdrun.sh .

echo "#!/bin/sh" > temp
echo "#PBS -W group_list=ku_10001 -A ku_10001" >> temp
echo "#PBS -N villin_${j}_relax" >> temp
cat relax_mdrun.sh >> temp
mv temp relax_mdrun.sh

qsub relax_mdrun.sh

cd ..

done

cd ..


