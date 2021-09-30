#!/bin/bash

for i in hSUMO_hnRNPA1 hnRNPA1 TIA1
do

cd $i

for j in 1.00 1.04 1.06 1.08 1.10 1.12 1.14
do

cd lambda_${j}
cp ../../prodrun_mdrun.sh .

echo "#!/bin/sh" > temp
echo "#PBS -W group_list=ku_10001 -A ku_10001" >> temp
echo "#PBS -N ${i}_${j}_md" >> temp
cat prodrun_mdrun.sh >> temp
mv temp prodrun_mdrun.sh

qsub prodrun_mdrun.sh

cd ..

done

cd ..

done

