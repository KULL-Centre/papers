#!/bin/bash

for i in Hst5_PME aSyn_PME
do

cd $i

for j in 1.00 1.04 1.06 1.08 1.10 1.12 1.14
do

cd lambda_${j}
cp ../../relax_mdrun_PME.sh .

echo "#!/bin/sh" > temp
echo "#PBS -W group_list=ku_10001 -A ku_10001" >> temp
echo "#PBS -N ${i}_${j}_relax" >> temp
cat relax_mdrun_PME.sh >> temp
mv temp relax_mdrun_PME.sh

qsub relax_mdrun_PME.sh

cd ..

done

cd ..

done

