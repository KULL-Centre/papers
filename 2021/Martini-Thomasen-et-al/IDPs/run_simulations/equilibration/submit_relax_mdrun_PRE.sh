#!/bin/bash

for i in htau40 OPN aSyn_PRE A2 #FUS
do

cd ${i}

for j in 1.00 1.10 1.12
do

cd lambda_${j}
cp ../../relax_mdrun.sh .

echo "#!/bin/sh" > temp
echo "#PBS -W group_list=ku_10001 -A ku_10001" >> temp
echo "#PBS -N ${i}_${j}_relax" >> temp
cat relax_mdrun.sh >> temp
mv temp relax_mdrun.sh

qsub relax_mdrun.sh

cd ..

done

cd ..

done

