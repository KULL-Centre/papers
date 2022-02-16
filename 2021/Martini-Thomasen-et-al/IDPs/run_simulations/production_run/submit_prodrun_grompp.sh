#!/bin/bash

for i in ACTR_helices aSyn_bigboxtest K25 A1 CoRNID ColNT FhuA Hst52 K19 PNt Sic1 aSyn Hst5 ACTR
do

cd $i

t=0

if [[ "$i" == "ColNT" ]]
then
        t=277
fi

if [[ "$i" == "ACTR" || "$i" == "ACTR_helices" ]]
then
        t=278
fi

if [[ "$i" == "K25" || "$i" == "K19" ]]
then
        t=288
fi

if [[ "$i" == "Hst5" || "$i" == "Sic1" || "$i" == "aSyn" || "$i" == "CoRNID" || "$i" == "aSyn_bigboxtest" ]]
then
        t=293
fi

if [[ "$i" == "A1" ]]
then
        t=296
fi

if [[ "$i" == "Hst52" || "$i" == "FhuA" || "$i" == "PNt" || "$i" == "FUS" || "$i" == "A2" ]]
then
        t=298
fi

echo "$i is at ${t}K"

for j in 1.00 1.04 1.06 1.08 1.10 1.12 1.14
do

cd lambda_${j}
cp ../../prodrun_grompp.sh .
qsub prodrun_grompp.sh -v temp=$t
cd ..

done

cd ..

done

