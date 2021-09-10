echo " CHI2_before CHI2_after PHI_eff" > GRID_fill

mkdir WEIGHTS

for i in `grep -v "#" GRID | awk '{print $1}'`
do
        if [ -f GP$i/gp$(($i))_ibme_19.log ]
        then
                echo `more GP$i/gp$(($i))_ibme_19.log | grep 'CHI2 before optimization:\|CHI2 after optimization:\|Fraction of effective frames:' | awk '{printf $NF" "}' ; echo " " | awk '{printf "\n"}'` >> GRID_fill
               mv GP$i/gp$(($i))_19.weights.dat WEIGHTS/w$i.dat 
        else
                echo "inf inf inf" >> GRID_fill
        fi
done

paste GRID GRID_fill > GRID_opt
rm GRID_fill
