#!/bin/sh
#

parm="../EXG.explicit.parm7"
cpin="../EXG.cpin"
beginstr="../EXG.equil.rst7"

Ngpu=4
Ntasks=$((Ngpu*5))

rm groupfile
Nrep=0
for pH in 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 9.5 
do
Nrep=`expr $Nrep + 1`
dir=pH${pH}
echo $i $dir
mkdir $dir

cp md_pH_template.mdin ${dir}/pH_${pH}.mdin
perl -pi -e"s/PH/${pH}/" ${dir}/pH_${pH}.mdin

cat <<EOF >>groupfile
# pH ${pH}
-O -i ${dir}/pH_${pH}.mdin -p $parm -c $beginstr -cpin $cpin -o ${dir}/EXG_39D66H.pHREMD.mdout -cpout ${dir}/EXG_39D66H.pHREMD.cpout -cprestrt ${dir}/EXG_39D66H.pHREMD.cpin -r ${dir}/EXG_39D66H.pHREMD.rst7 -inf ${dir}/EXG_39D66H.pHREMD.mdinfo -x ${dir}/EXG_39D66H.pHREMD.nc
EOF
done


cat <<EOF >submit.sh
#!/bin/bash
#SBATCH --job-name=EXG_pHREMD
#SBATCH --nodes=1
#SBATCH --ntasks=$Ntasks
#SBATCH --gres=gpu:$Ngpu
#SBATCH --time=168:00:00
#SBATCH --partition=gpu

# Amber, OpenMPI, Cuda exports
module load GCC/5.4.0-2.26
module load CUDA/8.0.44
module load OpenMPI/2.0.1
module load Amber/16-AT-17-gpu

mpirun -np $Nrep \$AMBERHOME/bin/pmemd.cuda.MPI -ng $Nrep -groupfile groupfile -rem 4 -remlog pHremd.log

EOF

chmod u+x submit.sh
