#!/bin/bash 
#SBATCH --job-name=10_AVG
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=2000
#SBATCH -t 72:00:00
#SBATCH -o opt.out
#SBATCH -e opt.err
#SBATCH --partition=sbinlab_ib

source /groups/sbinlab/giulio/.bashrc
conda activate hoomd

declare -a proteins_list=(OPN FUS FUS12E Sic1 aSyn A2 A1 Hst5 aSyn140 PNt Hst52 ACTR RNaseA p15PAF OPN220 Sic92 FhuA CoRNID ColNT hNL3cyt SH4UD K10 K27 K25 K32 K23 K44 M12FP12Y P7FM7Y M9FP6Y M8FP4Y M9FP3Y M10R M6R P2R P7R M3RP3K M6RP6K M10RP10K M4D P4D P8D P12D P12E P7KP12D P7KP12Db M12FP12YM10R M10FP7RP12D)
declare -a proteinsPRE_list=(OPN FUS FUS12E Sic1 aSyn A2)

# Clean-up the folder
rm log
for name in ${proteins_list[@]}
do
rm -r $name
mkdir $name
done

# Copy the experimental NMR PRE data in the respective folders
for name in ${proteinsPRE_list[@]}
do
cp -r ../expPREs/$name/expPREs $name
done

cp *.py $SCRATCH
cp *.pkl $SCRATCH

for name in ${proteins_list[@]}
do
cp -r $name $SCRATCH
echo $name
done

echo $SLURM_CPUS_PER_TASK

cd $SCRATCH

start=$(date +%s.%N)

echo $SLURM_CPUS_ON_NODE

# Perform 5 consecutive simulation-reweighting/simulated annealing cycles                                                                                                                    
for cycle in {o1,o2,o3,o4,o5}
do

for name in ${proteins_list[@]}
do
mkdir $name/$cycle
done

# Run all the simulations in parallel 
python ./simulate.py --outdir $cycle --num_cpus $SLURM_CPUS_PER_TASK

cp -r * $SLURM_SUBMIT_DIR

# Covert to all-atom trajectories (only for predicitons of NMR PRE data)
for name in ${proteinsPRE_list[@]}
do
python ./pulchra.py --name $name --num_cpus $SLURM_CPUS_PER_TASK --pulchra /groups/sbinlab/giulio/pulchra_306/pulchra
sleep 20
done

# Run the reweighting/simulated annealing routine
python ./optimize.py --log $SLURM_SUBMIT_DIR --cycle $cycle --num_cpus $SLURM_CPUS_PER_TASK

cp -r * $SLURM_SUBMIT_DIR

done

duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`

echo $execution_time
