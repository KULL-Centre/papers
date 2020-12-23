#for lambda in {1.00,1.01,1.02,1.03,1.04,1.05,1.06,1.07,1.08,1.09,1.10,1.20,1.50}
for c in 50 150 250 400 1000
do
    folder=${c}mM/stage
    program=run_pep_stage2_v4.py
    cp $program $folder
    cd $folder
    rm nohup.out 
    nohup /storage1/thomasen/software/miniconda3/bin/python3 $program FLA1_${c}mMNaCl &
    cd ../../
done
