#!/bin/bash

#  Input files : 
#  $1 : plmDCA statistal model of the MSA [J,h] in .mdl format ! OR ! MSA in FASTA format (Only UPPER CASES and - considered)
#     ! OR ! MSA in FASTA format (Only UPPER CASES and - considered, each sequence on a single line)
#  $2 : sequence(s) to score and to compute full mutagenesis from

#  Output files :
#  ON SCREEN : energies for sequences in $2
#  mut_screen_$1_$2 : full mutational screening for sequences in $2

#  EXAMPLES: 
#  a) ./run_saturation_mutagenesis.sh your_MSA.fa your_seq.fa       (to create also the .mdl)
#  b) ./run_saturation_mutagenesis.sh mdl_MSA.mdl your_seq.fa

i1=$1
i2=$2


if [ "${i1}" == "-h" ]  || [ $# -le 1 ] ; then
   echo ""
   echo "#  HELP MESSAGE for run_saturation_mutagenesis.sh : "
   echo ""

   echo "#  It compute the energy-score and runs full mutagenesis for the sequence(s) in input using"
   echo "#  a statistical model built by plmDCA from a MSA given in input"   

   head -n 15  /storage1/shared/software/SEQ_MODEL_DESIGN/run_saturation_mutagenesis.sh | tail -n 14

   exit
fi


i1_name=`echo $1 | awk 'BEGIN{FS="."}{print $1}'`
i1_ext=`echo $i1 | awk 'BEGIN{FS="."}{print $2}'` 

if [ ${i1_ext} != "mdl" ] ; then
   echo "" 
   
   i1_name=`echo $1 | awk 'BEGIN{FS="."}{print "mdl_"$1}'`

   if [ -f ${i1_name}.mdl ] ; then 
      echo "#  WARNING: statistical model for this MSA already exists. "
      echo "#           CHECK the input files, delete ${i1_name}.mdl or rename the MSA "
      echo ""
      exit
   fi      
   echo "" 
   echo "#  COMPUTING THE STATISTICAL MODEL WITH create_PLMDCA_model_from_MSA.sh " 
 
   /storage1/shared/software/SEQ_MODEL_DESIGN/create_PLMDCA_model_from_MSA.sh ${i1}   
   
   echo "#  DONE "
fi

echo "" 
echo "#  RUNNING SATURATION MUTAGENESIS ... "
echo "" 

./score_or_min_seqs.x -SEQ ${i2} -PAR ${i1_name}.mdl -ENE 

 


