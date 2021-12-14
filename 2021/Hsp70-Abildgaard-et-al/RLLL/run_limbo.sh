#!/bin/bash

p=/usr/bin/python2.7
limbo=/Users/lindorff/program/python/limbo/limbo/score.py
matrix=/Users/lindorff/program/python/limbo/limbo/mergedmatrix.mat

indir=ref
seqdir=sequences
inserts=inserts_with_X_for_stop.txt
NTER=CALLQS
CTER=SAPRRAAATARY

for insert in `cat ${indir}/${inserts}  `
do
  name=${insert}

  fasta=${seqdir}/${name}.fasta
  echo "> "${insert} > ${fasta}
  echo ${NTER}${insert}${CTER} >> ${fasta}
  sed -i -e 's/X.*$//g' ${fasta}
  echo "Running LIMBO on" ${name}

  #fasta=${name}.fasta
  binders=${name}.fasta_binders.txt
  boltzmann=${name}.fasta_boltzmann.txt
  binders_output=${name}.limbo_scores
  boltzmann_output=${name}.limbo_boltzmann


  ${p} ${limbo} ${matrix} ${fasta}
  awk '{print($2,$3)}' ${seqdir}/${binders} > binders/${binders_output}
  awk '{print($3,$4)}' ${seqdir}/${boltzmann} > boltzmann/${boltzmann_output}
  rm ${seqdir}/${binders}
  rm ${seqdir}/${boltzmann}

done

