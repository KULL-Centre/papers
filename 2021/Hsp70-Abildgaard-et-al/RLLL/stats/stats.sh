#!/bin/bash

# collects stats

indir=../ref
inserts=inserts_with_X_for_stop.txt
scoredir=../binders

for insert in `cat ${indir}/${inserts}  `
do
  name=${insert}

  limbodata=${scoredir}/${name}.limbo_scores
  scores=`awk -f stats.awk < ${limbodata}`
  echo ${name} ${scores}
done

