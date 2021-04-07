seqs=$1 
path=$(echo $seqs | rev | cut -d"/" -f2-  | rev)'/'
basen=$(echo $(basename $seqs) | cut -f1 -d"_")
title=$(echo $(basename $seqs) | cut -f1 -d".")
output=$path$basen
#output=$(echo $(basename $seqs) | cut -f1,2 -d"_")
weblogo -f $seqs -o $output'_seq_logo.pdf' -F pdf -U 'probability' -n 52 -t $title -c chemistry -s large
