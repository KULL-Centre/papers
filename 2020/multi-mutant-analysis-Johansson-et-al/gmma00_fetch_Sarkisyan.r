# Copyright (C) 2017-2020 Kristoffer Enoe Johansson <kristoffer.johansson@bio.ku.dk>
# This file (gmma00_fetch_Sarkisyan.r) is part of the GMMA project

options(width=200, digits=4, stringsAsFactors=F)

file = "../../original_data/amino_acid_genotypes_to_brightness.tsv"
print("")
print("")
print(paste("Reading",file))
print("")
d = read.csv(file, sep="\t", header=F, skip=1)

colnames(d) = c("seq","N","signal","sd")

mut_list = strsplit(as.character(d[,1]),split=":")

# In Sarkisyan, all substitutions for some reason starts with an 'S'
mut_list = sapply(mut_list, substring, 2)

# ################################################################################
# # Introduce 2-letter target amino acid notation, one letter + single digit number
# ################################################################################
# mut_list = lapply(mut_list, paste, "1", sep="")

# obs_cut = 250
# t = table(unlist(mut_list))
# ti = which(t>(2*obs_cut))
# reparam_mut = names(t[ti])
# reparam_n = t[ti] %/% obs_cut

# names(reparam_n) = reparam_mut
# c = 0
# print(sprintf("Reparam %d substitutions to have >%d observations",length(reparam_mut),obs_cut))
# for (m in reparam_mut) {
#     c = c+1
#     mi = which(sapply(mut_list, function(l) {m %in% l}))
#     mut_list = lapply(mut_list, function(l) {if (length(l)==0) return(l); i=which(m==l); if (length(i)==1) l[i] = paste(substr(m,1,nchar(m)-1), sample(1:(reparam_n[m]),1), sep=""); return(l)})
#     print(sprintf("%3d/%3d: Reparam %s with %d observations into %d unique substitutions",c,length(reparam_n),m,t[m],reparam_n[m]))
#     # print(table(unlist(mut_list))[paste(substr(m,1,nchar(m)-1),seq(reparam_n[m]),sep="")])
#     # print(mut_list[mi])
# }


################################################################################
# Dump
################################################################################
subst = sapply(mut_list, paste, collapse=":")
subst[1] = "WT"
df = data.frame(variant=subst, bright=d$signal, n_syn=d$N, std=d$sd)

outfile="amino_acid_genotypes_to_brightness_parsed.tsv"
write.table(data.frame(subst=subst, N=d$N, signal=d$signal, sd=d$sd), file=outfile, row.names=F, col.names=F, quote=F)

print(sprintf("Dumped %s with %d substitutions and %d variants",outfile,length(unique(unlist(mut_list))),length(mut_list)))
