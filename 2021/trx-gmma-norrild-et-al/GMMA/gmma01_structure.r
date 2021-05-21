# Copyright (C) 2017-2021 Kristoffer Enoe Johansson <kristoffer.johansson@bio.ku.dk>
# This file (gmma01_structure.r) is part of the GMMA project

options(width=200, digits=4, stringsAsFactors=F)

args = commandArgs(trailingOnly=TRUE)
if (interactive()) {
    file = "oligo4p5_c30binary.txt"
} else if (length(args) < 1) {
    print("")
    print("usage: Rscript gmma01_structure.r <data_table.txt>")
    quit(save="no")
} else {
    file = args[1]
}
print(paste("Reading",file))
d = read.table(file)

# First row should be wild type
stopifnot(tolower(d[1,1]) %in% c("","wt","wildtype","ref","reference"))
d[1,1] = ""

# Wild-type average log brightness, number of WT measurements and standard deviation
colnames(d) = c("seq","N","signal","sd")
wt = d[1,] # no reason not to make wt a regular data point in mutant data.frame

print("Making data structures")

################################################################################
# Settings
################################################################################
settings = list()
settings$taa_letters = 1

################################################################################
# Per mutant data
################################################################################
# Make a list where elements are a vector of substitutions for that variant
mut_list = strsplit(as.character(d[,1]),split=":")

# Columns: aaMutations	uniqueBarcodes	medianBrightness	std
mutant = data.frame(i=seq_along(mut_list), signal=d[,'signal'], N_sub=sapply(mut_list, length), N_obs=d[,'N'], sd_obs=d[,'sd'])

rownames(mutant) = sapply(mutant$i, function(i) {sprintf("mut%05d",i)})
mutant$subst = sapply(mut_list, paste, collapse=":")

################################################################################
# Per substitution data
################################################################################
subst_table = table(unlist(mut_list))
nsubst = length(subst_table)

# Make data structure of substitutions
subst = data.frame(subst_table)
subst =	data.frame(i=0, resi=0, taa="", obs=subst$Freq, signal=NA, row.names=subst$Var1, stringsAsFactors=F)
for (si in seq(nsubst)) {
    m = rownames(subst)[si]
    subst[si,'resi'] = as.numeric(substr(m,2,nchar(m)-settings$taa_letters))
    subst[si,'taa'] = substr(m,nchar(m)-settings$taa_letters+1,nchar(m))
}

# Re-order according to residue number and add index
subst = subst[order(subst$resi),]
subst$i = seq(nsubst)

################################################################################
# Per residue data
################################################################################
# Copy of mut_list only containing residue numbers as integers
mut_int_list = lapply(mut_list, function(l) {as.integer(substr(l,2,nchar(l)-settings$taa_letters))})
# mut_int_list = sapply(mut_list, substring, 2)
# mut_int_list = sapply(mut_int_list, gsub, pattern='.$', replacement='')
# mut_int_list = sapply(mut_int_list, as.integer)

# number of residues from mutations
nres = max(unlist(mut_int_list))

# Get and check wild-type sequence from mutations
print("Check reference amino acid sequence")
rl = unlist(mut_int_list)
ml = unlist(mut_list)
wt_seq=""
for (resi in seq(nres)) {
    mask = rl == resi
    aa_list = substr(ml[mask], 1, 1)
    aa = unique(aa_list)
    if (length(aa) < 1) {
        print(sprintf("WARNING: Position %d has no data",resi))
	aa = 'x'
    } else if (length(aa) > 1) {
        print(sprintf("WARNING: Position %d has more than 1 identity",resi))
        t = sort(table(aa_list), decreasing=T)
	print(t)
        aa = t[1]
	print(paste("using",aa))
    } else {
        aa = aa[1]
    }
    wt_seq = paste(wt_seq, aa, sep='')
}
print(sprintf("Found wild-type sequence of length %d with %d missing residues:",nchar(wt_seq),sum("x"==strsplit(wt_seq,"")[[1]])))
wt[1] = wt_seq
print(wt_seq)

residue = data.frame(strsplit(wt_seq,""))
colnames(residue) ="wt"
residue$subst = ""
residue$N_mut = NA
# Assign substitutions
for (ri in seq(nres)) {
    residue[ri,'subst'] = paste(subst[subst$resi==ri,'taa'], collapse="")
    si = which(subst$resi==ri)
    residue[ri,"N_mut"] = sum(subst[si,"obs"])
}

################################################################################
# Build index translation lists between data frames
################################################################################
# Find substitutions for each position
print("Building res_mut_indices")
res_mut_indices = list()
for (i in seq(nres)) {
    # print(paste("position",i,"of",nres))
    mask = sapply(mut_int_list, function(array, val){val %in% array}, i)
    res_mut_indices[[i]] = which(mask)
}

# Find indices of substitutions for each measured mutant
print("Building mut_subst_indices")
mut_subst_indices = list()
for (muti in seq(length(mut_list))) {
    mut = mut_list[[muti]]
    mutl = length(mut)
    subst_indices = rep(0,mutl)
    if (mutl > 0) {
        for (mi in seq(mutl)) {
            si = which(rownames(subst)==mut[mi]) # There can be only one...
	    if (length(si) > 1) {
	        print(sprintf("WARNING: Mutatant %d (%s) with %d substitutions has %d matches (%s) for substitution %d (%s)",
	            muti,paste(mut,collapse=" "),mutl,length(si),paste(si,collapse=" "),mi,paste(rownames(subst)[si],collapse=" ")))
            }
            subst_indices[mi] = si
	}
    }
    mut_subst_indices[[muti]] = subst_indices
}
print("done building mut_subst_indices")

print("")
print("Building subst_mut_indices")
# For each single aa substitution in subst, the mut_list indices of mutants in which the substitution occurs
subst_mut_indices = list()
for (si in seq(nsubst)) {
    # print(sprintf("Substitution %4d of %d",si,nsubst))
    m = rownames(subst)[si]
    subst_mut_indices[[si]] = which(sapply(mut_list, function(mut) {m %in% mut}))
}
print("done building subst_mut_indices")

################################################################################
# Assign single mutant signal
################################################################################
counter = 0
for (si in seq(nsubst)) {
    mi = subst_mut_indices[[si]]
    mii = which(sapply(mut_list[mi],length)==1)
    stopifnot(length(mii) < 2)
    if (length(mii) > 0) {
        subst[si,"signal"] = mutant[mi[mii[1]],"signal"]
	counter = counter + 1
    }
}
print(sprintf("%d of %d substitutions are oberserved as single-mutants",counter,dim(subst)[1]))

################################################################################
# Dump
################################################################################
subst$note = ""

# The data in structured format
# version=2 is necessary for R versions <3.5 to be able to read it (binf version 3.1)
save(settings, wt, mut_list, mutant, subst, residue, mut_subst_indices, subst_mut_indices, res_mut_indices, file="gmma_structured.rda", version=2)
print("Saved gmma_structured.rda")

print(sprintf("Successfully generated indexed data structures of %d variants carrying %d unique substitutions",length(mut_list),dim(subst)[1]))
