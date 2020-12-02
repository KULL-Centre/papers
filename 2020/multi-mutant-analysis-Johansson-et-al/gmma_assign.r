# Copyright (C) 2017-2020 Kristoffer Enoe Johansson <kristoffer.johansson@bio.ku.dk>
# This file (gmma_assign.r) is part of the GMMA project

options(width=200, digits=4, stringsAsFactors=F)

args = commandArgs(trailingOnly=TRUE)
if (interactive()) {
    file = "gmma_structured.rda"
} else if (length(args) < 1) {
    print("")
    print("usage: Rscript gmma_assign.r <*.rda>")
    quit(save="no")
} else {
    file = args[1]
}
print(sprintf("Read %s",file))
loaded_obj = load(file)

stopifnot( all(c("subst","mutant","residue") %in% loaded_obj) )

nmut = length(mutant[,1])
nsubst = length(subst[,1])
nres = length(residue[,1])

# rownames of the subst data.frame in single-letter target amino acid
subst_taa_one = substr(rownames(subst),1,nchar(rownames(subst))-settings$taa_letters+1)


################################################################################
# Assign secondary structure and burial to residues
################################################################################
# residue$ss= rep(NA,nres)
residue$hse = NA
residue$asa = NA
residue$burial2 = NA
residue$burial3 = NA

# Add exposure
expo = read.table("assignments/2b3p_pre_expo.dat", comment.char="#", stringsAsFactors=F)
colnames(expo) = c("Res","AA","HSE","ASA","3c","2c")
n2one = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
for (ri in seq(nres)) {
    if (ri+1 > length(expo$AA)) {
        residue[ri,'burial2'] = 'E'
        residue[ri,'burial3'] = 'E'
	next
    }
    expo_aa = n2one[expo[ri+1,'AA']+1]
    if (as.character(residue[ri,'wt']) != expo_aa) { print(sprintf("Sarkisyan res %d is %s but expo aa is %s",ri,residue[ri,'wt'],expo_aa)) }
    residue[ri,'hse'] = expo[ri+1,'HSE']
    residue[ri,'asa'] = expo[ri+1,'ASA']
    residue[ri,'burial2'] = expo[ri+1,'2c']
    residue[ri,'burial3'] = expo[ri+1,'3c']
}

################################################################################
# Known GFP mutatuions
################################################################################
residue$av = "-"
residue$sf = "-"
residue$ts = "-"
residue$do = "-"
residue$now = "-"
residue$pross = "-"
residue$pdb = "-"

# #Reset the note column if used to reassign. Only keep 'ifs' and 'nonsense' in subst$note
# subst$note = sapply(subst$note, function(s) {sl = strsplit(s," ")[[1]]; sr=""; if ("ifs" %in% sl) sr="ifs"; if ("nonsense" %in% sl) sr=paste("nonsense",sr,collapse=" "); return(sr) })

# Aequorea victoria uniprot P42212, 238 residues, removed NT MS(KGEE...)
# PDB: 1GFL 1.9Å L64F+Q80R immature, 1EMA 1.9Å L64F+S65T+Q80R mature. 
avGFP_seq = paste("KGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLV",
                  "NRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK", sep="")
# The early mutants
# Heim95         S65T
# Cormack96 mut1 F64L + S65T (codon optimized and called eGFP by Yang96)
#           mut2 S65A + V68L + S72A
#           mut3 S65G + S72A

# Pedelacq06: eGFP (F64L+S65T) and cycle-3 (F99S+M153T+V163A) additionally S30R+Y39N+N105T+Y145F+I171V+A206V.
# Mutations are tested indivitually, main ref of PDB 2B3P (also FR-GFP 2B3Q) which also includes Q80R?!
# SEQADV 2B3P ARG A   80  UNP  P42212    GLN    80 ENGINEERED
# Chudakov10 says that the improved solubility with fusions commes from dimerization facilitated by A206V
# Sequence from 2B3P, removed NT MS(KGEE...) and CT (...ELK)GSHHHHHH.
sfGFP_seq  = paste("KGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELK",
                   "GIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK", sep="")
# T-Sapphire, seq from fpbase.org [Zapata-Hommer03], cut MVS
tsGFP_seq = paste("KGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVMVFARYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLV",
                  "NRIELKGIDFKEDGNILGHKLEYNFNSHNVYIMADKQKNGIKANFKIRHNIEDGGVQLADHYQQNTPIGDGPVLLPDNHYLSIQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK", sep="")
# Sarkisyan15
nowGFP_seq = paste("KGEKLFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKMSLKFICTTGKLPVPWPTLKTTLTWGMQCFARYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELK",
                   "GVDFKEDGNILGHKLEYNAISGNANITADKQKNGIKAYFTIRHDVEDGSVLLADHYQQNTPIGDGPVLLPDNHYLSTQSKQSKDPNEKRDHMVLLEFVTAAGIPLGADELYK", sep="")

# Note this is a circular mutant that I have copy-pasted back to match the above
doGFP_seq  = paste("QGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATIGKLTLKFISTTGKLPVPWPTLVTTLSYGVQAFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGKYKTRAVVKFEGDTLVNRIELK",
	       	   "GTDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFTVRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQTVLSKDPNEKRDHMVLHEYVNAAGITHGMDELYG", sep="")
# Bandyopadhyay17 des11 designed using PROSS, plus S65T and H231L. 
# pross_seq  = paste("KGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGGQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGCYKTRAEVKFEGDTLVNRIELH", #this is des12 containing V66G
#                    "GIDFKEDGNILGHKLEYNFNSHNVYIMPDKQNNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDYHYLHTWSELSKDPNEKRDHMVLLEFVEAAGITLGMDELYK", sep="")
pross_seq  = paste("KGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGCYKTRAEVKFEGDTLVNRIELH",
                   "GIDFKEDGNILGHKLEYNFNSHNVYIMPDKQNNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDYHYLHTWSELSKDPNEKRDHMVLLEFVEAAGITLGMDELYK", sep="")

# Bandyopadhyay17 GroE increased dependence relative to eGFP 
GroE_dep   = list(c("N",23,"K"), c("D",234,"H"), c("N",144,"I"), c("P",211,"L"), c("N",212,"K"), c("T",97,"I"), c("F",114,"I"), c("L",178,"R"), c("R",215,"C"))
# Decreased dependence
GroE_indep = list(c("K",162,"N"), c("K",166,"M"), c("Y",39,"N"), c("P",58,"S"), c("S",30,"R"))

GroE_des11 = list(c("S",30,"R"), c("N",105,"C"), c("K",126,"H"), c("Y",145,"F"), c("A",154,"P"), c("K",158,"N"), c("N",198,"Y"), c("S",202,"H"), c("Q",204,"W"), c("A",206,"E"), c("T",225,"E"))
# my wt:  KGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTxxNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELY
# 2wur:   KGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKTRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHG
# Pross: SKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGGQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGCYKTRAEVKFEGDTLVNRIELHGIDFKEDGNILGHKLEYNFNSHNVYIMPDKQNNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDYHYLHTWSELSKDPNEKRDHMVLLEFVEAAGITLGMDELYK
# So, from my_wt to 2wur is Q78R, I165T but both are reverted in the pross sequence
# Pross not reported: S63T (considered eGFP mut but not in Saarkisyan WT), also H229L but L229 is sometimes considered av/eGFP. Uniprot P42212 and Sarkisyan have H229.
# K206E reported in paper is a typo - the translated DNA says A

# Question is: Does decreased GroE dependence correlate with the number of 'general mutations' that eGFP can tolerate?
# Answer is not clear: K166M seems really good from Sarkisyan data. K162N and P58S may be good but there's not enough data.
# The sfGFP mutations Y39N and Y145F are good from all perspectives.
# I got N198E from Rosetta and Jakob found N198T as one of 'the 5' - from PROSS it's N198Y. Position plot of 198 shows that most are better than Asn.
# I get S202T from Rosetta - from PROSS it's S202H
# eGFP has K206 but the crystal structure A206? Anyway, sfGFP has A206T and PROSS finds K206E. Position plot shows that most are better than A (K? Sarkisyan wt)

# for (m in GroE_des11) {
#     resi = as.numeric(m[2])-2
#     stopifnot(m[1] == residue[resi,"wt"])
#     mut=paste(m[1],resi,m[3], sep="")
#     residue[resi,'pross'] = m[3]
#     if (mut %in% rownames(subst)) {
#         subst[mut,'note'] = paste(subst[mut,'note'],'pross')
# 	print(subst[mut,])
#     } else {
#         print(sprintf("Pross mutation %s is not in Sarkisyan data",mut))
#     }
# }
si = which("V66G"==subst_taa_one)
for (i in si) subst[i,'assign'] = paste(subst[i,'assign'],'pross_n')

pdb_seq = read.table("assignments/pdb_trim.seq", stringsAsFactors=F)
for (resi in seq(nres)) {
    # Aequorea victoria sequence
    av_aa = substr(avGFP_seq, resi, resi)
    if (av_aa != residue[resi,'wt']) residue[resi,'av'] = av_aa

    # substitutions published in the PDB
    pdb_aa = c()
    for (iseq in seq(length(pdb_seq[,1]))) {
        aa = substr(pdb_seq[iseq,1],resi,resi)
        if (aa != residue[resi,'wt'] && aa != '-' && aa != 'X') pdb_aa = append(pdb_aa, aa)
    }
    if (length(pdb_aa) > 0) residue[resi,'pdb'] = paste(unique(pdb_aa), collapse="")

    # substiotutions from sfGFP
    sf_aa = substr(sfGFP_seq, resi, resi)
    if (sf_aa != residue[resi,'wt']) residue[resi,'sf'] = sf_aa

    # substiotutions from tsGFP
    ts_aa = substr(tsGFP_seq, resi, resi)
    if (ts_aa != residue[resi,'wt']) residue[resi,'ts'] = ts_aa

    # substiotutions from nowGFP
    now_aa = substr(nowGFP_seq, resi, resi)
    if (now_aa != residue[resi,'wt']) residue[resi,'now'] = now_aa

    # substiotutions from Do & Boxer split GFP
    do_aa = substr(doGFP_seq, resi, resi)
    if (do_aa != residue[resi,'wt']) residue[resi,'do'] = do_aa

    pross_aa = substr(pross_seq, resi, resi)
    if (pross_aa != residue[resi,'wt']) residue[resi,'pross'] = pross_aa
}

# Put a note in subst if the mutation is known
sf = c()
ts = c()
do = c()
now = c()
pdb = c()
pross = c()
sf_miss = 0
ts_miss = 0
do_miss = 0
now_miss = 0
pdb_miss = 0
pross_miss = 0
n_pdb = 0
subst$assign = ""
for (resi in seq(nres)) {
    for (taa in strsplit(residue$pdb[resi],"")[[1]]) {
        mut = paste(residue[resi,'wt'],resi,taa,sep="")
	n_pdb = n_pdb+1
	if (mut %in% subst_taa_one) {
	     pdb = append(pdb, mut)
	     si = which(mut==subst_taa_one)
	     # subst[mut,'assign'] = paste(subst[mut,'assign'], 'pdb')
	     si = which(mut==subst_taa_one)
	     for (i in si) subst[i,'assign'] = paste(subst[i,'assign'], 'pdb')
	} else {
	    # if (residue[resi,'wt']!='x' && taa!='-') { print(mut); pdb_miss = pdb_miss+1 }
	    if (residue[resi,'wt']!='x' && taa!='-') { pdb_miss = pdb_miss+1 }
        }
    }    
    mut = paste(residue[resi,'wt'],resi,residue[resi,'sf'],sep="")
    if (mut %in% subst_taa_one) {
        sf = append(sf, mut)
        si = which(mut==subst_taa_one)
	for (i in si) subst[i,'assign'] = paste(subst[i,'assign'], 'sf')
    } else {
        # if (residue[resi,'wt']!='x' && residue[resi,'sf']!='-') { print(mut); sf_miss = sf_miss+1 }
        if (residue[resi,'wt']!='x' && residue[resi,'sf']!='-') { sf_miss = sf_miss+1 }
    }
    
    mut = paste(residue[resi,'wt'],resi,residue[resi,'ts'],sep="")
    if (mut %in% subst_taa_one) {
        ts = append(ts, mut)
        si = which(mut==subst_taa_one)
	for (i in si) subst[i,'assign'] = paste(subst[i,'assign'], 'ts')
    } else {
        if (residue[resi,'wt']!='x' && residue[resi,'ts']!='-') { ts_miss = ts_miss+1 }
    }
    
    mut = paste(residue[resi,'wt'],resi,residue[resi,'do'],sep="")
    if (mut %in% subst_taa_one) {
        do = append(do, mut)
        si = which(mut==subst_taa_one)
	for (i in si) subst[i,'assign'] = paste(subst[i,'assign'], 'do')
    } else {
        # if (residue[resi,'wt']!='x' && residue[resi,'do']!='-') { print(mut); do_miss = do_miss+1 }
        if (residue[resi,'wt']!='x' && residue[resi,'do']!='-') { do_miss = do_miss+1 }
    }
    
    mut = paste(residue[resi,'wt'],resi,residue[resi,'now'],sep="")
    if (mut %in% subst_taa_one) {
        now = append(now, mut)
        si = which(mut==subst_taa_one)
	for (i in si) subst[i,'assign'] = paste(subst[i,'assign'], 'now')
    } else {
        # if (residue[resi,'wt']!='x' && residue[resi,'now']!='-') { print(mut); now_miss = now_miss+1 }
        if (residue[resi,'wt']!='x' && residue[resi,'now']!='-') { now_miss = now_miss+1 }
    }    
    
    mut = paste(residue[resi,'wt'],resi,residue[resi,'pross'],sep="")
    if (mut %in% subst_taa_one) {
        pross = append(pross, mut)
        si = which(mut==subst_taa_one)
	for (i in si) subst[i,'assign'] = paste(subst[i,'assign'], 'pross')
    } else {
        # if (residue[resi,'wt']!='x' && residue[resi,'now']!='-') { print(mut); now_miss = now_miss+1 }
        if (residue[resi,'wt']!='x' && residue[resi,'pross']!='-') { pross_miss = pross_miss+1 }
    }    
}
print(sprintf("Missed mutations: sf=%d/%d, ts=%d/%d, do=%d/%d, now=%d/%d, pross=%d/%d, pdb=%d/%d",sf_miss,sum(residue$sf!='-'),ts_miss,sum(residue$ts!='-'),do_miss,sum(residue$do!='-'),now_miss,sum(residue$now!='-'),pross_miss,sum(residue$pross!='-'),pdb_miss,n_pdb))

for (i in seq(nres)) {
    # for (tag in c("av","sf","do","now","pdb")) {
    for (tag in c("av","sf","ts","do","now","pross")) {
        if (residue[i,tag] != '-') {
            mut = paste(residue[i,"wt"],i,residue[i,tag], sep="")
            if (! mut %in% subst_taa_one) { print(sprintf("Known %4s mutant %5s is not in data",tag,mut)) } # else { print(sprintf("Found %s %s",tag,mut)) }
        }
    }
}

file_id = strsplit(file,".rda")[[1]]
outfile = paste(file_id, "assigned.rda", sep="_")
# version=2 is necessary for R versions <3.5 to be able to read it (binf version 3.1)
save(list=loaded_obj, file=outfile, version=2)
print(sprintf("Saved %s containing: %s",outfile,paste(loaded_obj,collapse=", ")))
