# Copyright (C) 2022  Kristoffer E. Johansson  <kristoffer.johansson@bio.ku.dk>

options(width=200)

ko_strain = c("Doa10", "Ltn1", "San1", "Ubr1", "Ubr1Doa10", "Ubr1San1")
nstrains = length(ko_strain)
ngates = 4

psudo_counts = 0.5


###################################################################################################################
###      Raw reads and reads-per-million (RPM) normalization
###################################################################################################################
# Load reads
raw = read.delim("seqrun_12May20__raw_all.tab")

# Make unique rownames and remove oligo related columns
rn = paste0("oligo",sprintf("%04d",raw$oligo))
rownames(raw) = rn

count_col = paste0("strain.",rep(c(ko_strain,"wt"),each=ngates+1),"_gate.",rep(c(seq(ngates),"no"),times=1+length(ko_strain)))

# Function that normalize all columns to RPM
reads_per_million = function(df) {
    col_sums = apply(df, MARGIN=2, sum)
    df * matrix(10^6/col_sums, nrow=dim(df)[1], ncol=dim(df)[2], byrow=TRUE)
}

rpm = reads_per_million(raw[,count_col]+psudo_counts)
rownames(rpm) = rn


###################################################################################################################
###      Protein stability index (PSI)
###################################################################################################################
# Function to calculate PSI based on all columns of a data frame
protein_stability_index = function(df, name) {
    # First gate (gate 1) is stable with index 4, and the late (gate 4) unstable with index 1
    psi = t(apply(df, MARGIN=1, function(v){vsum=sum(v); c(vsum,sum( v/vsum * c(4,3,2,1) ))}))
    # put in a data frame and give column name
    ret_df = data.frame(psi)
    colnames(ret_df) = paste0(name,c(".rpm_sum",".psi"))
    return(ret_df)
}

# Calculate PSI for all strains
psi = protein_stability_index(rpm[,paste0("strain.wt_gate.",seq(ngates))], "wt")
for (strain in ko_strain) {
    coln = paste0("strain.", strain, "_gate.", seq(ngates))
    psi = cbind(psi, protein_stability_index(rpm[,coln], strain))
}

# Function to calculate delta PSI = PSI_KO - PSI_WT
delta_psi = function(df) {
    ret_df = df[,paste0(ko_strain,".psi")] - df[,"wt.psi"]
    colnames(ret_df) = paste0(ko_strain,".dpsi")
    return(ret_df)
}

# Calculate delta PSI
psi = cbind(psi, delta_psi(psi))


###################################################################################################################
###      PSI errors
###################################################################################################################
# Relative errors of reads, sqrt(n)/n = 1/sqrt(n)
rel_err = 1/sqrt(raw[,count_col]+psudo_counts)

# Function to calculate PSI error
psi_error = function(df, df_rel_err, name) {
    p = t(apply(df, MARGIN=1, function(v){ v/sum(v) }))
    p_err = p * df_rel_err
    psi_err = apply(p_err, MARGIN=1, function(v){sqrt(sum( (v * c(1,1,1,1) )**2 ))})
    # put in a data frame and give column name
    ret_df = data.frame(psi_err)
    colnames(ret_df) = paste0(name,".psi_err")
    return(ret_df)
}

# Calculate PSI error for all strains
for (strain in c("wt",ko_strain)) {
    cn_gates = paste0("strain.", strain, "_gate.", seq(ngates))
    cn = paste0(strain, ".psi_err")
    psi[,cn] = psi_error(rpm[,cn_gates], rel_err[,cn_gates], strain)
}

# Function to propagate PSI errors to delta PSI
delta_psi_err = function(df) {
    ret_df = sqrt(df[,paste0(ko_strain,".psi_err")]**2 + df[,"wt.psi_err"]**2)
    colnames(ret_df) = paste0(ko_strain,".dpsi_err")
    return(ret_df)
}

# Calculate delta PSI error
psi = cbind(psi, delta_psi_err(psi))


###################################################################################################################
###      Dump
###################################################################################################################

# Calculate sum of raw counts over agtes
for (strain in c(ko_strain,"wt")) {
    cn_gates = paste0("strain.", strain, "_gate.", seq(ngates))
    cn = paste0(strain, ".raw_sum")
    psi[,cn] = apply(raw[,cn_gates], MARGIN=1, sum)
}

# order columns
wt_cn = paste0("wt",c(".raw_sum",".rpm_sum",".psi",".psi_err"))
ko_cn = paste0(rep(ko_strain,each=6), c(".raw_sum",".rpm_sum",".psi",".psi_err",".dpsi",".dpsi_err"))
psi = psi[,c(wt_cn,ko_cn)]

save(psi, file="psi.rda")



