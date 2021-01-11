#can be run as a different script only upon need, so it doesn't overcomplicate the main script for data cleaning and analysis (ArrayIIIb_Clean and plot data.Rmd)

#Creating a gathered substitution table with standard deviations instead of average fluorescence, to be used in Fig. S4

d4 <- d[,c(2,4)]
#The loop didn't work out 
#wt_sequence <- c("LPDNHYLSTQTVLSKDPN", "LPDNYYLSTQTVLSKDPN", "LPDNHYLSYQTVLSKDPN","LPDNHYLSIQTVLSKDPN","LPDNHYLSTQTVVSKDPN","LPDNHYLSTQTVRSKDPN","LPDNYYLSYQTVLSKDPN")
#context <- c("WT","H199Y","T203Y","T203I","L207V","L207R","H199Y/T203Y")

first_position_in_wt <- 198
substituting_aa <- "RHKDECGPAVILMFYWSTNQ"
#######################################################################WT
wt_sequence <- "LPDNHYLSTQTVLSKDPN"
context <- "WT"
# Extract sequence of the WT 
wtseq <- vector()
for(i in 1:nchar(as.character(wt_sequence))) {
  wtseq[i]<-substr(wt_sequence,i,i) 
}
# Subset dataframe to only include peptides of the same length as the wt
sameLengthAsWT <-subset(d4,nchar(as.character(d4$seq),type="chars") == nchar(as.character(wt_sequence),type="chars"))
# Create vector with single amino acid substitutions only
Boolean_string <- c()
for(m in 1:length(sameLengthAsWT$seq)) {  #for all peptides of same length as wt
  v <- unlist(strsplit(as.character(wt_sequence),character())) #wt sequence
  w <- unlist(strsplit(as.character(sameLengthAsWT$seq[m]),character()))  #sequence of the other peptides
  if(sum(v==w) == nchar(as.character(wt_sequence))-1) {  #if wt sequence and sequence of other peptide has all but one amino acids in common...
    Boolean_string[m]<-m        #returning a vector of the positions in sameLengthAsWT$Peptide meeting the criteria
  }
}
Boolean_string<-na.omit(Boolean_string)
SingleSubsW<-sameLengthAsWT$seq[Boolean_string]
# Creating a data frame with wt amino acids as column names and turning it into a vector
new_list <- list()
for(j in first_position_in_wt:(first_position_in_wt + nchar(as.character(wt_sequence)) - 1)) {
  new_list[[paste0(wtseq[j-first_position_in_wt+1],j,sep="")]] = NA
}
wt_amino_acids<-as.data.frame(new_list)
wt_amino_acids<-as.vector(names(wt_amino_acids))
# Vector with substituting amino acids
substituting_amino_acids<-unlist(strsplit(substituting_aa,character()))
# Create empty matrix with the substitution table format
no_substitutions_sd <- d4[d4$seq==wt_sequence,"sd"]
substitution_table <- matrix(data=c(no_substitutions_sd), nrow = length(substituting_amino_acids),ncol = length(wtseq))
colnames(substitution_table) <- wt_amino_acids
rownames(substitution_table) <- substituting_amino_acids
# Fill table with sd values
for(i in 1:length(SingleSubsW)) {
  for(j in 1:nchar(wt_sequence)) {
    if(substr(SingleSubsW[i],j,j) != wtseq[j]) {
      substitution_table[substr(SingleSubsW[i],j,j),j] <- d4[d4$seq==SingleSubsW[i],"sd"]
    }
  }
}
# Rearrange substitution_table for heatmap
st <- data.frame(substitution_table)
library(tidyr)
st_gathered <- gather(st)
names(st_gathered) <- c("Position", "SD") 
# Keep only the position number without the WT residue
st_gathered$Position <- substring(st_gathered$Position, 2)
# Add the substituting AA and context and rearrange the DF
st_gathered$Substitution <- substituting_amino_acids
st_gathered$Context <- context
st_gathered_WT <- st_gathered[,c(1,3,2,4)]


#######################################################################H199Y
wt_sequence <- "LPDNYYLSTQTVLSKDPN"
context <- "H199Y"
# Extract sequence of the WT 
wtseq <- vector()
for(i in 1:nchar(as.character(wt_sequence))) {
  wtseq[i]<-substr(wt_sequence,i,i) 
}
# Subset dataframe to only include peptides of the same length as the wt
sameLengthAsWT <-subset(d4,nchar(as.character(d4$seq),type="chars") == nchar(as.character(wt_sequence),type="chars"))
# Create vector with single amino acid substitutions only
Boolean_string <- c()
for(m in 1:length(sameLengthAsWT$seq)) {  #for all peptides of same length as wt
  v <- unlist(strsplit(as.character(wt_sequence),character())) #wt sequence
  w <- unlist(strsplit(as.character(sameLengthAsWT$seq[m]),character()))  #sequence of the other peptides
  if(sum(v==w) == nchar(as.character(wt_sequence))-1) {  #if wt sequence and sequence of other peptide has all but one amino acids in common...
    Boolean_string[m]<-m        #returning a vector of the positions in sameLengthAsWT$Peptide meeting the criteria
  }
}
Boolean_string<-na.omit(Boolean_string)
SingleSubsh<-sameLengthAsWT$seq[Boolean_string]
# Creating a data frame with wt amino acids as column names and turning it into a vector
new_list <- list()
for(j in first_position_in_wt:(first_position_in_wt + nchar(as.character(wt_sequence)) - 1)) {
  new_list[[paste0(wtseq[j-first_position_in_wt+1],j,sep="")]] = NA
}
wt_amino_acids<-as.data.frame(new_list)
wt_amino_acids<-as.vector(names(wt_amino_acids))
# Vector with substituting amino acids
substituting_amino_acids<-unlist(strsplit(substituting_aa,character()))
# Create empty matrix with the substitution table format
no_substitutions_sd <- d4[d4$seq==wt_sequence,"sd"]
substitution_table <- matrix(data=c(no_substitutions_sd), nrow = length(substituting_amino_acids),ncol = length(wtseq))
colnames(substitution_table) <- wt_amino_acids
rownames(substitution_table) <- substituting_amino_acids
# Fill table with sd values
for(i in 1:length(SingleSubsh)) {
  for(j in 1:nchar(wt_sequence)) {
    if(substr(SingleSubsh[i],j,j) != wtseq[j]) {
      substitution_table[substr(SingleSubsh[i],j,j),j] <- d4[d4$seq==SingleSubsh[i],"sd"]
    }
  }
}
# Rearrange substitution_table for heatmap
st <- data.frame(substitution_table)
library(tidyr)
st_gathered <- gather(st)
names(st_gathered) <- c("Position", "SD") 
# Keep only the position number without the WT residue
st_gathered$Position <- substring(st_gathered$Position, 2)
# Add the substituting AA and context and rearrange the DF
st_gathered$Substitution <- substituting_amino_acids
st_gathered$Context <- context
st_gathered_H199Y <- st_gathered[,c(1,3,2,4)]

#######################################################################T203Y
wt_sequence <- "LPDNHYLSYQTVLSKDPN"
context <- "T203Y"
# Extract sequence of the WT 
wtseq <- vector()
for(i in 1:nchar(as.character(wt_sequence))) {
  wtseq[i]<-substr(wt_sequence,i,i) 
}
# Subset dataframe to only include peptides of the same length as the wt
sameLengthAsWT <-subset(d4,nchar(as.character(d4$seq),type="chars") == nchar(as.character(wt_sequence),type="chars"))
# Create vector with single amino acid substitutions only
Boolean_string <- c()
for(m in 1:length(sameLengthAsWT$seq)) {  #for all peptides of same length as wt
  v <- unlist(strsplit(as.character(wt_sequence),character())) #wt sequence
  w <- unlist(strsplit(as.character(sameLengthAsWT$seq[m]),character()))  #sequence of the other peptides
  if(sum(v==w) == nchar(as.character(wt_sequence))-1) {  #if wt sequence and sequence of other peptide has all but one amino acids in common...
    Boolean_string[m]<-m        #returning a vector of the positions in sameLengthAsWT$Peptide meeting the criteria
  }
}
Boolean_string<-na.omit(Boolean_string)
SingleSubsty<-sameLengthAsWT$seq[Boolean_string]
# Creating a data frame with wt amino acids as column names and turning it into a vector
new_list <- list()
for(j in first_position_in_wt:(first_position_in_wt + nchar(as.character(wt_sequence)) - 1)) {
  new_list[[paste0(wtseq[j-first_position_in_wt+1],j,sep="")]] = NA
}
wt_amino_acids<-as.data.frame(new_list)
wt_amino_acids<-as.vector(names(wt_amino_acids))
# Vector with substituting amino acids
substituting_amino_acids<-unlist(strsplit(substituting_aa,character()))
# Create empty matrix with the substitution table format
no_substitutions_sd <- d4[d4$seq==wt_sequence,"sd"]
substitution_table <- matrix(data=c(no_substitutions_sd), nrow = length(substituting_amino_acids),ncol = length(wtseq))
colnames(substitution_table) <- wt_amino_acids
rownames(substitution_table) <- substituting_amino_acids
# Fill table with sd values
for(i in 1:length(SingleSubsty)) {
  for(j in 1:nchar(wt_sequence)) {
    if(substr(SingleSubsty[i],j,j) != wtseq[j]) {
      substitution_table[substr(SingleSubsty[i],j,j),j] <- d4[d4$seq==SingleSubsty[i],"sd"]
    }
  }
}
# Rearrange substitution_table for heatmap
st <- data.frame(substitution_table)
library(tidyr)
st_gathered <- gather(st)
names(st_gathered) <- c("Position", "SD") 
# Keep only the position number without the WT residue
st_gathered$Position <- substring(st_gathered$Position, 2)
# Add the substituting AA and context and rearrange the DF
st_gathered$Substitution <- substituting_amino_acids
st_gathered$Context <- context
st_gathered_T203Y <- st_gathered[,c(1,3,2,4)]

#######################################################################T203I
wt_sequence <- "LPDNHYLSIQTVLSKDPN"
context <- "T203I"
# Extract sequence of the WT 
wtseq <- vector()
for(i in 1:nchar(as.character(wt_sequence))) {
  wtseq[i]<-substr(wt_sequence,i,i) 
}
# Subset dataframe to only include peptides of the same length as the wt
sameLengthAsWT <-subset(d4,nchar(as.character(d4$seq),type="chars") == nchar(as.character(wt_sequence),type="chars"))
# Create vector with single amino acid substitutions only
Boolean_string <- c()
for(m in 1:length(sameLengthAsWT$seq)) {  #for all peptides of same length as wt
  v <- unlist(strsplit(as.character(wt_sequence),character())) #wt sequence
  w <- unlist(strsplit(as.character(sameLengthAsWT$seq[m]),character()))  #sequence of the other peptides
  if(sum(v==w) == nchar(as.character(wt_sequence))-1) {  #if wt sequence and sequence of other peptide has all but one amino acids in common...
    Boolean_string[m]<-m        #returning a vector of the positions in sameLengthAsWT$Peptide meeting the criteria
  }
}
Boolean_string<-na.omit(Boolean_string)
SingleSubsti<-sameLengthAsWT$seq[Boolean_string]
# Creating a data frame with wt amino acids as column names and turning it into a vector
new_list <- list()
for(j in first_position_in_wt:(first_position_in_wt + nchar(as.character(wt_sequence)) - 1)) {
  new_list[[paste0(wtseq[j-first_position_in_wt+1],j,sep="")]] = NA
}
wt_amino_acids<-as.data.frame(new_list)
wt_amino_acids<-as.vector(names(wt_amino_acids))
# Vector with substituting amino acids
substituting_amino_acids<-unlist(strsplit(substituting_aa,character()))
# Create empty matrix with the substitution table format
no_substitutions_sd <- d4[d4$seq==wt_sequence,"sd"]
substitution_table <- matrix(data=c(no_substitutions_sd), nrow = length(substituting_amino_acids),ncol = length(wtseq))
colnames(substitution_table) <- wt_amino_acids
rownames(substitution_table) <- substituting_amino_acids
# Fill table with sd values
for(i in 1:length(SingleSubsti)) {
  for(j in 1:nchar(wt_sequence)) {
    if(substr(SingleSubsti[i],j,j) != wtseq[j]) {
      substitution_table[substr(SingleSubsti[i],j,j),j] <- d4[d4$seq==SingleSubsti[i],"sd"]
    }
  }
}
# Rearrange substitution_table for heatmap
st <- data.frame(substitution_table)
library(tidyr)
st_gathered <- gather(st)
names(st_gathered) <- c("Position", "SD") 
# Keep only the position number without the WT residue
st_gathered$Position <- substring(st_gathered$Position, 2)
# Add the substituting AA and context and rearrange the DF
st_gathered$Substitution <- substituting_amino_acids
st_gathered$Context <- context
st_gathered_T203I <- st_gathered[,c(1,3,2,4)]

#######################################################################L207V
wt_sequence <- "LPDNHYLSTQTVVSKDPN"
context <- "L207V"
# Extract sequence of the WT 
wtseq <- vector()
for(i in 1:nchar(as.character(wt_sequence))) {
  wtseq[i]<-substr(wt_sequence,i,i) 
}
# Subset dataframe to only include peptides of the same length as the wt
sameLengthAsWT <-subset(d4,nchar(as.character(d4$seq),type="chars") == nchar(as.character(wt_sequence),type="chars"))
# Create vector with single amino acid substitutions only
Boolean_string <- c()
for(m in 1:length(sameLengthAsWT$seq)) {  #for all peptides of same length as wt
  v <- unlist(strsplit(as.character(wt_sequence),character())) #wt sequence
  w <- unlist(strsplit(as.character(sameLengthAsWT$seq[m]),character()))  #sequence of the other peptides
  if(sum(v==w) == nchar(as.character(wt_sequence))-1) {  #if wt sequence and sequence of other peptide has all but one amino acids in common...
    Boolean_string[m]<-m        #returning a vector of the positions in sameLengthAsWT$Peptide meeting the criteria
  }
}
Boolean_string<-na.omit(Boolean_string)
SingleSubslv<-sameLengthAsWT$seq[Boolean_string]
# Creating a data frame with wt amino acids as column names and turning it into a vector
new_list <- list()
for(j in first_position_in_wt:(first_position_in_wt + nchar(as.character(wt_sequence)) - 1)) {
  new_list[[paste0(wtseq[j-first_position_in_wt+1],j,sep="")]] = NA
}
wt_amino_acids<-as.data.frame(new_list)
wt_amino_acids<-as.vector(names(wt_amino_acids))
# Vector with substituting amino acids
substituting_amino_acids<-unlist(strsplit(substituting_aa,character()))
# Create empty matrix with the substitution table format
no_substitutions_sd <- d4[d4$seq==wt_sequence,"sd"]
substitution_table <- matrix(data=c(no_substitutions_sd), nrow = length(substituting_amino_acids),ncol = length(wtseq))
colnames(substitution_table) <- wt_amino_acids
rownames(substitution_table) <- substituting_amino_acids
# Fill table with sd values
for(i in 1:length(SingleSubslv)) {
  for(j in 1:nchar(wt_sequence)) {
    if(substr(SingleSubslv[i],j,j) != wtseq[j]) {
      substitution_table[substr(SingleSubslv[i],j,j),j] <- d4[d4$seq==SingleSubslv[i],"sd"]
    }
  }
}
# Rearrange substitution_table for heatmap
st <- data.frame(substitution_table)
library(tidyr)
st_gathered <- gather(st)
names(st_gathered) <- c("Position", "SD") 
# Keep only the position number without the WT residue
st_gathered$Position <- substring(st_gathered$Position, 2)
# Add the substituting AA and context and rearrange the DF
st_gathered$Substitution <- substituting_amino_acids
st_gathered$Context <- context
st_gathered_L207V <- st_gathered[,c(1,3,2,4)]

#######################################################################L207R
wt_sequence <- "LPDNHYLSTQTVRSKDPN"
context <- "L207R"
# Extract sequence of the WT 
wtseq <- vector()
for(i in 1:nchar(as.character(wt_sequence))) {
  wtseq[i]<-substr(wt_sequence,i,i) 
}
# Subset dataframe to only include peptides of the same length as the wt
sameLengthAsWT <-subset(d4,nchar(as.character(d4$seq),type="chars") == nchar(as.character(wt_sequence),type="chars"))
# Create vector with single amino acid substitutions only
Boolean_string <- c()
for(m in 1:length(sameLengthAsWT$seq)) {  #for all peptides of same length as wt
  v <- unlist(strsplit(as.character(wt_sequence),character())) #wt sequence
  w <- unlist(strsplit(as.character(sameLengthAsWT$seq[m]),character()))  #sequence of the other peptides
  if(sum(v==w) == nchar(as.character(wt_sequence))-1) {  #if wt sequence and sequence of other peptide has all but one amino acids in common...
    Boolean_string[m]<-m        #returning a vector of the positions in sameLengthAsWT$Peptide meeting the criteria
  }
}
Boolean_string<-na.omit(Boolean_string)
SingleSubslr<-sameLengthAsWT$seq[Boolean_string]
# Creating a data frame with wt amino acids as column names and turning it into a vector
new_list <- list()
for(j in first_position_in_wt:(first_position_in_wt + nchar(as.character(wt_sequence)) - 1)) {
  new_list[[paste0(wtseq[j-first_position_in_wt+1],j,sep="")]] = NA
}
wt_amino_acids<-as.data.frame(new_list)
wt_amino_acids<-as.vector(names(wt_amino_acids))
# Vector with substituting amino acids
substituting_amino_acids<-unlist(strsplit(substituting_aa,character()))
# Create empty matrix with the substitution table format
no_substitutions_sd <- d4[d4$seq==wt_sequence,"sd"]
substitution_table <- matrix(data=c(no_substitutions_sd), nrow = length(substituting_amino_acids),ncol = length(wtseq))
colnames(substitution_table) <- wt_amino_acids
rownames(substitution_table) <- substituting_amino_acids
# Fill table with sd values
for(i in 1:length(SingleSubslr)) {
  for(j in 1:nchar(wt_sequence)) {
    if(substr(SingleSubslr[i],j,j) != wtseq[j]) {
      substitution_table[substr(SingleSubslr[i],j,j),j] <- d4[d4$seq==SingleSubslr[i],"sd"]
    }
  }
}
# Rearrange substitution_table for heatmap
st <- data.frame(substitution_table)
library(tidyr)
st_gathered <- gather(st)
names(st_gathered) <- c("Position", "SD") 
# Keep only the position number without the WT residue
st_gathered$Position <- substring(st_gathered$Position, 2)
# Add the substituting AA and context and rearrange the DF
st_gathered$Substitution <- substituting_amino_acids
st_gathered$Context <- context
st_gathered_L207R <- st_gathered[,c(1,3,2,4)]

#######################################################################H199Y/T203Y
wt_sequence <- "LPDNYYLSYQTVLSKDPN"
context <- "H199Y/T203Y"
# Extract sequence of the WT 
wtseq <- vector()
for(i in 1:nchar(as.character(wt_sequence))) {
  wtseq[i]<-substr(wt_sequence,i,i) 
}
# Subset dataframe to only include peptides of the same length as the wt
sameLengthAsWT <-subset(d4,nchar(as.character(d4$seq),type="chars") == nchar(as.character(wt_sequence),type="chars"))
# Create vector with single amino acid substitutions only
Boolean_string <- c()
for(m in 1:length(sameLengthAsWT$seq)) {  #for all peptides of same length as wt
  v <- unlist(strsplit(as.character(wt_sequence),character())) #wt sequence
  w <- unlist(strsplit(as.character(sameLengthAsWT$seq[m]),character()))  #sequence of the other peptides
  if(sum(v==w) == nchar(as.character(wt_sequence))-1) {  #if wt sequence and sequence of other peptide has all but one amino acids in common...
    Boolean_string[m]<-m        #returning a vector of the positions in sameLengthAsWT$Peptide meeting the criteria
  }
}
Boolean_string<-na.omit(Boolean_string)
SingleSubsht<-sameLengthAsWT$seq[Boolean_string]
# Creating a data frame with wt amino acids as column names and turning it into a vector
new_list <- list()
for(j in first_position_in_wt:(first_position_in_wt + nchar(as.character(wt_sequence)) - 1)) {
  new_list[[paste0(wtseq[j-first_position_in_wt+1],j,sep="")]] = NA
}
wt_amino_acids<-as.data.frame(new_list)
wt_amino_acids<-as.vector(names(wt_amino_acids))
# Vector with substituting amino acids
substituting_amino_acids<-unlist(strsplit(substituting_aa,character()))
# Create empty matrix with the substitution table format
no_substitutions_sd <- d4[d4$seq==wt_sequence,"sd"]
substitution_table <- matrix(data=c(no_substitutions_sd), nrow = length(substituting_amino_acids),ncol = length(wtseq))
colnames(substitution_table) <- wt_amino_acids
rownames(substitution_table) <- substituting_amino_acids
# Fill table with sd values
for(i in 1:length(SingleSubsht)) {
  for(j in 1:nchar(wt_sequence)) {
    if(substr(SingleSubsht[i],j,j) != wtseq[j]) {
      substitution_table[substr(SingleSubsht[i],j,j),j] <- d4[d4$seq==SingleSubsht[i],"sd"]
    }
  }
}
# Rearrange substitution_table for heatmap
st <- data.frame(substitution_table)
library(tidyr)
st_gathered <- gather(st)
names(st_gathered) <- c("Position", "SD") 
# Keep only the position number without the WT residue
st_gathered$Position <- substring(st_gathered$Position, 2)
# Add the substituting AA and context and rearrange the DF
st_gathered$Substitution <- substituting_amino_acids
st_gathered$Context <- context
st_gathered_DM <- st_gathered[,c(1,3,2,4)]

#Create the data frame for all contexts

st_gathered_SD <- rbind(st_gathered_WT,st_gathered_H199Y,st_gathered_T203Y,st_gathered_T203I,st_gathered_L207V,st_gathered_L207R,st_gathered_DM)

# Getting objects ready for the plot
Context <- factor(c("WT", "H199Y", "T203Y", "T203I", "L207V", "L207R","H199Y/T203Y"), levels = c("WT", "H199Y", "T203Y", "T203I", "L207V", "L207R", "H199Y/T203Y"))

Position <- factor(Position, levels = c("195", "196", "197", "198", "199", "200", "201", "202", "203", "204", "205", "206", "207", "208", "209", "210", "211", "212"))

SD <- st_gathered_SD$SD
