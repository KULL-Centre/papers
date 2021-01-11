# Copyright (C) 2017-2020 Kristoffer Enoe Johansson <kristoffer.johansson@bio.ku.dk>
# This file (gmma05_analysis.r) is part of the GMMA project

options(width=300, digits=4, stringsAsFactors=F)
pcol=c(seq(6),9,10,11,12,17,19,20,21,27,28,29)

args = commandArgs(trailingOnly=TRUE)
if (interactive()) {
    file = "gmma_fit_global.rda"
    # file = "gmma_fit_global_assigned.rda"
} else if (length(args) < 1) {
    print("")
    print("usage: Rscript gmma05_analysis.r <gmma_fit_global.rda>")
    quit(save="no")
} else {
    file = args[1]
}
print(sprintf("Read %s",file))
load(file)

# Wild-type stability
cp$dG_wt_init = cp$dG_wt
if ('dG_wt' %in% names(global_fit$par)) cp$dG_wt = global_fit$par['dG_wt']

# Maximum fluorescence base line
cp$B_max_init = cp$B_max
if ('B_max' %in% names(global_fit$par)) cp$B_max = global_fit$par['B_max']

# Minimum fluorescence base line
cp$B_D_init = cp$B_D
if ('B_D' %in% names(global_fit$par)) cp$B_D = global_fit$par['B_D']

print(sprintf("Base line initial fit: dG_wt %6.2f, B_max %6.2f, B_D %6.2f",cp$dG_wt_init,cp$B_max_init,cp$B_D_init))
print(sprintf("Base line global fit : dG_wt %6.2f, B_max %6.2f, B_D %6.2f",cp$dG_wt,cp$B_max,cp$B_D))

nmut = length(mutant[,1])
nsubst = length(subst[,1])
nres = length(residue[,1])

print("=== GMMA Settings ===")
print(settings)

################################################################################
# Global fitted ddG in subst and mutant
################################################################################
print(sprintf("Global fit result (%d): %s",global_fit$info,global_fit$message))
global_dof = length(global_fit$fvec)-length(global_fit$par)
print(sprintf("Estimated %d parameters from %d data points (%d degrees of freedom)",length(global_fit$par),length(global_fit$fvec),global_dof))
print(sprintf("Sum-of-squares %.2f in %d iterations",global_fit$deviance,global_fit$niter))

chi_sq = global_fit$deviance/wt$sd^2
print(sprintf("Chi-squared: %.3f (using WT uncertainty %.2f)",chi_sq,wt$sd))
chi_sq_red = chi_sq/global_dof
print(sprintf("Reduced chi-squared: %.3f",chi_sq_red))

subst$ddG_glob = NA
mutant$dG_glob = NA
mutant$residual = NA
si = which(subst$gmma == "use")
subst[si,'ddG_glob'] = global_fit$par[rownames(subst[si,])]
mi = which(mutant$gmma == "use")
mutant[mi,'dG_glob'] = sapply(mut_subst_indices[mi], function(l) {sum(subst[l,'ddG_glob'], na.rm=F)}) + cp$dG_wt
res = global_fit$fvec
names(res) = names(mut_ddG_indices)
mutant[mi,'residual'] = res[rownames(mutant[mi,])]
wt$residual = res["WT"]

B_pred_init = function(dG) { e = exp(dG/settings$RT); (cp$B_max_init + cp$B_D_init*e) /(1.0+e) }
B_pred_glob = function(dG) { e = exp(dG/settings$RT); (cp$B_max + cp$B_D*e) /(1.0+e) }

# Compare initial and global fit
quartz(width=10, height=4)
par(mfcol=c(1,2))
x=seq(-20,20,.1)
plot(mutant$dG_init, mutant$signal, pch=".", xlim=c(-8,12),
    main=sprintf("Initial fit, dG_wt %.2f",cp$dG_wt_init), xlab="Mutant stability [kcal/mol]", ylab="Brightness")
lines(x, B_pred_init(x), col="red", lwd=2)
plot(mutant$dG_glob, mutant$signal, pch=".", xlim=c(-8,12), main=sprintf("Global fit, dG_wt %.2f",cp$dG_wt), xlab="Mutant stability [kcal/mol]", ylab="Brightness")
lines(x, B_pred_glob(x), col="red", lwd=2)
quartz.save("global_fit.png",type="png")

get_si = function(tag, col="i", verbose=F) {
    ri = which(residue[,tag] != "-")
    if (residue[117,tag] == "L") ri = ri[ri!=117]
    if (residue[118,tag] == "V") ri = ri[ri!=118]
    sn = unname(unlist(mapply(function(wt,i,s) { paste(wt,i,unlist(strsplit(s,"")),sep="") }, residue[ri,'wt'], rownames(residue)[ri], residue[ri,tag])))
    si = which(substr(rownames(subst),1,nchar(rownames(subst))-settings$taa_letters+1) %in% sn)
    return(si)
}

# Plot ddG global vs init correlation w marked known mutations
quartz(width=6, height=6)
plot(subst$ddG_glob, subst$init_ddG, pch='.')
points(subst[get_si('pdb'),'ddG_glob'], subst[get_si('pdb'),'init_ddG'], pch='x', col=5)
points(subst[get_si('do'),'ddG_glob'], subst[get_si('do'),'init_ddG'], pch='x', col=3)
points(subst[get_si('now'),'ddG_glob'], subst[get_si('now'),'init_ddG'], pch='x', col=4)
points(subst[get_si('ts'),'ddG_glob'], subst[get_si('ts'),'init_ddG'], pch='x', col=6)
points(subst[get_si('sf'),'ddG_glob'], subst[get_si('sf'),'init_ddG'], pch='x', col=2)
legend("topleft", c("SF","TS","DoBox","Now","PDB"), pch='x', col=c(2,6,3,4,5))
quartz.save("ddG_global_vs_init.png", type="png")


################################################################################
# Uncertainties from global fit
################################################################################
h = global_fit$hessian
print(sprintf("Reciprocal condition number of Hessian: %.2g",rcond(h)))
print("If this is very low, consider reparamerizing some variables to avoid covarying parameters")

# The covariance matrix may be approximated by the inverse of the Hessian
ih = tryCatch(chol2inv(chol(h)), error = function(e) {print("Cannot calculate inverse of full hessian"); print(e); NA})

# Parameter derivatives, i.e. d-param/d-residual, in a list that match the Hessian 
pdh = sqrt(diag(ih))
names(pdh) = names(global_fit$par)
subst$se = pdh[rownames(subst)]

# Matrix inversion may result in very large numbers
subst[which(subst$se > 100),'se'] = Inf

print("Calculating standard errors of WT fit")
stderr_wt = sqrt(abs(diag(solve(fit_wt$hessian))))
n_glob_par = 0
se_glob = list()
for (tag in c('dG_wt','B_max','B_D')) {
    if (tag %in% names(pdh)) {
        se_glob[tag] = pdh[tag]
        n_glob_par = n_glob_par+1
	print(sprintf("Std error of %s estimated from global fit",tag))
    } else if (tag %in% names(stderr_wt)) {
        se_glob[tag] = stderr_wt[tag]
	print(sprintf("Std error of %s estimated from WT fit",tag))
    } else {
        print(sprintf("Cannot estimate error for %s",tag))
    }
}
print("Std. error of global parameters")
print(se_glob)

# Remove global parameters and check
if (n_glob_par > 0) pdh = pdh[-seq(length(pdh)+1-n_glob_par,length(pdh))]
i_pdh = which(sapply(names(pdh), function(n) {! n %in% rownames(subst)}))
if (length(i_pdh) > 0) { print("WARNING: Globally fitted parameter(s) not in subst nor recognised as global parameters:"); print(names(pdh)[i_pdh]) }

si_pds  = which(sapply(rownames(subst), function(n) {! n %in% names(pdh)}))
print(sprintf("Substitutions not in global fit: %d",length(si_pds)))

si_use = which(subst$gmma=="use")
mi_use = which(mutant$gmma=="use")
obs_use = sapply(si_use, function(si) { length(intersect(subst_mut_indices[[si]], mi_use)) })
obserr = 1/sqrt(obs_use)

# make sure that hessian rows match subst rows before assigning new columns in subst
stopifnot(all( rownames(subst)[si_use] == names(obserr) ))

# The uncertainty based to the variantion of the measurement of the reference
# Includes a factor 3 to get the 99.7% percentile (2 would be 95%)
# The effect categories (subst$eff) depends on this with neutral being zero +- stderrmeas
subst$stderr_meas = subst$se * wt$sd *3 

# Uncertainty of estimation
subst$stderr_meas_est = NA
subst[si_use,"stderr_meas_est"] = subst$stderr_meas[si_use] * obserr

# The uncertainty based on the fit of the model
subst$stderr_fit = subst$se * sqrt(global_fit$deviance/global_dof)
# Uncertainty of estimation
subst$stderr_fit_est = NA
subst[si_use,"stderr_fit_est"] = subst$stderr_fit[si_use] * obserr

res_sq = global_fit$fvec^2
names(res_sq) = names(mut_ddG_indices)

dpf = 1 - length(si_use)/length(mi_use)
print(sprintf("data-to-parameters factor: %.3f (should be close to one but smaller)",dpf))

# Instead of scaling the curvature with the global average residual, scale with the average residual of the points that are
#   directly coupled to the parameter: 'subset error'
calc_suberr = function(si) {
    mi = intersect(subst_mut_indices[[si]], mi_use)
    mn = rownames(mutant)[mi]
    if (length(mn) <= 2) return(NA)
    dof = length(mn)*dpf-2
    stopifnot(subst[si,"obs"] > dof)
    # Relative mean residual: if 1 the subset of data has the same residuals as the entire model
    return( sqrt( sum(res_sq[mn])/dof / (global_fit$deviance/global_dof) ) )
}

print("Calculate data subset uncertainties")
subst$suberr = NA
subst$suberr[si_use] = sapply(si_use, calc_suberr)

# Uncertainty from number of observations. This should compensate for the case where a suberr is low by chance simply because there are few points to fit.
# E.g. a suberr may be ~0.5 whereas obserr is 0.5 for n=4 and 0.25 for n=16, i.e. suberr*obserr is 0.25 for a lucky fit of 4 points and 0.25 for a normal fit of 16 points
# 4 is arbitrary because I believe that obs=16 is where this factor should be one (could also be 5...). This factor is not in the '_est' unceratinties
subst$obserr = NA
subst$obserr[si_use] = 4.0/sqrt(obs_use)

# An uncertainty to indicate if the ddG_glob value is well estimated
subst$stderr_subfit_est = subst$stderr_fit_est * subst$suberr


################################################################################
# Calculate ddG for hanging nodes 
################################################################################
# NOTE: substitutions without error estimate are discarded later!!

print(sprintf("Estimating stability and std errors for %d substitutions with gmma=='hanging'",sum(subst$gmma=='hanging')))
for (mi in which(mutant$gmma == "hanging")) {
    si = mut_subst_indices[[mi]]
    # which does not have ddG assigned
    sii = which(is.na(subst[si,'ddG_glob']))
    nest = length(sii)
    if (nest < 1) {
        # Subst in hanging variants may in rare cases be in more than one
	print(sprintf("WARNING: No substitutions in hanging variant %s needs to estimate hanging ddG",rownames(mutant)[mi]))
        next
    }
    # ddG sum of those not assigned
    ddG_other = sum(subst[si,'ddG_glob'], na.rm=T)
    # If there's a NA in the error of the ddG's the error of the calculated ddG is also NA (remove hanging subst error which is always NA)
    sum_var = sqrt(sum(subst[si[which(!is.na(subst[si,'ddG_glob']))],'stderr_meas']^2))
    # Only assign non-zero ddG if activity and the sign of dG does not match
    if (sign(0.5-mutant[mi,'active']) == sign(cp$dG_wt+ddG_other)) {
        ddG_missing = 0.0
    # Only calculate ddG in switching region
    } else if (mutant[mi,'signal'] <= cp$B_D+3*wt$sd) {
        ddG_deactivating = -2*cp$dG_wt  # Say it takes twice the stability of the protein to kill it properly 
        ddG_missing = max(ddG_deactivating - ddG_other, 0.0) # if others destabilize enough for deactivation assume no effect
    } else if (mutant[mi,'signal'] > cp$B_max-3*wt$sd) {
        ddG_missing = min(-1.0*ddG_other, 0.0)               # if others substitutions are stabilizing assume no effect
    } else {
        ddG_missing = settings$RT*log((cp$B_max-mutant[mi,'signal'])/(mutant[mi,'signal']-cp$B_D)) - cp$dG_wt - ddG_other
    }
    
    # If more substitutions are unknown, simply split missing ddG among them. In principle the two are reparametrized to one effect that cannot be split.
    subst[si[sii],'ddG_glob'] = ddG_missing/nest
    # trim to allowed ddG region
    if (ddG_missing/nest > settings$glob_ddG_max) {
        subst[si[sii],'ddG_glob'] = settings$glob_ddG_max
    } else if (ddG_missing/nest < settings$glob_ddG_min) {
        subst[si[sii],'ddG_glob'] = settings$glob_ddG_min
    }
    
    subst[si[sii],'stderr_meas'] = sum_var
    mutant[mi,"dG_glob"] = cp$dG_wt + sum(subst[si,'ddG_glob']) 
    mutant[mi,"residual"] = mutant[mi,"signal"] - B_pred_glob(mutant[mi,"dG_glob"])
    print(sprintf("Added ddG = %6.2f pm %.2f (original %.2f) to subst %4d %5s based on mutant %5d %s of brightness %.2f, dG_wt+ddG_other = %.2f, dG = %.2f and residual %.2f",
        subst[si[sii],'ddG_glob'], sum_var, ddG_missing/nest, si[sii], rownames(subst)[si[sii]], mi, paste(mut_list[[mi]], collapse=","),
	mutant[mi,'signal'], cp$dG_wt+ddG_other, mutant[mi,"dG_glob"], mutant[mi,"residual"]))
	
    # Report if more subst for this mutant was not determined in global fit (independent of how many was estimated here)
    if (sum(subst[si,'gmma']=='hanging') > 1) {
        i = which(subst[si,'gmma']=='hanging')
        print(sprintf("WARNING: Cluster of %d coupled substitutions outside global fit: %s",length(i),paste(rownames(subst[si[i],]),collapse=" ")))
    }
    # stopifnot(abs(subst[si[sii],'ddG_glob']) < 1e-9 | sign(mutant[mi,"dG_glob"]) == sign(0.5-mutant[mi,'active']) )
}
print(sprintf("Total estimated ddG's %d out of %d. Total missing %d", sum(! is.na(subst$ddG_glob)), nsubst, sum(is.na(subst$ddG_glob))))
print(sprintf("Total estimated std. errors %d out of %d. Total missing %d", sum(! is.na(subst$stderr_meas)), nsubst, sum(is.na(subst$stderr_meas))))


################################################################################
# Substitutions with reliable ddG fit
################################################################################
# Substitutions for which we believe the fit
max_sfe = 0.05
i_sig = which(!is.na(subst$stderr_meas_est) & subst$stderr_subfit_est < max_sfe)

subst_sig = subst[i_sig,]
nsubst_sig = length(subst_sig$i)
subst_sig = subst_sig[order(subst_sig$ddG_glob),]
subst_sig$rank = seq(nsubst_sig)
subst$rank = rep(NA,nsubst)
subst[rownames(subst_sig),'rank'] = subst_sig$rank


################################################################################
# Substitution effect
################################################################################
# Stability effect of each substitution: unknown, neutral, stab, destab, (very_destab, ifs)
# To be used in the position categories
stability_effect = rep("neu",nsubst)
stability_effect[which(is.na(subst$rank))] = "unknown"  # also covers ddG_glob==NA

# Having fixed boundaries on the neutral category is similar to defining a flat-rate uncertainty that we believe more than other uncertainties
stability_effect[which(subst$ddG_glob+subst$stderr_meas < 0 & ! is.na(subst$rank))] = "stab"
stability_effect[which(subst$ddG_glob-subst$stderr_meas > 0 & ! is.na(subst$rank))] = "destab"

stability_effect[which(subst$gmma=="hanging" & subst$ddG_glob+subst$stderr_meas < 0.0)] = "stab"
stability_effect[which(subst$gmma=="hanging" & subst$ddG_glob-subst$stderr_meas > 0.0)] = "destab"
stability_effect[which(subst$ddG_glob > abs(cp$dG_wt))] = "destab"  # Also covers subst with Inf or NA errors - consider limit 1.5 or abs(cp$dG_wt)

subst$eff = as.factor(stability_effect)
subst_sig$eff = subst[rownames(subst_sig),"eff"]

# Stability effect categories
ssi_stab = which(subst_sig$eff == "stab")
ssi_neut = which(subst_sig$eff == "neu")
ssi_dest = which(subst_sig$eff == "destab")
# check that all subst_sig should be in the above 3 categories
ssi = c(ssi_stab,ssi_neut,ssi_dest)
stopifnot(ssi[order(ssi)] == seq(nsubst_sig))
# stopifnot(ssi == seq(nsubst_sig))
print(sprintf("Of %d sig. subst., %d (%.1f%%)  %d (%.1f%%)  %d (%.1f%%)  are stabilizing, neutral and destabilizing respectively.",
    nsubst_sig, length(ssi_stab),length(ssi_stab)*100.0/nsubst_sig, length(ssi_neut),length(ssi_neut)*100.0/nsubst_sig, length(ssi_dest),
    length(ssi_dest)*100.0/nsubst_sig))

for (tag in c("sf","ts","do","pross","now","pdb")) {
    ssi = which(sapply(subst_sig$assign, function(s) {tag %in% strsplit(s," ")[[1]]}))
    print(sprintf("=== %s stabilizing, neutral, destabilizing and uncertain ===",tag))
    ssii = intersect(ssi,ssi_stab)
    print(paste(strsplit(wt$seq, "")[[1]][subst_sig[ssii,"resi"]], subst_sig[ssii,"resi"]+2, subst_sig[ssii,"taa"], sep="", collapse=", "))
    ssii = intersect(ssi,ssi_neut)
    print(paste(strsplit(wt$seq, "")[[1]][subst_sig[ssii,"resi"]], subst_sig[ssii,"resi"]+2, subst_sig[ssii,"taa"], sep="", collapse=", "))
    ssii = intersect(ssi,ssi_dest)
    print(paste(strsplit(wt$seq, "")[[1]][subst_sig[ssii,"resi"]], subst_sig[ssii,"resi"]+2, subst_sig[ssii,"taa"], sep="", collapse=", "))
    si = intersect(get_si(tag), which(is.na(subst$rank)))
    print(paste(strsplit(wt$seq, "")[[1]][subst[si,"resi"]], subst[si,"resi"]+2, subst[si,"taa"], sep="", collapse=", "))
}


################################################################################
# Position categories
################################################################################
residue$fit = rep("",nres)    # Which observed substitutions (residue$subst) could be fitted reliably
residue$ifs = rep("",nres)    # Which observed substitutions (residue$subst) are IFS
residue$opt = rep('-',nres)   # Optimal identity from GMMA
residue$opt0 = rep('-', nres) # Optimal identity from intial fit
residue$top = rep(NA,nres)    # Min rank among substitutions
residue$n_max = rep(0,nres)   # Max number of substitutions among active variants
residue$cat = rep(NA,nres)    # Category of position

# Assign residue cathegory
si_ifs = which(sapply(subst$note, function(s) {"ifs" %in% strsplit(s," ")[[1]]} ))
for (resi in seq(nres)) {
    # Assignments
    si = which(subst$resi==resi)
    residue[resi,'fit'] = paste(subst[intersect(si,which(!is.na(subst$rank))),'taa'], collapse="")
    residue[resi,'ifs'] = paste(subst[intersect(si,si_ifs),'taa'], collapse="")
    if (any(!is.na(subst[si,'rank']))) residue[resi,'top'] = min(subst[si,'rank'], na.rm=T)

    mi = unique(unlist(subst_mut_indices[si]))
    mii_active = which(mutant[mi,'active'] > 0.5)
    residue[resi,'n_max'] = max(c(0, mutant[mi[mii_active],'N_sub']), na.rm=T)
    si = which(subst$resi==resi)
    if (length(si) > 0) {
        if (min(subst[si,'init_ddG'], na.rm=T) < 0.0) residue[resi,'opt0'] = subst[ si[which.min(subst[si,'init_ddG'])], 'taa' ]
        if (any(subst[si,'eff'] == "stab")) {
	    si_stab = si[which(subst[si,'eff']=="stab")]
	    residue[resi,'opt'] = subst[ si_stab[which.min(subst[si_stab,'ddG_glob'])], 'taa' ]
	}
    }

    t = table(subst[si,"eff"])
    t_known_sum = sum(t[c("stab","destab","neu")])
    if      (t["stab"]   > 1)                              residue[resi,'cat'] = 9  # more than one stabilizing
    else if (t["stab"]   > 0)                              residue[resi,'cat'] = 8  # any stabilizing
    else if (t["neu"]    > 0 & t["neu"]    >  t["destab"]) residue[resi,'cat'] = 6  # tolerant position with most neutral
    else if (t["destab"] > 0 & t["destab"] >= t["neu"])    residue[resi,'cat'] = 0  # stable position with most destabilizing
    else if (t["unknown"] == sum(t))                       residue[resi,'cat'] = 5  # this covers sum(t)==0
    else                                                   residue[resi,'cat'] = 4  # what didn't I think about?
    
}

f = file("residue.csv", "wt")
write('"# Global Multi-Mutant Analysis per position"',f)
write('"# wt: Wild-type (reference) amino acid at position"',f)
write('"# subst: Attempted substitutions at position"',f)
write('"# N_mut: Sum of observations for each substitution at position"',f)
write('"# HSE and asa: Solvent accesibility given as half-sphere exposure (HSE) and solvent accessible surface area (ASA)"',f)
write('"# burial2 and burial3: Solvent accesibility category with 2 or 3 categories based on ASA percentage (DSSP)"',f)
write('"# av: Amino acid in Aequorea victoria (even-wilder type)"',f)
write('"# av,sf,ts,do,now: Amino acid in Aequorea victoria, SuperFolderGFP, T-Sapphire, Do & Boxer split-GFP and Sarkisyan & Kondrashov NowGFP"',f)
write('"# pdb: Amino acid observed in homologs from the PDB"',f)
write('"# active and inactive: Number of active and inactive variants with this position mutated"',f)
write('"# ifs: Amino acids (among attempted) that inactivates irreversibly at this position "',f)
write('"# fit: Amino acids (among attempted) that resulted in an reliable GMMA fit"',f)
write('"# opt: Optimal amino acid according to GMMA"',f)
write('"# opt0: Optimal amino acid according to initial fit"',f)
write('"# top: Min rank among substitutions"',f)
write('"# n_max: Max number of substitutions among active variants"',f)
write('"# cat: Mutational category of position. Zero is only destabilizing substitutions and 9 is only stabilizing substitutions"',f)
residue_named = cbind(resi=seq(nres), residue)

write.table(residue_named, sep=";", row.names=F, file=f)
close(f)
print("Dumped residue.csv")

for (ri in seq(nres)) {
    si = which(subst$resi == ri)
    ssi = which(subst_sig$resi == ri)
    if (residue[ri,"cat"] == 0) {
        subst[si,"note"] = paste(subst[si,"note"], "posS", sep=" ")
        subst_sig[ssi,"note"] = paste(subst_sig[ssi,"note"], "posS", sep=" ")
    } else if (residue[ri,"cat"] %in% c(7,8,9)) {
        subst[si,"note"] = paste(subst[si,"note"], "posU", sep=" ")
        subst_sig[ssi,"note"] = paste(subst_sig[ssi,"note"], "posU", sep=" ")
    } else if (residue[ri,"cat"] %in% c(6)) {
        subst[si,"note"] = paste(subst[si,"note"], "posT", sep=" ")
        subst_sig[ssi,"note"] = paste(subst_sig[ssi,"note"], "posT", sep=" ")
    }
}

# Dump subst and mut_list in excel raedable format
subst_named = cbind(s=rownames(subst), subst)
write.table(subst_named, "subst.csv", sep=";", row.names=F)
write.table(mutant, "mutant.csv", sep=";", row.names=F)
print("Dumped subst.csv and mutant.csv")


################################################################################
# Dump results of post processing
################################################################################
# version=2 is necessary for R versions <3.5 to be able to read it (binf version 3.1)
save(global_fit, mut_ddG_indices, cp, settings, fit_wt, wt, mut_list,
     mutant, subst, residue, mut_subst_indices, subst_mut_indices, res_mut_indices,
     file="gmma_result.rda", version=2)
print("Post-processing done, dumped gmma_result.rda.")
