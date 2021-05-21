# Copyright (C) 2017-2021 Kristoffer Enoe Johansson <kristoffer.johansson@bio.ku.dk>
# This file (gmma05_analysis.r) is part of the GMMA project

options(width=300, digits=4, stringsAsFactors=F)
pcol=c(seq(6),9,10,11,12,17,19,20,21,27,28,29)

args = commandArgs(trailingOnly=TRUE)
if (interactive()) {
    file = "gmma_fit_global.rda"
    # file = "gmma_fit_global_assigned.rda"
} else if (length(args) < 1) {
    print("")
    print("usage: Rscript gmma06_analysis.r <gmma_fit_global.rda>")
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

# chi_sq_red = 1/(length(global_fit$fvec)-length(global_fit$par)) * global_fit$deviance/wt$sd^2
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

# I may change this!
B_pred_init = function(dG) { e = exp(dG/settings$RT); (cp$B_max_init + cp$B_D_init*e) /(1.0+e) }
B_pred_glob = function(dG) { e = exp(dG/settings$RT); (cp$B_max + cp$B_D*e) /(1.0+e) }

# Compare initial and global fit
quartz(width=10, height=4)
par(mfcol=c(1,2))
x=seq(-20,20,.1)
jitter = rnorm(length(mutant$signal), 0, .05)
# if (! is.null(mutant$dG_init))  {
plot(mutant$dG_init, mutant$signal+jitter, pch=".", xlim=c(-8,12),
         main=sprintf("Initial fit, dG_wt %.2f",cp$dG_wt_init), xlab="Mutant stability [kcal/mol]", ylab="Brightness")
# } else if (! is.null(mutant$ddG_init))  {
#     plot(cp$dG_wt_init+mutant$ddG_init, mutant$signal, pch=".", xlim=c(-8,12),
#          main=sprintf("Initial fit, dG_wt %.2f",cp$dG_wt_init), xlab="Mutant stability [kcal/mol]", ylab="Brightness")
# } else {
#     plot(0,0,main="Missing mutant$ddG_init or mutant$dG_init")
# }
lines(x, B_pred_init(x), col="red", lwd=2)
plot(mutant$dG_glob, mutant$signal+jitter, pch=".", xlim=c(-8,12), main=sprintf("Global fit, dG_wt %.2f",cp$dG_wt), xlab="Mutant stability [kcal/mol]", ylab="Brightness")
lines(x, B_pred_glob(x), col="red", lwd=2)
quartz.save("global_fit.png",type="png")


################################################################################
# Uncertainties from global fit
################################################################################
h = global_fit$hessian
print(sprintf("Reciprocal condition number of Hessian: %.2g",rcond(h)))
print("If this is very low, consider reparamerizing some variables to avoid covarying parameters")

# The covariance matrix may be approximated (I hope!) by the inverse of the Hessian
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
# print(rownames(subst)[si_pds])

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
    # return( sqrt(sum(res_sq[mn]) /dof /chi_sq_red) )
    # Relative mean residual: if 1 the subset of data has the same residuals as the entire model
    return( sqrt( sum(res_sq[mn])/dof / (global_fit$deviance/global_dof) ) )
}

print("Calculate data subset uncertainties")
subst$suberr = NA
subst$suberr[si_use] = sapply(si_use, calc_suberr)

# Uncertainty from number of observations. 
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
    # } else if (mutant[mi,'signal'] > cp$B_max) {
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
    
    # subst[si[sii],'stderr_meas'] = max_se
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
# print(sprintf("Total estimated std. errors %d out of %d. Total missing %d", sum(! is.na(subst$se)), nsubst, sum(is.na(subst$se))))
print(sprintf("Total estimated std. errors %d out of %d. Total missing %d", sum(! is.na(subst$stderr_meas)), nsubst, sum(is.na(subst$stderr_meas))))


################################################################################
# Substitutions with reliable ddG fit
################################################################################
# Substitutions for which we believe the fit
# se cut between 1 and 2 makes a difference
min_obs = 40
max_se = 1.5
# original selection
i_sig = which(!is.na(subst$se) & subst$se < max_se & subst$obs > min_obs)

# Brings discards from 962 to 886 with max_see=0.2
max_see = 0.2
# i_sig = which(!is.na(subst$se) & subst$se/sqrt(subst$obs) < max_see)

subst_sig = subst[i_sig,]
nsubst_sig = length(subst_sig$i)
subst_sig = subst_sig[order(subst_sig$ddG_glob),]
# subst_sig = subst_sig[order(subst_sig$ddG_glob+subst_sig$stderr_meas),]
subst_sig$rank = seq(nsubst_sig)
subst$rank = rep(NA,nsubst)
subst[rownames(subst_sig),'rank'] = subst_sig$rank

quartz(height=6, width=6)
plot(subst$se, 1/sqrt(subst$obs), col=(subst$ddG_glob<0)+1, xlim=c(0,4), ylim=c(0,.6), pch=20, main="Significant substitutions")
abline(v=max_se, h=1/sqrt(min_obs), col=8)
x=seq(0,4,.1)
lines(x,max_see/x, col=8)
points(subst_sig$se, 1/sqrt(subst_sig$obs), col=3, pch=1)
legend("bottomright", c("ddG_glob < 0","ddG_glob > 0","subst_sig","decision bounds"), col=c(2,1,3,8), pch=c(20,20,1,NA), lty=c(NA,NA,NA,1), bg="white")
quartz.save("err_obs.png",type="png")

print(sprintf("Discarding %d substitutions (out of %d) with high uncertainties",nsubst-length(subst_sig$i),nsubst))
print("Top25 substitutions")
print(head(subst_sig, n=25))


################################################################################
# Substitution effect
################################################################################
# Stability effect of each substitution: unknown, neutral, stab, destab, (very_destab, ifs)
stability_effect = rep("neu",nsubst)
stability_effect[which(is.na(subst$rank))] = "unknown"  # also covers ddG_glob==NA

# Having fixed boundaries on the neutral category is similar to defining a flat-rate uncertainty that we believe more than other uncertainties
stability_effect[which(subst$ddG_glob+subst$stderr_meas < 0 & ! is.na(subst$rank))] = "stab"
stability_effect[which(subst$ddG_glob-subst$stderr_meas > 0 & ! is.na(subst$rank))] = "destab"

# stability_effect[which(subst$ddG_glob-subst$stderr_meas > abs(cp$dG_wt) & ! is.na(subst$rank))] = "very_destab"
stability_effect[which(subst$gmma=="hanging" & subst$ddG_glob+subst$stderr_meas < 0.0)] = "stab"
stability_effect[which(subst$gmma=="hanging" & subst$ddG_glob-subst$stderr_meas > 0.0)] = "destab"
stability_effect[which(subst$ddG_glob > abs(cp$dG_wt))] = "destab"                 # Also covers subst with Inf or NA errors - consider limit 1.5 or abs(cp$dG_wt)

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
    nsubst_sig, length(ssi_stab),length(ssi_stab)*100.0/nsubst_sig, length(ssi_neut),length(ssi_neut)*100.0/nsubst_sig, length(ssi_dest),length(ssi_dest)*100.0/nsubst_sig))


################################################################################
# Plot fit
################################################################################
print("Method of ddG_glob vs. method of ddG_init")
print(xtabs(~ subst$gmma + subst$init_m))
print("Substantial fit vs. method of ddG_init")
sig = c("sig","un-sig")[is.na(subst$rank)+1]
print(xtabs(~ sig + subst$init_m))

# # Histogram of ddG_glob
# quartz(height=5, width=6)
# hist(subst$ddG_glob, breaks=seq(-5,10,.2), xlim=c(-1,5), col="gray90", main="Global ddG fit")
# h=hist(subst_sig$ddG_glob, breaks=seq(-5,10,.2), plot=F)
# plot(h, add=T, col="gray70")
# legend("topright", c("Global fit","Low uncertainty"), col=c("gray90","gray70"), pch=15)
# quartz.save("ddG_hist.png",type="png")

# How much does N average mutations destabilize?
max_nmut = max(mutant$N_sub)
x = seq(max_nmut)
ddG_nmut = sapply(x, function(n) {mean(mutant[which(mutant$N_sub==n),'dG_glob'], na.rm=T)})
ddG_nmut_sd = sapply(x, function(n) {sd(mutant[which(mutant$N_sub==n),'dG_glob'], na.rm=T)})
ddG_nmut_fit = lm(ddG_nmut ~ x)
cp$ddG_avg_init = cp$ddG_avg
cp$ddG_avg = ddG_nmut_fit$coefficients[2]
print(sprintf("ddG_avg_init: %.2f,  ddG_avg_glob: %.2f",cp$ddG_avg_init,cp$ddG_avg))

# quartz(width=6, height=6)
# plot(x, ddG_nmut, ylim=c(-5,17), xlab="Number of substitutions", ylab=expression(paste("Average ",Delta,italic(G)," [kcal/mol]")), main="Average N-mutant Stability")
# arrows(x0=x, y0=ddG_nmut-ddG_nmut_sd, y1=ddG_nmut+ddG_nmut_sd, code=3, angle=90, length=0.05)
# abline(ddG_nmut_fit)
# axis(2, at=cp$dG_wt, label="dG_wt", las=2)
# points(x, cp$ddG_avg_init*x + cp$dG_wt_init, col=2)
# legend("topleft",c(sprintf("ddG_avg fit from ddG_glob %.2f",cp$ddG_avg),sprintf("ddG_avg from dG_wt fit %.2f",cp$ddG_avg_init)), col=c(1,2), pch=1, bg="white")
# legend("bottomright",c(sprintf("Subst <ddG> = %.2f",mean(subst$ddG_glob, na.rm=T)),sprintf("Subst_sig <ddG> = %.2f",mean(subst_sig$ddG_glob))))
# quartz.save("ddG_avg.png", type="png")

# # 1-mutant plot
# require("TeachingDemos")
# quartz(height=5, width=7)
# # tiff(file="single_mut.tiff", width=90, height=60, units="mm", res=300, compression="lzw", pointsize=32)
# si = which(! is.na(subst$signal))
# plot(subst[si,'ddG_glob'], subst[si,'signal'], pch=1, xlab=expression(paste("GMMA ",Delta,Delta,italic(G)," [kcal/mol]")), ylab="Fluorescence", col=8)
# ssi = which(! is.na(subst_sig$signal))
# points(subst_sig[ssi,'ddG_glob'], subst_sig[ssi,'signal'])
# # abline(v=c(0,cp$ddG_avg), col="lightgray")
# abline(v=0, h=wt$signal, col="lightgray")
# # arrows(x0=subst[si,'ddG_glob']-subst[si,'se'], x1=subst[si,'ddG_glob']+subst[si,'se'], y0=mutant[subst[si,'signal'],'signal'],
# #        code=3, angle=90, length=0.02, col="gray")
# #subst$note = rep("", nsubst)

# plot_subst = c("S173T","R71L","P73H")
# for (mut in plot_subst) {
#     mut_split = strsplit(mut,"")
#     aa_prev = mut_split[[1]][1]
#     resi = as.numeric(paste(mut_split[[1]][2:(nchar(mut)-1)], collapse=""))
#     aa_post = mut_split[[1]][nchar(mut)]
#     mut_name = paste(aa_prev,resi+2,aa_post,sep="")
#     if (mut %in% c()) pos=1 else if (mut %in% c("R71L")) pos=2 else if (mut %in% c()) pos=4 else pos=3
#     shadowtext(subst[mut,'ddG_glob'], subst[mut,'signal'], mut_name, pos=pos, cex=.7, col=1, bg="white", font=2, r=0.1)
# }

# legend("bottomleft", c("Low uncertainty","High uncertainty","sfGFP","tsGFP","splitGFP","des11","nowGFP","PDB"), pch=c(1,1,20,20,20,20,20,20), col=c(1,8,2,6,3,4,5,8), bg="white")
# # legend("bottomleft", c("All single mutants","Uncertainty <0.2 kcal/mol","sfGFP","splitGFP","des11"), pch=c(1,1,20,20,20), col=c(8,1,2,3,4), bg="white")
# quartz.save("single_mut.png", type="png")
# # dev.off()

# quartz(height=2.5, width=4)
# par(mar=c(2,2,1,1)+.1)
# plot(subst[si,'ddG_glob'], subst[si,'signal'], pch=1, xlab="", ylab="", col=8, xlim=c(-1,.2), ylim=c(3.61,3.92))
# ssi = which(! is.na(subst_sig$signal))
# points(subst_sig[ssi,'ddG_glob'], subst_sig[ssi,'signal'])

# plot_subst = c(rownames(subst_sig[1:11,]), "N142G", "N103C", "Y37N", "Y37S", "D17E", "A108S", "L176V")
# for (mut in plot_subst) {
#     mut_split = strsplit(mut,"")
#     aa_prev = mut_split[[1]][1]
#     resi = as.numeric(paste(mut_split[[1]][2:(nchar(mut)-1)], collapse=""))
#     aa_post = mut_split[[1]][nchar(mut)]
#     mut_name = paste(aa_prev,resi+2,aa_post,sep="")
#     if (mut %in% c("S169V","V161G","I169V"))
#         pos=1
#     else if (mut %in% c("N103C","K164Q","E170A","S203T"))
#     	 pos=2
#     else if (mut %in% c("N142G","E4K","Y37S","A108S","H23Q","T60S"))
#         pos=4
#     else
#         pos=3
#     shadowtext(subst[mut,'ddG_glob'], subst[mut,'signal'], mut_name, pos=pos, cex=.7, col=1, bg="white", font=2, r=0.1)
# }
# quartz.save("single_mut_insert.png", type="png")

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

#     wt           subst na_mut nd_mut av sf do now      pdb opt opt0 ifs     fit n_max cat top
# 1    K         *EMNQRT    734    595  -  -  Q   -        A   -    -       EMQRT     9   2 144
f = file("residue.txt", "wt")
write('"# Global Multi-Mutant Analysis per position"',f)
write('"# wt: Wild-type amino acid at position"',f)
write('"# subst: Attempted substitutions at position"',f)
write('"# active: Number of active variants with this position mutated"',f)
write('"# inactive: Number of inactive variants with this position mutated"',f)
write('"# av: Amino acid in Aequorea victoria (even-wilder type)"',f)
write('"# sf,do,now: Amino acid in SuperFolderGFP, Do & Boxer split-GFP and Sarkisyan & Kondrashov NowGFP"',f)
write('"# pdb: Amino acid observed in homologs from the PDB"',f)
write('"# opt: Optimal amino acid according to GMMA"',f)
write('"# opt0: Optimal amino acid according to initial fit"',f)
write('"# ifs: Amino acids (among attempted) that inactivates irreversibly at this position "',f)
write('"# fit: Amino acids (among attempted) that resulted in an reliable GMMA fit"',f)
write('"# top: Min rank among substitutions"',f)
write('"# n_max: Max number of substitutions among active variants"',f)
write('"# cat: Mutational category of position:"',f)
write('"#     9   Hot-spot: All stabilizing, no IFS"',f)
write('"#     8   Half or more stabilizing, no IFS"',f)
write('"#     7   Less than half stabilizing and possible some IFS but most variant are active"',f)
write('"#     6   At lest one stabilizing substitution"',f)
write('"#     5   Nothing is known"',f)
write('"#     4   No fitted substitutions but some IFS"',f)
write('"#     3   Nothing is known but few active"',f)
write('"#     2   Less than half are very destabilizing and some variants are active"',f)
write('"#     1   Not all very destabilizing"',f)
write('"#     0   Conserved: All fitted substitutions are very detabilizing"',f)
residue_named = cbind(resi=seq(nres), residue)

# clean fields without data
ri = which(residue_named$wt=="x")
residue_named = residue_named[-ri,]
print(sprintf("Removed %d rows from residue data frame without substitutions: %s",length(ri),paste(ri,collapse=" ")))

write.table(residue_named,f)
close(f)
print("Dumped residue.txt")

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

print("")
ri = which(residue$cat == 0)
print("Stabil positions (resi+2)")
print(paste(residue[ri,"wt"], ri+2, sep="", collapse=", "))
print("PyMol notation (resi+2):")
print(paste(ri+2, collapse="+"))

print("")
ri = which(residue$cat == 6)
print("Tolerant positions (resi+2)")
print(paste(residue[ri,"wt"], ri+2, sep="", collapse=", "))
print("PyMol notation (resi+2):")
print(paste(ri+2, collapse="+"))

print("")
ri = which(residue$cat %in% c(7,8,9))
print("Positions with engineering potential (resi+2)")
print(paste(residue[ri,"wt"], ri+2, sep="", collapse=", "))
print("PyMol notation (resi+2):")
print(paste(ri+2, collapse="+"))

print("")
print("stabilizing substitution at unstable positions")
si = which(sapply(subst$note, function(s){length(grep(s,pattern="posU"))>.5}) & subst$eff=="stab")
si = si[order(subst[si,"ddG_glob"])]
print(paste(strsplit(wt$seq, "")[[1]][subst[si,"resi"]], subst[si,"resi"]+2, subst[si,"taa"], sep="", collapse=", "))
print("stabilizing substitution at possible unstable positions")
si = which(sapply(subst$note, function(s){length(grep(s,pattern="posPU"))>.5}) & subst$eff=="stab")
si = si[order(subst[si,"ddG_glob"])]
print(paste(strsplit(wt$seq, "")[[1]][subst[si,"resi"]], subst[si,"resi"]+2, subst[si,"taa"], sep="", collapse=", "))

# print("Effect of substitutions from Gly:")
# si = which(strsplit(wt$seq, "")[[1]][subst$resi] == "G")
# print(table(subst[si,"eff"]))

# print("Efect of substitutions from His:")
# si = which(strsplit(wt$seq, "")[[1]][subst$resi] == "H")
# print(table(subst[si,"eff"]))

print("Effect of sig. substitutions from Gly:")
si = which(strsplit(wt$seq, "")[[1]][subst_sig$resi] == "G")
print(table(subst_sig[si,"eff"]))

print("Efect of sig. substitutions from His:")
si = which(strsplit(wt$seq, "")[[1]][subst_sig$resi] == "H")
print(table(subst_sig[si,"eff"]))

################################################################################
# Dump results of post processing
################################################################################
# version=2 is necessary for R versions <3.5 to be able to read it (binf version 3.1)
save(global_fit, mut_ddG_indices, cp, settings, fit_wt, wt, mut_list,
     mutant, subst, subst_sig, residue, mut_subst_indices, subst_mut_indices, res_mut_indices,
     file="gmma_result.rda", version=2)
print("Post-processing done, dumped gmma_result.rda.")


################################################################################
# Curve plot for individual substitutions
################################################################################
plot_subst = function(si, mi_use=which(mutant$gmma=='use'), mark=c()) {
    if (is.character(si)) si = which(rownames(subst)==si)
    mi = intersect(subst_mut_indices[[si]], mi_use)
    
    s = sprintf("%d: ddG_%s = %.2f (rank=%d, init %.1f, err %.2f,%.1f) kcal/mol. Obs=%d, m=%s, gmma=%s",
      si,rownames(subst)[si],subst[si,'ddG_glob'],subst[si,'rank'],subst[si,'init_ddG'],subst[si,'stderr_meas'],subst[si,'stderr_subfit_est']*100,length(mi),subst[si,'init_m'],subst[si,'gmma'])
    plot(mutant[mi,'dG_init'], mutant[mi,'signal'], xlim=c(min(dG,-3, na.rm=T),max(dG,3, na.rm=T)), ylim=c(min(mutant[mi,'signal'],cp$B_D),max(mutant[mi,'signal'],cp$B_max)), main=s, col="red")
    points(mutant[subst_mut_indices[[si]],'dG_glob'], mutant[subst_mut_indices[[si]],'signal'], col='gray')
    points(mutant[mi,'dG_glob'], mutant[mi,'signal'])
    x = seq(-25,25,.1)
    lines(x, B_pred_glob(x), col="red")
    legend("topright",c("gmma==use","init==use","gmma!=use"), pch=1, col=c("black","red","gray"))
    for (i in seq_along(mark)) {
        sn = NA
        if (is.numeric(mark[i])) sn = rownames(subst)[mark[i]] else sn = mark[i]
	mask = sapply(mut_list[mi], function(l) {sn %in% l})
	points(dG[mask], mutant[mi,'signal'][mask], pch='x', col=i+1)
	print(sprintf("Mark %d %s has %d markings",i,sn,sum(mask)))
    }
    mii = which(mutant[mi,'residual'] > 1.0)
    if (length(mii) > 0) { text(mutant[mi[mii],'dG_glob'], mutant[mi[mii],'signal'], rownames(mutant[mi[mii],]), pos=3, cex=.8) }
    print("High residual mutants:")
    print(paste(mutant[mi,'i'][mii], collapse=","))
}

subst$sig = ! is.na(subst$rank)
table(subst[,c('gmma','eff')])
table(subst[,c('sig','eff')])
table(subst[which(subst$eff=='unknown'),'obs'])

