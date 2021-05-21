# Copyright (C) 2017-2021 Kristoffer Enoe Johansson <kristoffer.johansson@bio.ku.dk>
# This file (gmma04_global_estimation.r) is part of the GMMA project

options(width=200, digits=4, stringsAsFactors=F)

file="gmma_graph.rda"
print(sprintf("Read %s",file))
load(file)
subst_graph = subst
mutant_graph = mutant

file="gmma_fit_init.rda"
print(sprintf("Read %s",file))
load(file)

# All we need from the graph
subst$gmma = subst_graph$gmma
subst$gmma_reparam = subst_graph$gmma_reparam
mutant$gmma = mutant_graph$gmma
# Clean up
rm(subst_graph,mutant_graph)

# Load the more robust Levenberg-Marquardt non-linear least squares algorithm, nls.lm
require("minpack.lm")

nmut = length(mutant[,1])
nsubst = length(subst[,1])
nres = length(residue[,1])

print("Force cp$B_max = 1.0 and cp$B_D = 0.0")
cp$B_max = 1.0
cp$B_D = 0.0

settings$glob_ddG_max =  10.0
settings$glob_ddG_min = -5.0


################################################################################
# Data structures for global fit
################################################################################
print("Re-structuring cleaned data for global fit")

si = which(subst$gmma=="use" & subst$init_m=="ifs" & is.na(subst$init_ddG))
if (length(si) > 0) { print(sprintf("Setting %d IFS ddG's to %.2f (settings$glob_ddG_max)",length(si),settings$ddG_ifs)); subst[si,"init_ddG"] = settings$glob_ddG_max }
# the rest is set to ddG_avg
si = which(subst$gmma=="use" & is.na(subst$init_ddG))
if (length(si) > 0) { print(sprintf("Setting %d un-initialized substitution ddG's to %.2f (ddG_avg)",length(si),cp$ddG_avg)); subst[si,"init_ddG"] = cp$ddG_avg }

# initial values of fitting parameters
si = which(subst$gmma == "use")
param_init = subst[si,'init_ddG']
names(param_init) = rownames(subst[si,])
nparam = length(si)

# data to fit
mutant["mut00001","gmma"] = "use" # WT will always end up as disconnected in the gaph analysis
mi = which(mutant$gmma == "use")
# B = c(mutant[mi,'signal'],wt$signal)
B = mutant[mi,'signal']
nB = length(B)

# reverse lookup list for subst indices
si_rev = rep(NA,nsubst)
si_rev[si] = seq(length(si))

# List that gives parameter(subst ddG) indices to add per data(variant)
mut_ddG_indices = list()
for (i in seq(length(mi))) { mut_ddG_indices[[i]] = si_rev[mut_subst_indices[[mi[i]]]] }
# # no parameters for the WT
# mut_ddG_indices[[length(mut_ddG_indices)+1]] = numeric(0)
# names(mut_ddG_indices) = c(rownames(mutant[mi,]),"WT")
names(mut_ddG_indices) = rownames(mutant[mi,])


################################################################################
# Re-parameterize
################################################################################
nparam_reparam = 0
groups = na.omit(unique(subst$gmma_reparam))
for (gi in groups) {
    si = which(subst$gmma_reparam==gi)
    if (all(subst[si,'gmma'] != "use")) next
    # All parameters to be reparameterized should belong to the same global fit
    stopifnot(all(subst[si,'gmma'] == subst[si[1],'gmma']))

    name = paste(rownames(subst)[si], collapse="_")
    init_ddG = sum(subst[si,'init_ddG'])
    param_init[name] = init_ddG
    pi = sapply(si, function(i) {which(names(param_init)==rownames(subst)[i])})
    param_init = param_init[-pi]
    
    # Adapt mut_ddG_indices!
    
    nparam = nparam-2
    nparam_reparam = nparam_reparam+1
}

################################################################################
# Add epistasis parameters
################################################################################
print("Adding custom parameters")
nparam_add = 0

# # # Random example
# # > which(sapply(mut_ddG_indices, function(l) {690 %in% l & 884 %in% l}))
# # [1]  3643 41037
# # > names(param_init)[690]
# # [1] "D131G"
# # > names(param_init)[884]
# # [1] "G158V"
# # > names(param_init)[mut_ddG_indices[[3643]]]
# # [1] "D74G"  "D131G" "G158V"
# # > names(param_init)[mut_ddG_indices[[41037]]]
# # [1] "V10E"  "V27A"  "D131G" "G158V" "V191A" "L193R"
# param_init['D131GxG158V'] = 0.0
# i = which(names(param_init)=='D131GxG158V')
# mut_ddG_indices[[3643]] = append(mut_ddG_indices[[3643]], i)
# mut_ddG_indices[[41037]] = append(mut_ddG_indices[[41037]], i)
# nparam_add = nparam_add+1

print(sprintf("Added %d parameter(s)",nparam_add))


################################################################################
# Explicit Jacobian
################################################################################
# The Jacobian can be heavily simplified by the following observations
# The Jacobian d_mut/d_subst is sparse i.e. the predicted fluorescence of a mutant only changes if the involved substitutions ddG's changes. I call this delta_subst
# The Jacobian elements of a given mutant are either zero or one number because all substitutions changes dG the same amount d_mut/d_dG. Global baseline parameters are exeptions.
# Vector containing the dG_wt gradients which are identical to substitution gradients

B_expr = expression((B_max+B_D*exp(dG/RT_kcal))/(1+exp(dG/RT_kcal)))
B_DB_max = D(B_expr, "B_max")
B_DB_D = D(B_expr, "B_D")
# B_DdG = D(B_expr, "dG")

B_max_grad = function(dG) {
    local_cp = cp
    local_cp$dG = dG
    eval(B_DB_max, envir=local_cp)
}
B_D_grad = function(dG) {
    local_cp = cp
    local_cp$dG = dG
    eval(B_DB_D, envir=local_cp)
}
ddG_grad = function(param, dG) {
    e = exp(dG/settings$RT)
    ( cp$B_D*e/(1+e) - (cp$B_max*e+cp$B_D*e*e)/((1+e)*(1+e)) )/settings$RT
    # ( param['B_D']*e/(1+e) - (param['B_max']*e+param['B_D']*e*e)/((1+e)*(1+e))  )/settings$RT
}

# This is fast!
# mi/mj is row/col matrix index (par/data) and Ji is the vector index
ndata = length(B)
npar = length(param_init)
i = mapply(function(mi,mj) {mj+(mi-1)*ndata}, mut_ddG_indices, seq(ndata))
Ji = unlist(i)
gradi = unlist(mapply(function(Ji,mj) {rep(mj,length(Ji))} ,i,seq(ndata)))

jacobian = function(param, observed) {
    if ('dG_wt' %in% names(param)) {
        dG = sapply(mut_ddG_indices, function(l) {param['dG_wt']+sum(param[l])})
    } else {
        dG = sapply(mut_ddG_indices, function(l) {cp$dG_wt+sum(param[l])})
    }
    grad = -ddG_grad(param, dG)
    ndata = length(observed)
    npar = length(param)
    # only grad need s to be updated!
    J_vec = rep(0.0, npar*ndata) # A sparse vector would be good here but that's not supported by nls.lm
    J_vec[Ji] = grad[gradi]
    if ('dG_wt' %in% names(param)) { Ji = (which(names(param)=='dG_wt')-1)*ndata; J_vec[(Ji+1):(Ji+ndata)] = grad}    
    # if ('B_max' %in% names(param)) { Ji = (which(names(param)=='B_max')-1)*ndata; J_vec[(Ji+1):(Ji+ndata)] = -B_max_grad(dG)}
    # if ('B_D' %in% names(param))   { Ji = (which(names(param)=='B_D')-1)*ndata;   J_vec[(Ji+1):(Ji+ndata)] = -B_D_grad(dG)}
    return(J_vec)
}


################################################################################
# Global fit
################################################################################
# Get initial parameter values within upper and lower bounds of fit
print(sprintf("Shifted %d initial ddG's to below %.2f with a max shift of %.2f",sum(subst$init_ddG > settings$glob_ddG_max, na.rm=T), settings$glob_ddG_max, max(subst$init_ddG, na.rm=T)-settings$glob_ddG_max))
print(sprintf("Shifted %d initial ddG's to above %.2f with a max shift of %.2f",sum(subst$init_ddG < settings$glob_ddG_min, na.rm=T), settings$glob_ddG_min, -min(subst$init_ddG, na.rm=T)+settings$glob_ddG_min))
param_init = sapply(param_init, max, settings$glob_ddG_min+0.001)
param_init = sapply(param_init, min, settings$glob_ddG_max-0.001)

# Add the global parameters
# param_init['dG_wt'] = cp$dG_wt
# param_init['B_max'] = cp$B_max
# param_init['B_slope'] = 0.0
# param_init['B_D'] = cp$B_D

# Bounds for the fit
# bound_up  = c(rep(settings$glob_ddG_max, nparam), rep(settings$glob_ddG_max*2, nparam_reparam), rep(+5,nparam_add), +10.0, 100.0, +1.0, 2.5)
# bound_low = c(rep(settings$glob_ddG_min, nparam), rep(settings$glob_ddG_min*2, nparam_reparam), rep(-5,nparam_add), -10.0,   2.5, -1.0, 1.0)
# bound_up  = c(rep(settings$glob_ddG_max, nparam), rep(settings$glob_ddG_max*2, nparam_reparam), rep(+5,nparam_add), +10.0)
# bound_low = c(rep(settings$glob_ddG_min, nparam), rep(settings$glob_ddG_min*2, nparam_reparam), rep(-5,nparam_add), -10.0)
bound_up  = c(rep(settings$glob_ddG_max, nparam), rep(settings$glob_ddG_max*2, nparam_reparam), rep(+5,nparam_add))
bound_low = c(rep(settings$glob_ddG_min, nparam), rep(settings$glob_ddG_min*2, nparam_reparam), rep(-5,nparam_add))

# cp$slope = -0.02

# Only optimize dG_wt if this was optimized during param. init.
if ('dG_wt' %in% names(fit_wt$par)) {
    param_init['dG_wt'] = cp$dG_wt
    bound_up  = append(bound_up,   10.0)
    bound_low = append(bound_low, -10.0)
}

# Loss function
B_pred = function(param) {
    if ('dG_wt' %in% names(param)) {
        dG = sapply(mut_ddG_indices, function(l) {param['dG_wt']+sum(param[l])})
    } else {
        dG = sapply(mut_ddG_indices, function(l) {cp$dG_wt+sum(param[l])})
    }
    e = exp(dG/settings$RT)
    (cp$B_max + cp$B_D*e) /(1.0+e)
    # (cp$B_max+cp$slope*dG + (cp$B_D+cp$slope*dG)*e) /(1.0+e)
    # (param['B_max'] + param['B_slope']*dG + param['B_D']*e) /(1.0+e)
}

B_residual = function(param, observed) {observed - B_pred(param)}
ctrl = nls.lm.control(nprint=5, maxiter=1000)

res = B_residual(param_init, B)
res_cut = (cp$B_max-cp$B_D)/2
print(sprintf("Observations with high (>%.1f) initial residual: %d",res_cut,sum(res > res_cut)))

# do the global fit
print("Global fit")
#global_fit = nls.lm(par=param_init, fn=B_residual,               observed=B, lower=bound_low, upper=bound_up, control = ctrl)
#global_fit = nls.lm(par=param_init, fn=B_residual, jac=jacobian, observed=B,                                  control = ctrl)
global_fit = nls.lm(par=param_init, fn=B_residual, jac=jacobian, observed=B, lower=bound_low, upper=bound_up, control = ctrl)

# save
# version=2 is necessary for R versions <3.5 to be able to read it (binf version 3.1)
save(global_fit, mut_ddG_indices, cp, settings, fit_wt, wt, mut_list, mutant, subst, residue, mut_subst_indices, subst_mut_indices, res_mut_indices, file="gmma_fit_global.rda", version=2)
print(sprintf("Global fit done. info=%d",global_fit$info))


