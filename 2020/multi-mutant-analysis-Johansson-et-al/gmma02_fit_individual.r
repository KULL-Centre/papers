# Copyright (C) 2017-2020 Kristoffer Enoe Johansson <kristoffer.johansson@bio.ku.dk>
# This file (gmma02_fit_individual.r) is part of the GMMA project

options(width=200, digits=4, stringsAsFactors=F)

args = commandArgs(trailingOnly=TRUE)
# Allow fixed dG_wt to enable gmma along a line of dG_wt
if (interactive()) {
    dG_wt = NA
} else if (length(args) > 0) {
    dG_wt = as.numeric(args[1])
    print(sprintf("Using fixed dG_wt = %.2f",dG_wt))
} else {
    dG_wt = NA
    print("Will estimate dG_wt")
}

file="gmma_structured_assigned.rda"
print(sprintf("Read %s",file))
load(file)

nmut = length(mutant[,1])
nsubst = length(subst[,1])
nres = length(residue[,1])

################################################################################
# Settings
################################################################################
settings$RT = 1.985*300/1000.0   
settings$RT_unit = "kcal/mol (at 300 K)"
# settings$RT = 8.314*300/1000.0
# settings$RT_unit = "kJ/mol (at 300 K)"

settings$dG_wt_input = dG_wt

# WT Curve parameters in cp
cp = list(dG_wt=settings$dG_wt_input)

################################################################################
# Assign variants as functional or non-functional
################################################################################
mutant$active = NA    # is variant active
min_sig = min(mutant$signal)
max_sig = max(mutant$signal)
cp$B_mid = min_sig + (max_sig-min_sig)/2.0
print(sprintf("Signal below B_mid=%.2f considered inactive, the rest are considered active",cp$B_mid))

mutant$active = 1.0
mutant[which(mutant$signal < cp$B_mid),"active"] = 0.0

# For subst$active and subst$inactive counts, interpret mutant$active as active if >= active_cut and inactive if < inactive_cut
settings$active_cut = 0.5
settings$inactive_cut = settings$active_cut

subst$active = NA     # number of mutants measured active
subst$inactive = NA   # number of mutants measured inactive
for (si in seq(nsubst)) {
    mi = subst_mut_indices[[si]]
    subst[si,'active'] = sum(mutant[mi,'active'] >= settings$active_cut)
    subst[si,'inactive'] = sum(mutant[mi,'active'] < settings$inactive_cut)
}

residue$active = NA   # Number of active mutants involving this position
residue$inactive = NA # Number of active mutants involving this position
for (i in seq(nres)) {
    residue[i,'active'] = sum(subst[subst$resi==i, 'active'])
    residue[i,'inactive'] = sum(subst[subst$resi==i, 'inactive'])
}

################################################################################
# Clean data for initial fit
################################################################################
mutant$init = factor("use")
subst$init_m = factor("none")
subst$init_ddG = NA

# Function to set mutant$init_use, subst$init_m and possible subst$init_ddG
exclude_data = function(si, tag, overwrite=FALSE) {
    if (length(si)==0) return(0) # do not add level if not used
    # Make local copies of objects for use in this function
    subst_loc = subst
    mutant_loc = mutant
    # Update factor levels if necessary
    if (! tag %in% levels(subst_loc$init_m)) {
        levels(subst_loc$init_m) = c(levels(subst_loc$init_m), tag)
        stopifnot(! tag %in% levels(mutant_loc$init))
        levels(mutant_loc$init) = c(levels(mutant_loc$init), tag)	
    }
    mutants_removed = 0
    if (overwrite) {
        subst_loc[si,"init_m"] = tag
        mi = unique(unlist(subst_mut_indices[si]))
	mutant_loc[mi,'init'] = tag
	mutants_removed = length(mi)
    } else {
        sii = which(subst_loc[si,"init_m"] == "none")
        subst_loc[si[sii],"init_m"] = tag
        mi = unique(unlist(subst_mut_indices[si[sii]]))
	mii = which(mutant_loc[mi,"init"] == "use")
	mutant_loc[mi[mii],"init"] = tag
	mutants_removed = length(mi)
    }
    # Append the tag to the note-string of the substitutions, init_m may be overwritten later
    subst_loc[si,'note'] = paste(subst_loc[si,'note'], tag)
    # Update global objectes
    assign("subst", subst_loc, envir=.GlobalEnv)
    assign("mutant", mutant_loc, envir=.GlobalEnv)
    return(mutants_removed)
}

si_nonsense = which(substr(subst$taa,1,1)=="*")
mutants_removed = exclude_data(si_nonsense, "nonsense")
print(sprintf("Excluded %d nonsense substitutions and %d dependent variants", length(si_nonsense), mutants_removed))

settings$init_min_obs = 0
si_lowobs = which(subst$obs < settings$init_min_obs)
mutants_removed = exclude_data(si_lowobs, "low_obs")
print(sprintf("Excluded %d substitutions with fewer than %d observations and %d dependent variants",length(si_lowobs), settings$init_min_obs, mutants_removed))

settings$init_min_active = 0
si_lowact = which(subst$active < settings$init_min_active)
mutants_removed = exclude_data(si_lowact, "low_act")
print(sprintf("Excluded %d substitutions with fewer than %d active observations and %d dependent variants",length(si_lowact), settings$init_min_active, mutants_removed))

################################################################################
# Probability of substitution being irreversible fatal
################################################################################
subst$p_mn = NA      # multinomial probability to observe set of active and inactive N-mutants
subst$p_renom = NA   # p_mn times the number of states of p_mn

# Background probability in data of N-mutants to be inactive
max_n_mut = max(mutant$N_sub) # There are 2 15-mutants, 6 14-mutants, and 25 13-mutants all dead. 2 12-mutants are active.
p_inact = rep(0,2*max_n_mut)
for (i in seq(max_n_mut)) {
    mi = which(mutant$N_sub == i)
    f_active = sum(mutant[mi,'active'] >= settings$active_cut)
    f_inact = sum(mutant[mi,'active'] < settings$inactive_cut)
    p_inact[2*i-1] = f_active/nmut
    p_inact[2*i] = f_inact/nmut
}

for (si in seq(nsubst)) {
    mi = subst_mut_indices[[si]]

    # Propability to observe set of mutants
    x = c()
    for (n in seq(max_n_mut)) {
        i = which(mutant[mi,'N_sub']==n)
        x = append(x,c( sum(mutant[mi[i],'active'] >= settings$active_cut), sum(mutant[mi[i],'active'] < settings$inactive_cut)))
        # x_old = append(x,c( sum(mutant[mi[i],'signal'] > cp$B_mid), sum(mutant[mi[i],'signal'] < cp$B_mid)))
	# stopifnot(all(x==x_old))
    }
    subst[si,'p_mn'] = dmultinom(x,sum(x),p_inact)

    # Re-normalization
    # The source of the problem is that the number of states i very different for different n and the multinomial is normalized so the
    # sum of density over all states =1. So an alternative normalization is sum_states mn = number of states. p(n) = 1/n_states(n)
    # The number of states for a multinomial is the number of combinations for k indistinguishable balls (observed mutants) into n distinguishable
    # boxes (mutant chategory: active 1-mut, inactive 1-mut, active 2-mut, ...) without exclusion (no restrictions on the categories):
    k = sum(x) # subst[si,'obs']   # Careful, this was previously called 'n'
    n = length(p_inact)
    n_states = choose(n+k-1,k)
    subst[si,'p_renom'] = subst[si,'p_mn']*n_states
}

# Mark data for irreversible fatal IFS
settings$ifs_min_active = 1
settings$ifs_min_prob = 1e-3

si_ifs = which(subst$active < settings$ifs_min_active & subst$p_renom < settings$ifs_min_prob)
mutants_removed = exclude_data(si_ifs, "ifs")
print(sprintf("Marked %d IFS with less than %d active variants observed and renormalized multinomial probability less than %.1e and %d dependent variants",
    length(si_ifs), settings$ifs_min_active, settings$ifs_min_prob, mutants_removed))

# Dump to stdout
print("subst$ini_m:")
print(table(subst$init_m))
print("mutant$init:")
print(table(mutant$init))

# How is p_renom distributed on substitution observations and active observations
# Try detecting IFS with a minimum _fraction_ of active, e.g. W55R and L218P seems IFS but are not with 589 and 345 obs and 3 and 2 active.
#   The poor ddG of these means that the active variants of these have high residuals, i.e. misclassified as inactive.
quartz(height=6, width=7)
plot(subst$obs, subst$p_renom, log="y", pch='.', main="Likelihood of observations per substitution")#, xlim=c(0,50), ylim=c(1e-5,1e+5))
si = which(subst$init_ddG > 9.9)
points(subst[si,'obs'], subst[si,'p_renom'], col=2, pch='.', cex=2)
low = 3
si = which(subst$active < low)
n_low = length(si)
points(subst[si,'obs'], subst[si,'p_renom'], col=2)
n = c()
i = seq(low,10)
for (ii in i) {
    si = which(subst$active == ii)
    n = append(n,length(si))
    points(subst[si,'obs'], subst[si,'p_renom'], col=ii-low+3)
}
abline(h=settings$ifs_min_prob, col="lightgray")
legend("bottomright", c("All","High bound",paste("active < ",low," (",n_low,")", sep=""),paste(i," active (",n,")", sep="")), 
       pch=c(".",".",rep('o',1+length(i))), col=c(1,2,2,i-low+3), pt.cex=c(2,3,rep(1,1+length(i))), bg="white")


################################################################################
# Methods and re-parameterization for the individual estimations
################################################################################
# No more data exclusion from this point so check if any substitutions have all of their data excluded
levels(subst$init_m) = c(levels(subst$init_m), "data_excluded")
n_use = sapply(subst_mut_indices, function(l) {sum(mutant[l,"init"] == "use")})
subst[which(n_use==0 & subst$init_m=="none"),"init_m"] = "data_excluded"

# methods and error signals for the ddG estimations
init_m_excluded_levels = levels(subst$init_m)[-1]                                               # All levels until now is excluding data except the first "none" level
init_m_fit_levels = c("1_curve", "2_curve_ddG1", "3_curve_ddG2",                                # Levels for fitted ddG's
                      "4_all_flat", "5_most_flat", "6_avg_individal")
init_m_err_levels = c("ce_few_Bavg", "ce_no_inact", "ce_no_act", "ce_fit_info", "ce_high_dev",   # Fitting error signals for development: ce_ is "curve error" and ee_ is "estimation error"
                     "ee_poor_fit", "ee_bad_data", "ee_few_data")                               
levels(subst$init_m) = c("none", init_m_excluded_levels, init_m_fit_levels, init_m_err_levels)   # Careful with the order here since some subst are already assigned

subst$init_reparam = NA
print("Looking for substitutions that always co-occur")
reparam_group = 0
for (si in seq(nsubst)) {
    # reparam may have been set in previous iterations
    if (! is.na(subst[si,"init_reparam"])) { next }
    # if single mutant is observed, a ddG can always be estimated from dG_wt (but may be isolated)
    if (! is.na(subst[si,"signal"])) { next }

    mi = subst_mut_indices[[si]]
    mi = mi[which(mutant[mi,"init"]=="use")]
    si_neighbor_list = setdiff(unlist(mut_subst_indices[mi]), si)
    for (si_neighbor in si_neighbor_list) {
        si_neighbor_mi = subst_mut_indices[[si_neighbor]]
        si_neighbor_mi = si_neighbor_mi[which(mutant[si_neighbor_mi,"init"]=="use")]
	if (length(mi) != length(si_neighbor_mi)) { next }
	if (all(mi == si_neighbor_mi)) {
	    if (is.na(subst[si,"init_reparam"])) {
	        reparam_group = reparam_group+1
		subst[si,"init_reparam"] = reparam_group
		print(sprintf("New reparam group %d with %s",reparam_group,rownames(subst[si,])))
	    }
	    subst[si_neighbor,"init_reparam"] = reparam_group
	    print(sprintf("Add %s to reparam group %d",rownames(subst[si_neighbor,]),reparam_group))
	}
    }
}
print(sprintf("Found %d groups of substitutions to set for reparameterization",reparam_group))


################################################################################
# First fit wild-type deactivation free energy, dG_wt, and the average single mutation ddG, ddG_avg
################################################################################
print("")
print("Fit WT stability to all data assuming an average destabilizing effect of substitutions")

# Load the more robust Levenberg-Marquardt non-linear least squares algorithm, nls.lm
require("minpack.lm")

if (is.na(settings$dG_wt_input)) {
    stopifnot(is.na(cp$dG_wt))
    init_wt = list(B_D=1.5, dG_wt=0, ddG_avg=1.7)
    stability_wt = function(param, n_mut) {
        e=exp((param$dG_wt + param$ddG_avg*n_mut)/settings$RT)
        B_max = wt$signal + (wt$signal-param$B_D)*exp(param$dG_wt/settings$RT)
	(B_max+param$B_D*e)/(1+e)
    }
} else {
    init_wt = list(B_D=1.5, ddG_avg=1.7)
    stability_wt = function(param, n_mut) {
        e=exp((cp$dG_wt + param$ddG_avg*n_mut)/settings$RT)
        B_max = wt$signal + (wt$signal-param$B_D)*exp(cp$dG_wt/settings$RT)
	(B_max+param$B_D*e)/(1+e)
    }
}

# Multi-mutant data structure
use_mask = mutant$init == "use"
mi_use = which(use_mask)
max_n_mut = max(mutant[mi_use,"N_sub"])
d = data.frame(table(mutant[mi_use,"N_sub"]), freq_act=0)
colnames(d) = c("n_mut","obs","freq_act")
rownames(d) = d$n_mut
d$n_mut = as.numeric(levels(d$n_mut)[d$n_mut])
t = table(mutant[which(mutant$active > .9),'N_sub']) #assume mutant$active is either 1 or 0
d[names(t),'freq_act'] = t
d$dens_lib = d$obs/sum(d$obs)
d$frac_act = d$freq_act/d$obs
d$signal_avg = sapply(d$n_mut, function(n) { mean(mutant[which(mutant$N_sub == n & use_mask),"signal"]) })
# The median is rather sensitive for a bimodal distribution. On the other hand, it forces the curve onto an
#     observed value which is what the global analysis sees
# Why average here at all? If I fit to individual data, 1000 4-mutants will weight 1000 times more than
#     one 3-mutant which is currently not the case?! Also in individual curve fits.
d$signal_med = sapply(d$n_mut, function(n) { median(mutant[which(mutant$N_sub == n & use_mask),"signal"]) })
stopifnot(all((d$obs==0) == (is.na(d$signal_avg))))
if (any(is.na(d$signal_avg))) d = d[-which(is.na(d$signal_avg)),]

# Plot library composition versus fraction of active mutants
dp = d[as.character(seq(max_n_mut)),] # skip n_mut zero, i.e. wt
quartz(height=4, width=6)
plot(dp$n_mut, dp$dens_lib, ylim=c(0,1), xlab="N subst", ylab="", type="o", col=2)
points(dp$n_mut, dp$frac_act, pch=2, type="o", col=4)
legend("topright",c("Library composition","Fraction active"),pch=c(1,2),col=c(2,4))
quartz.save("library.png", type="png")

# Plot signal distributions
quartz(height=5, width=9)
breaks = seq(min(mutant[mi_use,"signal"])-.1, max(mutant[mi_use,"signal"])+.1, .02)
plot(0,0, xlim=c(1,4.2), ylim=c(0,4), xlab="signal", ylab="density")
d$signal_mode = NA
for (nmut in d$n_mut) {
    mi = which(mutant$N_sub==nmut & use_mask)
    h = hist(mutant[mi,'signal'], breaks=breaks, plot=F)
    d[as.character(nmut),"signal_mode"] = mean(h$mids[which(h$counts==max(h$counts))])
    if (nmut %in% seq(8)) {
    lines(h$mids, h$density, col=nmut, lwd=1)
    points(d[as.character(nmut),"signal_avg"], 1, cex=1.5, pch=1, col=nmut)
    points(d[as.character(nmut),"signal_med"], .7, cex=1.5, pch=2, col=nmut)
    points(d[as.character(nmut),"signal_mode"], .4, cex=1.5, pch=3, col=nmut)
    }
}
# legend("top",c(paste(seq(8),"mut"),"Avg"), col=c(seq(8),1), pch=c(rep(NA,8),"|"), pt.cex=1.5, lty=c(rep(1,8),NA), lwd=2, ncol=3)
legend("top",paste(seq(8),"mut"), col=seq(8), lty=1, lwd=2, ncol=3)
legend("topright", c("Average","Median","Mode"), pch=seq(3))

# Fit B_max, B_D, dG_wt and <ddG>to all data
residual_stability_wt = function(param, observed, x) {observed - stability_wt(param, x)}
fit_wt = nls.lm(par=init_wt, fn=residual_stability_wt, observed=d$signal_avg, x=d$n_mut)
print(summary(fit_wt))

# The fitted values to use in the following
if (is.na(settings$dG_wt_input)) cp$dG_wt = fit_wt$par$dG_wt
cp$ddG_avg = fit_wt$par$ddG_avg
cp$B_D = fit_wt$par$B_D
cp$B_max = d["0","signal_avg"] + (d["0","signal_avg"]-cp$B_D)*exp(cp$dG_wt/settings$RT)
cp$B_mid = cp$B_D + (cp$B_max-cp$B_D)/2.0     # Define from curve

print(sprintf("WT fit: dG_wt=%.4f, ddG_avg=%.4f, B_max=%.4f, B_D=%.4f ",cp$dG_wt,cp$ddG_avg,cp$B_max,cp$B_D))

# Take a look at the fit
quartz(height=5, width=9)
plot(d$n_mut, d$signal_avg, ylab="Brightness", ylim=c(cp$B_D,cp$B_max), xlim=c(0,10), xlab="N_mut", col="white",
     main=sprintf("WT fit dG_wt = %.2f kcal/mol, <ddG> = %.2f, B_max=%.2f, B_D=%.2f",cp$dG_wt,cp$ddG_avg,cp$B_max,cp$B_D))
jitter = rnorm(length(mutant$N_sub),0,.1)
points(mutant$N_sub+jitter, mutant$signal, pch=".", col="red")
mask = mutant$init=="use"
points(mutant$N_sub[mask]+jitter[mask], mutant$signal[mask], pch=".")
lines(seq(0,11,.1), stability_wt(fit_wt$par, seq(0,11,.1)), col="red", lwd=2)
points(d$n_mut, d$signal_avg, pch=19, col="green")
points(d$n_mut, d$signal_med, pch=19, col="cyan")
points(d$n_mut, d$frac_act*(wt$signal-cp$B_D)+cp$B_D, pch=19, col="magenta")
abline(h=c(cp$B_max,cp$B_mid,cp$B_D), col="lightgray")
legend("topright", c("Average","Median","Fraction active","used variant","discarede var."), pch=c(19,19,19,46,46), col=c("green","cyan","magenta","black","red"), bg="white")
quartz.save("fit_wt.png", type="png")

################################################################################
# N-mut curve fits of individual substitutions
################################################################################
print("")
print("N-mut curve fit of individual stabilities")

# Fit ddG_tm for a single target substitution
# Fix ddG_avg (typically more difficult to estimate) and rely on midtpoint in ddG estimation. When substitutions are
#   independent, ddG_avg should not change (in practice it is very library dependent but thats a sampling problem)
init = list(ddG_tm=0)
stability = function(param, n_mut) {
    e = exp((cp$dG_wt + param$ddG_tm*(n_mut>0) + cp$ddG_avg * (n_mut-1)*(n_mut>1))/settings$RT)
    (cp$B_max + cp$B_D * e)/(1 + e)
}

residual_stability = function(param, observed, x) {observed - stability(param, x)}

# par(ask=T)
settings$curve_max_deviance = 5.0
settings$curve_max_stderr = 1.0

# Discard hanging node substitutions (gmma=0) because these always have obs=1
# Remember to use a consistant set of substitutions and mutants, e.g. gmma>0
nmut_seq = seq(0,max(mutant$N_sub))
mi_use = which(mutant$init == "use")
for (si in which(subst$init_m == "none")) {
    m = rownames(subst)[si]
    if (subst[si,'obs'] < 10) {
        # print(sprintf("SKIP %s %5d: Too few variant observations %d < 10",m,si,subst[si,'obs'])); subst[si,'m'] = -1
	next
    }

    # Make average fluorescence of data subset contaning substitution si
    mi = intersect(mi_use, subst_mut_indices[[si]])
    # Average is suboptimal because the signal distribution is scwrd towards medium signal, find a smart error model...
    Bavg = data.frame(n_mut = nmut_seq,
                      signal_avg = sapply(nmut_seq, function(n) {mean(mutant[mi[which(mutant[mi,'N_sub']==n)],'signal'])}),
		      n_data = sapply(nmut_seq, function(n) {sum(mutant[mi,'N_sub']==n)}))
    # The reference variant is not returned by subst_mut_indices so add is manually
    Bavg[1,c("signal_avg","n_data")] = c(wt$signal,1)
    # Remove rows for missing multi-mutants
    i = which(is.nan(Bavg$signal_avg))
    if (length(i) > 0) Bavg = Bavg[-i,]

    if (length(Bavg[,1]) < 3) {
        # print(sprintf("SKIP %s %5d:Too few N_mutant points for fitting (%d < 3)",m,si,length(Bavg[,1])))
	subst[si,'init_m'] = "ce_few_Bavg"
    } else if (all(Bavg$signal_avg > cp$B_mid)) {
        # print(sprintf("SKIP %s %5d: Missing deactivated data for",m,si))
	subst[si,'init_m'] = "ce_no_inact"
    } else if (all(Bavg[seq(2,length(Bavg[,1])),'signal_avg'] < cp$B_mid)) { 
        # print(sprintf("SKIP %s %5d: Missing (on average) active N_mutants data. N of first N-mut: %d", m, si, Bavg[2,'n_mut']))
        # IFS requres all mutants to be inactive, here we use average flourescence so this may happen
    	subst[si,'init_m'] = "ce_no_act"
    	if (Bavg[2,'n_mut']==1) {
    	    # if the single mutant is in the data assume deactivationg 
            # subst[si,'init_ddG'] = ddG_deactivating
    	    subst[si,'init_m'] = "ce_no_act"
    	}
    } else {
        # Here, it could make a lot of sence to use weighted least squares to avoid over-emphasis on Bavg means of few data poinst
        fit = nls.lm(par=init, fn=residual_stability, observed=Bavg$signal_avg, x=Bavg$n_mut, control = nls.lm.control(maxiter=200))
        fit_err = sqrt(1/fit$hessian)
        fit_max_res = max(abs(fit$fvec))
        if (fit$info > 3) {
            print(sprintf("SKIP %s %5d: Fitting problem with info=%d",m,si,fit$info))
            print(fit) 
            subst[si,'init_m'] = "ce_fit_info"
        } else if (fit_err > settings$curve_max_stderr | fit$deviance > settings$curve_max_deviance) {
            print(sprintf("SKIP %s %5d: Poor fit with deviance=%.2f or fit error %.2f (from Hessian)",m,si,fit$deviance,fit_err))
            print(fit)
            subst[si,'init_m'] = "ce_high_dev"
        } else {
            subst[si,'init_ddG'] = fit$par$ddG_tm
	    subst[si,'init_m'] = "1_curve"	    
	    print(sprintf("Estimated %s %5d from N-mut curve fit: Fit ddG=%.2f kcal/mol, deviance=%.2f, fit_err=%.2f, fit_max_res=%.2f",m,si,fit$par$ddG_tm,fit$deviance,fit_err,fit_max_res))
	}
        # s = sprintf("%d ddG_%s = %.2f (%.2f) kcal/mol. Obs %d, dev %.2f, max_res %.2f, %s",si,m,subst[si,'init_ddG'],fit_err,subst[si,'obs'],fit$deviance,fit_max_res,subst[si,'init_m'])
        # plot(Bavg$n_mut, Bavg$signal_avg, ylab="Brightness", ylim=c(cp$B_D,cp$B_max), xlim=c(0,10), xlab="N_mut", main=s)
        # abline(h=cp$B_mid, col="lightgray")
        # text(Bavg$n_mut, Bavg$signal_avg, labels=Bavg$n_data, pos=1)
        # lines(c(0,seq(1,10,.1)), stability(list(ddG_tm=subst[si,'init_ddG'], ddG_other=coef(fit)['ddG_other']), c(0,seq(1,10,.1))))
	# # quartz.save("fit_subst.png", type="png")
    }
    # if (subst[si,'m'] != 1) {
    #     plot(Bavg$n_mut, Bavg$signal_avg, ylab="Brightness", ylim=c(1,4), xlim=c(0,10), xlab="N_mut",
    # 	     main=sprintf("%d ddG_%s - No fit",si,m))
    # }
}
print(sprintf("Done Nmut-curve fitting. Estimated %d of %d ddG's",sum(subst$init_m=="1_curve", na.rm=T),nsubst))
print(table(subst$init_m))

################################################################################
# Where curve fits are not possible, estimate if only one ddG is missing in a variant
################################################################################
# Estimate gmma>0 from curve fits only
B_pred = function(par, dG) { e=exp((dG)/settings$RT); (cp$B_max+cp$B_D*e)/(1+e)}
residuals = function(par, observed, ddG_other) {observed - B_pred(par, cp$dG_wt+par$ddG_tm+ddG_other)}
x = seq(-10,10,.1)

settings$init_estimate_excluded = TRUE

# par(ask=T)
# ploti = 0
iterations = 10
prev_est = 0
for (pass in seq(iterations)) {
    # for (si in which(is.na(subst$init_ddG) & ! subst$init_m %in% init_m_excluded_levels)) {
    for (si in which(is.na(subst$init_ddG))) {
        # May have been assigned since the list was generated
        if (!is.na(subst[si,'init_ddG'])) next
	
        # print(sprintf("Subst %4d %5s with %3d obs and m=%d:", si, rownames(subst)[si], subst[si,'obs'], subst[si,'m']))
        # mi = intersect(subst_mut_indices[[si]], mi_use)

	# For each mutant carrying si, how many substitutions are missing curve fitted init_ddG
	if (pass==1) {
	    if (subst[si,'init_m'] %in% init_m_excluded_levels ) next
	    method = "2_curve_ddG1"
	    mi = intersect(subst_mut_indices[[si]], which(mutant$init == "use")) # only use non-excluded data in first pass
	    # missing = sapply(mut_subst_indices[mi], function(l) {sum(subst[l,'m']!=1)})
	    missing = sapply(mut_subst_indices[mi], function(l) {sum(subst[l,'init_m']!="1_curve")})
	} else {
	    if (subst[si,'init_m'] %in% init_m_excluded_levels & ! settings$init_estimate_excluded) next
	    method = "3_curve_ddG2"
	    mi = subst_mut_indices[[si]] # enable estimation of excluded substitutions (ifs and nonsense) by including all data
	    # missing = sapply(mut_subst_indices[mi], function(l) {sum(subst[l,'m']<1)})
	    missing = sapply(mut_subst_indices[mi], function(l) {sum(subst[l,'init_m'] %in% c("none",init_m_err_levels))})
	}
	
	# Data structure for fitting
        dd = data.frame(mi=mi, B=mutant[mi,'signal'], active=mutant[mi,'active'],
                   missing = missing,                
                   ddG_other = sapply(mut_subst_indices[mi], function(l) {sum(subst[l,'init_ddG'], na.rm=T)}), # Sum of non-NA init_ddG's
                   reparam = sapply(mut_subst_indices[mi], function(l) {length(unique(subst[l,'init_reparam']))-1}),
                   row.names=rownames(mutant)[mi])
        # Only keep the rows with one missing substitution which should be si
        d = dd[which((dd$missing-dd$reparam)==1),]
	n_act = sum(d$active >= settings$active_cut)
	n_ina = sum(d$active < settings$inactive_cut)
        # Try to fit
        if (dim(d)[1] > 0) {
            if (n_act > 0 & n_ina > 0) {
                fit = nls.lm(par=list(ddG_tm=0.0), fn=residuals, observed=d$B, ddG_other=d$ddG_other)
	        fit_err = sqrt(1/fit$hessian)
	        fit_max_res = max(abs(fit$fvec))
	        if (fit_err < 1 | fit_max_res < 1) {
                    subst[si,'init_ddG'] = fit$par$ddG_tm
                    subst[si,'init_m'] = method
		    print(sprintf("Estimated dG_%s = %.2f (%.2f) kcal/mol. Max res %.2f, using %d/%d mut and dG-curve fit (m=%s)",
		                  rownames(subst)[si],fit$par$ddG_tm,fit_err,fit_max_res,length(d$mi),length(mi),method))
                } else if (n_act > 10*n_ina  & pass > 2) {
	            # Poor fit may be due to most active except for a few outliers
	            subst[si,'init_ddG'] = -1*max(max(d$ddG_other)+cp$dG_wt+1, 0.0)
		    subst[si,'init_m'] = "5_most_flat"
                    print(sprintf("Estimated dG_%s = %.2f kcal/mol using %d (%d) active (inactive) of %d (m=5)",
		          rownames(subst)[si],subst[si,'init_ddG'],n_act,n_ina,length(mi)))
                } else if (n_ina > 10*n_act & pass > 2) {
	            # Poor fit may be due to most inactive except for a few outliers
	            subst[si,'init_ddG'] = -1*min(min(d$ddG_other)+cp$dG_wt-1, 0.0)
		    subst[si,'init_m'] = "5_most_flat"
		    print(sprintf("Estimated dG_%s = %.2f kcal/mol using %d (%d) inactive (active) of %d (m=5)",
		          rownames(subst)[si],subst[si,'init_ddG'],n_act,n_ina,length(mi)))
                } else if (pass > 3) {
		    # This is typically when there's too few data point to say that something is an outlier
		    # The spread of the average is typically unreasonable so try to avoid this
		    ddG_est = c()
		    for (di in seq(length(d$mi))) {
		        if (d[di,'B'] > .9*cp$B_max) {
			    ddG_est = append(ddG_est, -1*max(d[di,'ddG_other']+cp$dG_wt+1, 0.0))
			} else if (d[di,'B'] <= 1.1*cp$B_D) {
			    # If a binary activity is used as B, B_D may be zero which is why d[di,'B'] should be smaller than _or equal to_ B_D
			    ddG_est = append(ddG_est, -1*min(d[di,'ddG_other']+cp$dG_wt-1, 0.0))
			} else {
		            ddG_est = append(ddG_est, settings$RT*log((cp$B_max-d[di,'B'])/(d[di,'B']-cp$B_D)) - cp$dG_wt - d[di,'ddG_other'])
			}
		    }
		    subst[si,'init_ddG'] = mean(ddG_est)
                    subst[si,'init_m'] = "6_avg_individal"
		    print(ddG_est)
		    print(sprintf("Estimated dG_%s = %.2f (%.2f) kcal/mol using %d/%d individual calculations (m=6_avg_individal)",
		                   rownames(subst)[si],subst[si,'init_ddG'],sd(ddG_est),length(d$mi),length(mi)))
                } else {
		    # only for iterations < 4
                    # subst[si,'m'] = -7
                    subst[si,'init_m'] = "ee_poor_fit"
		}
		
	        # # Plot
	        # dG = d$ddG_other+fit$par$ddG_tm+cp$dG_wt
	        # s = sprintf("dG_%s = %.2f (%.2f) kcal/mol. Max res %.2f, used %d/%d, m=%d",rownames(subst)[si],fit$par$ddG_tm,fit_err,fit_max_res,length(d$mi),length(mi),subst[si,'m'])
	        # plot(dG, d$B, xlim=c(min(dG,-3),max(dG,3)), ylim=c(min(d$B,cp$B_D),max(d$B,cp$B_max)), main=s)
	        # lines(x, B_pred(fit$par,x))
	    # } else if (all(d$B > cp$B_mid) & pass > 2) {
	    } else if (all(d$active >= settings$active_cut) & pass > 2) {
	        # All active cannot fit a curve
		# Consider to use m=6 if there's only 1 data point
	        subst[si,'init_ddG'] = -1*max(max(d$ddG_other)+cp$dG_wt+1, 0.0)
		subst[si,'init_m'] = "4_all_flat"
		print(sprintf("Estimated dG_%s = %.2e kcal/mol using %d/%d active mutants (m=4_all_flat)",rownames(subst)[si],subst[si,'init_ddG'],length(d$mi),length(mi)))
	    # } else if (all(d$B < cp$B_mid) & pass > 2) {
	    } else if (all(d$active < settings$inactive_cut) & pass > 2) {
	        # All inactive cannot fit a curve
		# Consider to use m=6 if there's only 1 data point
	        subst[si,'init_ddG'] = -1*min(min(d$ddG_other)+cp$dG_wt-1, 0.0)
		subst[si,'init_m'] = "4_all_flat"
		print(sprintf("Estimated dG_%s = %.2e kcal/mol using %d/%d inactive mutants (m=4_all_flat)",rownames(subst)[si],subst[si,'init_ddG'],length(d$mi),length(mi)))
		# if (rownames(subst)[si] == "M86S") {print(d$ddG_other)}
            } else {
		# only for iterations < 3
                # subst[si,'m'] = -8
                subst[si,'init_m'] = "ee_bad_data"
            }
        } else {
	    # Not enough data to fit ddG - this should be re-parameterized
            # subst[si,'m'] = -9
            subst[si,'init_m'] = "ee_few_data"
	    if (pass>4) {print(sprintf("=== %s %d estimation error: %s ===",rownames(subst)[si],si,subst[si,'init_m'])); print(dd)}
        }
	# If a ddG is assigned to a subst with reparam then split the ddG amoung the reparam group
	if (!is.na(subst[si,'init_ddG']) & !is.na(subst[si,'init_reparam'])) {
	    group_si = which(subst$init_reparam == subst[si,'init_reparam'])
	    ddG = subst[si,'init_ddG']
	    subst[group_si,'init_ddG'] = ddG/length(group_si)
            method = subst[si,'init_m']
            subst[group_si,'init_m'] = method
	    s=paste(rownames(subst)[group_si], collapse=" ")
	    print(sprintf("Estimated reparam group %d (%s) to %.2f based on %s %d (%.2f, m=%d) ", subst[si,'init_reparam'], s, subst[si,'init_ddG'], rownames(subst)[si], si, ddG, method))
	}
	# # Plot
	# if (!is.na(subst[si,'init_ddG'])) {
	#     ploti = ploti + 1
        #     dG = d$ddG_other + cp$dG_wt + subst[si,'init_ddG']
	#     s = sprintf("Plot %d: ddG_%s = %.2f kcal/mol. used %d/%d, m=%d, gmma=%d",ploti,rownames(subst)[si],subst[si,'init_ddG'],length(d$mi),length(mi),subst[si,'m'],subst[si,'gmma'])
	#     plot(dG, d$B, xlim=c(min(dG,-3),max(dG,3)), ylim=c(min(d$B,cp$B_D),max(d$B,cp$B_max)), main=s)
	#     lines(x, B_pred(list(ddG_tm=subst[si,'m']),x))
	# }
    }
    est = sum(subst$init_m %in% init_m_fit_levels)
    print(sprintf("Estimated init_ddG for %d substitutions after pass %d",est,pass))
    print(table(subst$init_m))
    if (est == prev_est) { print("No more estimation are being made, break"); break } else { prev_est = est}
    print("")
}

################################################################################
# Stability of each mutant
################################################################################
# If some subst$init_ddG are not assigned the NA is propagated
mutant$dG_init = sapply(mut_subst_indices, function(l) {sum(subst[l,'init_ddG'])}) + cp$dG_wt
print(sprintf("Found %d data points used in initial estimations but without estimated dG (due to some missing ddG's)",sum(mutant$init=="use" & is.na(mutant$dG_init))))

B_pred = function(dG) {
    e = exp((dG)/settings$RT)
    return((cp$B_max+cp$B_D*e)/(1+e))
}

dG_pred = function(B) {
    eps = 0.1
    ok_mask = B<(cp$B_max-eps) & B>(cp$B_D+eps)
    ret = rep(NA,length(B))
    ret[which(ok_mask)] = settings$RT*log((cp$B_max-B[which(ok_mask)])/(B[which(ok_mask)]-cp$B_D))
    return(ret)
}

# These are predicted very unstable but measured to beactive
dG_diff = mutant$dG_init - dG_pred(mutant$signal)
B_diff = mutant$signal - B_pred(mutant$dG_init)
mi = which(dG_diff > 2 & B_diff > 1)
print("Active variants with high dG:")
print(mutant[mi,])
print("")

################################################################################
# Dump
################################################################################
# Check if all init_ddG has been assigned to subst with $gmma >0
stopifnot(sum(is.na(subst[which(subst$gmma>0),'init_ddG'])) == 0)

# Put in residue data.frame
residue$init = ""
residue$ifs = ""
init_mask = subst$init_m %in% init_m_fit_levels
ifs_mask = subst$init_m=="ifs" | sapply(subst$note, function(s){"ifs" %in% strsplit(s, " ")[[1]]})
for (ri in seq(nres)) {
    resi_mask = subst$resi==ri
    si = which(init_mask & resi_mask)
    residue[ri,"init"] = paste(subst[si,"taa"], collapse="")
    si = which(ifs_mask & resi_mask)
    residue[ri,"ifs"] = paste(subst[si,"taa"], collapse="")
}

print(sprintf("Estimated %d of %d individual ddG's from %d of %d variants (include all after pass 1: %d) using settings$init_min_obs=%d, settings$init_min_active=%d, settings$ifs_min_active=%d and settings$ifs_min_prob=%.1e",
              sum(! is.na(subst$init_ddG)), dim(subst)[1], sum(mutant$init=="use"), dim(mutant)[1], settings$init_estimate_excluded, settings$init_min_obs, settings$init_min_active, settings$ifs_min_active, settings$ifs_min_prob))
print(sprintf("WT: dG=%.2f, minB=%.1f, maxB=%.1f, ddG_avg=%.2f",cp$dG_wt,cp$B_D,cp$B_max,cp$ddG_avg))

# version=2 is necessary for R versions <3.5 to be able to read it (binf version 3.1)
save(cp, fit_wt, settings, wt, mut_list, mutant, subst, residue, mut_subst_indices, subst_mut_indices, res_mut_indices, file="gmma_fit_init.rda", version=2)
print("Dumped gmma_fit_individual.rda")

