#!/usr/bin/env Rscript
#SBATCH -c 20

##########################
# Reanalyse prior groups #
##########################
# Based on postprocess-mcmc-again.R
# R code to merge MCMC chains and produce final results
# Re-processing using groups defined by prior variable correlation structure
##########################

help = paste(
"Usage: Rscript reanalyse-prior-groups.R config_file",
sep="\n")

# Argument
args = commandArgs(trailingOnly = TRUE)
print(args)
if(length(args)!=1) {
	cat(help,sep="\n")
	stop("\nIncorrect usage\n")
}

# Get the variables from the config file
source(args[1])
setwd(file.path(workdir,outname))

# Load dependencies
library(FMStable, lib.loc=R_lib_loc)
library(harmonicmeanp, lib.loc=R_lib_loc)
# Version of the HMP sampler (1 or 2)
version = 2

# Input/output files: names incorporate updated burning and thinning interval
SUFFIX = paste0(".prior-groups")
infile.a.rsq = paste0(outname,"_corr.csv")
results.one_minus_rsq = paste0(outname,"_one_minus_rsq.tsv")
results.group.posteriors.inclusion.probs = paste0("results.group-posterior-inclusion-probs",SUFFIX,".tsv")

##################
# Load functions #
##################

# Get the folder with the doublethink functions
sourcefile = paste0(sourcedir,"/doublethink-functions.R")
source(sourcefile)

#####################
# Define the groups #
#####################

# Read the squared correlation between variables, and output one minus that quantity
a.rsq = as.matrix(read.csv(infile.a.rsq, row.names=1, check.names=FALSE))
write.table(as.data.frame(1-a.rsq), results.one_minus_rsq, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

# Move current version of results.clusters.out, if it exists
if(file.exists("results.clusters.out")) {
	i = 1
	while(file.exists(paste0("results.clusters.out.",i))) i = i+1
	stopifnot(system(paste0("mv results.clusters.out results.clusters.out.",i))==0)
}

# Cluster based on the prior variable correlation
stopifnot(system("bash -l", input=c("shopt -s expand_aliases", paste0(python_loc, " ", sourcedir, "/cluster_posterior_probabilities.py -cov ", results.one_minus_rsq)))==0)
# Check it worked
stopifnot(file.exists("results.clusters.out"))

# Create factors for groups
gp = as.integer(readLines("results.clusters.out"))
# Force negative values to NA
fgp = factor(gp, levels=0:max(gp))
# Ensure gp is a vector of length ncol(mvec) with integer group
# assignments between 1 and ngp, where ngp is the number of groups
gp = as.integer(fgp)
ngp = max(gp,na.rm=TRUE)
stopifnot(nlevels(factor(gp))==ngp)

# Describe the groups
# Use delimiter \037 (ASCII 31) known as Unit [i.e. field] Separator (US)
who.group = sapply(1:ngp,function(GP) paste0(rownames(a.rsq)[!is.na(gp) & gp==GP],collapse="\037"))
who.group.ct = sapply(1:ngp,function(GP) sum(gp==GP, na.rm=TRUE))

# Create a binary matrix where mgp[i,j] indicates whether variable i
# is a member of group j (1) or not (0).
mgp = 1*outer(gp, 1:ngp, function(x,y) !is.na(x) & x==y)

##############################################
# Reprocess the chains to compute group PoPs #
##############################################
string2mvec = function(s) as.integer(substring(s,1:nchar(s),1:nchar(s)))

chainfiles = sort(dir(".",glob2rx("doublethink-mcmc-out.ID*.RDS")))
i = 1
chain.i = readRDS(chainfiles[i])
string.mvec = with(chain.i, names(loglik.rec)[state])[-(1:(burnin+1))]
mvec = t(sapply(string.mvec, string2mvec))

# NB: mgp might be missing some rows because mvec.gd==0
#     Adjust for that using information in chain.i
mvec.gd = chain.i$mvec.gd

# Compute a group inclusion matrix
# NB: this is a count of the number of variables in each group
#     included in the model in each iteration of the MCMC.
#     However, the focus on whether mvec.gp==0 vs >0.
mvec.gp = mvec[,mvec.gd] %*% mgp
mcmc.PoP.gp = colMeans(mvec.gp>0)
E2.mcmc.PoP.gp = colMeans(mvec.gp>0)^2

# Repeat for the remaining chains
for(i in 2:length(chainfiles)) {
	chain.i = readRDS(chainfiles[i])
	string.mvec = with(chain.i, names(loglik.rec)[state])[-(1:(burnin+1))]
	mvec = t(sapply(string.mvec, string2mvec))
	mvec.gp = mvec[,mvec.gd] %*% mgp
	# mcmc.PoP.gp[j] is the posterior inclusion probability of group j
	# This enables calculation of a standard error for mcmc.PoP.gp across chains
	mcmc.PoP.gp = mcmc.PoP.gp+colMeans(mvec.gp>0)
	E2.mcmc.PoP.gp = E2.mcmc.PoP.gp+colMeans(mvec.gp>0)^2
	print(paste("Reprocessed chain",i))
}
mcmc.PoP.gp = mcmc.PoP.gp/length(chainfiles)
E2.mcmc.PoP.gp = E2.mcmc.PoP.gp/length(chainfiles)
V.mcmc.PoP.gp = E2.mcmc.PoP.gp - mcmc.PoP.gp^2
SE.mcmc.PoP.gp = sqrt(V.mcmc.PoP.gp/length(chainfiles))

names(mcmc.PoP.gp) <- who.group

tb = data.frame("Group"=who.group, "PoP"=mcmc.PoP.gp, "SE"=SE.mcmc.PoP.gp, "df"=who.group.ct)
od = order(tb$PoP, decreasing=TRUE)
write.table(tb[od,], file=results.group.posteriors.inclusion.probs, quote=F, sep="\t")
