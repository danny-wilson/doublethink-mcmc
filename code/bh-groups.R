#!/usr/bin/env Rscript

#####################
# Compute BH groups #
#####################

# Based on postprocess-mcmc.R

help = paste(
"Usage: Rscript bh-groups.R config_file",
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

infile.posterior.inclusion.probs = "results.posterior-inclusion-probs.tsv"

results.table.bh.groups = "results.table.bh.groups.tsv"

# Specify the posterior inclusion probability threshold for analysis
pp.thresh = 0.01

##################
# Load functions #
##################

# Get the folder with the doublethink functions
sourcefile = paste0(sourcedir,"/doublethink-functions.R")
source(sourcefile)

###########################
# Load MCMC summary files #
###########################

# Locate the MCMC summary files
summaryfiles = dir(".",glob2rx("doublethink-mcmc-summary-out.ID*.RDS"))
od = order(as.numeric(gsub("ID","",sapply(summaryfiles,function(s)unlist(strsplit(s,".",fixed=TRUE))[2]))))
if(any(is.na(od))) stop("Non-sensical task IDs in output workdir")
summaryfiles = summaryfiles[od]
n.M = length(summaryfiles)
print(paste("Found",n.M,"summary files"))

# Load them
rds = list()
for(i in 1:n.M) rds[[i]] = readRDS(summaryfiles[i])
print(paste("Read",n.M,"summary files"))

# Load one chain directly to obtain parameter details
chainfiles = gsub("summary-","",summaryfiles)
if(!file.exists(chainfiles[1])) stop("Could not find chain file")
system.time((chain = readRDS(chainfiles[1])))

# In case need to load all the chains directly:
#if(!any(file.exists(chainfiles))) stop("Could not find chain files")

##################################
# Check consistency of the files #
##################################

# Check if one of the chains was prematurely ended
L = length(rds[[1]]$mvec.mean)
regressor.names = names(rds[[1]]$mvec.mean)
if(n.M>1) for(i in 2:n.M) {
	if(length(rds[[i]]$mvec.mean)!=L) stop("Inconsistent number of variable lengths across chains")
	if(!all(names(rds[[i]]$mvec.mean)==regressor.names)) stop("Inconsistent variable names across chains")
}
print(paste("Found",L,"regressors, beginning:"))
print(head(regressor.names))

##################################
# Load posterior inclusion probs #
##################################

posterior.inclusion.probs = read.delim(infile.posterior.inclusion.probs, check.names=FALSE)
if(!all(rownames(posterior.inclusion.probs)==regressor.names)) stop("Inconsistent variable names in posterior inclusion probs file")

# Define an ordering of variables
od = order(posterior.inclusion.probs$PoP, decreasing=TRUE)
# And the number with posterior inclusion probability above pp.thresh
ntop = sum(posterior.inclusion.probs$PoP>=pp.thresh) # 67

# Create a binary matrix where mgp[i,j] indicates whether variable i
# is a member of group j (1) or not (0).
# Each includes the bottom however-many variables
mgp = matrix(1, length(regressor.names), ntop)
for(i in 1:ntop) mgp[od[1:i], i] = 0

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


mcmc.PoP.gp
posterior.inclusion.probs$PoP[od][1:ntop]

# Posterior distribution of the number of variables for inclusion
mean(sapply(chain.i$state[burnin:niter], function(i) chain.i$loglik.rec.n[[i]]))
quantile(sapply(chain.i$state[burnin:niter], function(i) chain.i$loglik.rec.n[[i]]), c(0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95))
# 5% 10% 20% 50% 80% 90% 95%
# 18  19  19  20  21  22  22
nlo = quantile(sapply(chain.i$state[burnin:niter], function(i) chain.i$loglik.rec.n[[i]]), c(0.1)) # 19

# Cumulative FDR
cummean = function(x) cumsum(x)/1:length(x)
fdr = cummean(1-posterior.inclusion.probs$PoP[od])

# Target, say tau=10
tau = 10; fdr.thresh = 1/(1+tau) # 0.09090909

# Size of the BH group
which(fdr > fdr.thresh)[1]-1 # 14 (cf 19 for BH based on the p-values)

# FDR based on the lower 90% confidence interval on number of variables for inclusion
fdr[nlo] # 0.1813351

# The top nlo
regressor.names[od][1:nlo]

# Sex only just makes it, age does not: limitations of analysing individual variables
# whether with local FDR, BH, or whatever
#[1] "41214 Carer support indicators : 1 : Yes"
#[2] "Z86.4 Personal history of psychoactive substance abuse"
#[3] "F03 Unspecified dementia"
#[4] "137 Number of treatments/medications taken"
#[5] "J22 Unspecified acute lower respiratory infection"
#[6] "R29.6 Tendency to fall, not elsewhere classified"
#[7] "41218 History of psychiatric care on admission : 8 : Not applicable"
#[8] "6138 Qualifications : 3 : O levels/GCSEs or equivalent"
#[9] "6138 Qualifications : 1 : College or University degree"
#[10] "48 Waist circumference (cm)"
#[11] "J18.1 Lobar pneumonia, unspecified"
#[12] "Z50.1 Other physical therapy"
#[13] "26413 Health score (England)"
#[14] "21000 Ethnic background : 1001 : British"
#[15] "K59.0 Constipation"
#[16] "N39.0 Urinary tract infection, site not specified"
#[17] "31 Sex : 0 : Female"
#[18] "31 Sex : 1 : Male"
#[19] "34 Year of birth (years)"
