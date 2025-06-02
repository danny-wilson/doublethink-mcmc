# screen -x 15143
# module add R/3.6.2-foss-2019b
# R --args /filepath/B-2023-05-26-I/B-2023-05-26-I.cfg

############
# Calc FDP #
############
# Based on calc-fdp.R
# R code to merge MCMC chains and produce final FDP
#
###################

help = paste(
"Usage: Rscript calc-fdp.R config_file",
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

# Load dependencies
library(FMStable, lib.loc=R_lib_loc)                           
library(harmonicmeanp, lib.loc=R_lib_loc) 

# Version of the HMP sampler (1 or 2)
version = 2

# Output files
results.fdp = "results.fdp.tsv"

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

#~~~~~~~~~~#~~~~~~~~~~#~~~~~~~~~~#~~~~~~~~~~#~~~~~~~~~~#~~~~~~~~~~#~~~~~~~~~~#~~~~~~~~~~

n.s = matrix(0, n.M, length(rds[[1]]$beta.sim[[1]]))
for(i in 1:n.M) {
	for(j in 1:length(rds[[1]]$beta.sim)) {
		n.s[i,] = n.s[i,] + 1*(rds[[i]]$beta.sim[[j]]!=0)
	}
	cat("Done", i, "\n")
}
table(n.s)
#n.s
#	 16      17      18      19      20      21      22      23      24      25
#	465   50769  521532 1719547 2469244 1785084  725209  192133   31897    3715
#	 26      27      28
#	339      62       4

quantile(n.s, c(0.025, 0.05, 0.1, 0.2))
#2.5%   5%  10%  20%
#  18   18   19   19

# Compare to BH approach
#pp =
#od =
