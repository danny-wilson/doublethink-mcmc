########################
# Based on run-mcmc.R  #
########################
# R code to summarize the MCMC using different burn-in or thinning
#
###################

##################
# Read config file #
##################

help = paste(
"Usage: Rscript summarize-again.R config_file taskid burnin thinning",
sep="\n")
# Argument
args = commandArgs(trailingOnly = TRUE)

# Should only have 2 command line arguments
print(args)
if(length(args)!=4) {
	cat(help,sep="\n")
	stop("\nIncorrect usage\n")}

# Get all the variables from the Config file
source(args[1])

# Now check that all the variables have the correct format
taskid = as.numeric(args[2])
if(is.na(taskid)) stop("Nonsensical taskid")
if(taskid<=0) stop("taskid must be positive")
set.seed(taskid)

# Change burnin
old.burnin = burnin
burnin = as.integer(args[3])
cat("Changing burnin from ",old.burnin," to ",burnin,"\n")
stopifnot(burnin>=0)

# Change thinning interval
old.thinning = thinning
thinning = as.integer(args[4])
cat("Changing thinning interval from ",old.thinning," to ",thinning,"\n")
stopifnot(thinning>=1)

##################
# Load functions #
##################
# Get the functions 
sourcefile = paste0(sourcedir,"/doublethink-functions.R")
source(sourcefile)

####################
# Wrangle the data #
####################
# Input file
mcmc.infilename = paste0("doublethink-mcmc-out.ID", taskid ,".RDS")
# Output file
mcmc.summary.outfilename = paste0("doublethink-mcmc-summary-out.ID", taskid ,".burnin",burnin,".thinning",thinning,".RDS")

stopifnot(file.exists(mcmc.infilename))
stopifnot(!file.exists(mcmc.summary.outfilename))

# Load mcmc.result
mcmc.result = readRDS(mcmc.infilename)

# Summarize the MCMC again
mcmc.summary = summarize.mcmc(
	mcmc.result = mcmc.result,
	burnin = burnin,
	thinning = thinning,
	version = 2
)

# Save the summary
saveRDS(object = mcmc.summary, file = mcmc.summary.outfilename)

# Display top hits
tp = data.frame(substr(names(mcmc.summary$mvec.mean),1,50),"%"=round(1000*mcmc.summary$mvec.mean)/10,check.names=F); head(tp[order(-mcmc.summary$mvec.mean),],n=pmin(30,sum(mcmc.summary$mvec.mean>0)))

# Done
