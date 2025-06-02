# screen -x 23306
#####################################
# Low memory version of run-mcmc.R  #
#####################################
# R code to run an MCMC on input data
#
# Assume variables are already
# entirely numeric - avoids copying data
#####################################

####################
# Read config file #
####################

help = paste(
"Usage: Rscript run-lowmem.R config_file taskid",
sep="\n")
# Argument
args = commandArgs(trailingOnly = TRUE)

# Should only have 2 command line arguments
print(args)
if(length(args)!=2) {
	cat(help,sep="\n")
	stop("\nIncorrect usage\n")}

# Get all the variables from the Config file
source(args[1])

# For backward compatibility with config files, set default values
if(!exists("seed_subsample")) {
	seed_subsample = 1
	warning("In config, seed_subsample not specified. Setting to default value of 1.")
}
if(!exists("seed_offset")) {
	seed_offset = 0
	warning("In config, seed_offset not specified. Setting to default value of 0.")
}

# Now check that all the variables have the correct format
taskid = as.numeric(args[2])
if(is.na(taskid)) stop("Nonsensical taskid")
if(taskid<=0) stop("taskid must be positive")
set.seed(taskid + seed_offset)

if(is.na(niter)) stop("Nonsensical number of iterations")
if(niter<=0) stop("iterations must be positive")
niter = as.numeric(niter)

chain_numbers = as.numeric(chain_numbers)
if(is.na(chain_numbers)) stop("Nonsensical number of chains")
if(chain_numbers<=0) stop("Number of chains must be positive")


f_prior = as.numeric(f_prior)
if(is.na(f_prior)) stop("Nonsensical f_prior")
if(f_prior<=0) stop("f_prior must be positive")

MU = as.numeric(MU)
if(is.na(MU)) stop("Nonsensical mu")
if(MU<=0 | MU>1) stop("mu must be between 0 and 1")

chain_numbers = as.numeric(chain_numbers)
nlevmax = as.numeric(nlevmax)
burnin = as.numeric(burnin)
thinning = as.numeric(thinning)

if(burnin>=niter) stop("iterations must exceed burnin")

##################
# Load functions #
##################
# Get the functions 
sourcefile = paste0(sourcedir,"/doublethink-functions.R")
source(sourcefile)

####################
# Wrangle the data #
####################
# Output files
mcmc.outfilename = paste0("doublethink-mcmc-out.ID", taskid ,".RDS")
mcmc.summary.outfilename = paste0("doublethink-mcmc-summary-out.ID", taskid ,".RDS")
traceplot.outfilename = paste0("doublethink-mcmc-trace-out.ID", taskid,".png")
system(paste("cp", exclude_columns, "." ))

stopifnot(!file.exists(mcmc.summary.outfilename))
stopifnot(!file.exists(traceplot.outfilename))

# Identify the input files
load.data = function(columns_filename, input_filename, cases_filename, exclude_filename, seed_subsample=1, seed_offset=0) {
	# Define the input column types for faster reading
	col = rbind(cbind("row","character"),read.csv(columns_filename,h=F,stringsAsFactors=F),stringsAsFactors=F); colnames(col) = c("column","class")

	# Assert that all input columns are already doubles
	# Consequently, some of what follows will be redundant
	stopifnot(all(col$class[-1]=="double"))

	# Summarize column types
	table(col$class)

	# R doesn't have class float, change to double
	col$class[col$class=="float"] = "double"

	# Read the data
	a = read.csv(input_filename,check.names=F,colClasses=as.character(col$class))

	# Get the ids of participants to remove and remove them
	exclude.eids = scan(exclude_filename)
	a=a[!a$eid %in% exclude.eids,]

	# Factors: replace any negative factor level with -999: next line may issue warning
	# In FUN(if (length(d.call) < 2L) newX[, 1] else array(newX[, 1L],  :
	#   NAs introduced by coercion"
	wh = which(apply(a[,col$class=="factor",drop=FALSE],2,function(y) as.numeric(as.character(y)))<0,arr=TRUE)
	a[,col$class=="factor"][wh] = "-999"

	# Check there are no missing factor values
	stopifnot(sum(is.na(a[,col$class=="factor",drop=FALSE]))==0)

	# Drop redundant levels in factors (i.e. the negative levels other than -999)
	for(j in which(col$class=="factor")) a[,j] = droplevels(a[,j])

	# Covariates: crudely impute missing values with the mean
	for(j in which(col$class!="factor")) {
		IS.NA = is.na(a[,j])
		if(sum(IS.NA)>0) a[IS.NA,j] = mean(a[,j],na.rm=TRUE)
	}

	# Load the phenotypes as a table of IDs and pheno
	case.eids = read.csv(cases_filename)

	# Define the phenotype matrix
	case.eids = case.eids[!case.eids$eid %in% exclude.eids,]
	pheno = data.frame(case.eids[match(a[,"eid"], case.eids[,"eid"]),])
	if(controls_subsample !=0){
		controls = pheno[pheno[,"pheno"] == 0,]
		cases = pheno[pheno[,"pheno"] == 1,]
		# Temporarily fix the seed at a pre-specified level for consistently
		# sub-sampling the controls across chains. Then restore the state of the RNG
		saved.seed = .Random.seed
		set.seed(seed_subsample)
		subsampled_eids = sample(controls$eid, controls_subsample, replace=FALSE)
		write(subsampled_eids, paste0("subsampled_control_eids.", taskid, ".txt"), ncol=1)
		.Random.seed = saved.seed
		all_eids = c(cases$eid, subsampled_eids)
		a=a[a$eid %in% all_eids,]
		pheno=pheno[pheno$eid %in% all_eids,]
	}

	# Check eids match between pheno and a
	stopifnot(all(a$eid==pheno$eid))
	
	# Now drop them from a
	rownames(a) = a$eid
	a$eid = NULL
	
	# Also drop first (unnamed) column of row numbers
	a[,match("",colnames(a))] = NULL

	# Type cast object a from data frame to matrix since all entries are doubles
	ret = list("pheno"=pheno, "a"=as.matrix(a))
}

# Apply the load function
data = load.data(
	columns_filename = columns_filename,
	input_filename = input_filename,
	cases_filename = cases_filename,
	exclude_filename = exclude_filename,
	seed_subsample = seed_subsample,
	seed_offset = seed_offset
)

# Get the number of levels from factor
data$nlev = unlist(lapply(data$a,nlevels))

# Table of the number of levels of each column of the data
table(data$nlev)

# Apply the maximum number of levels
data$nlevmax = nlevmax


##################################################################
# No longer: Create matrix b, a purely numeric counterpart to a  #
# Instead scale a in situ.                                       #
##################################################################

print("All chains now create standardised matrix")

# Create standardised matrix in situ
{
	colnames.remove = c()
	mtc = match(colnames.remove,colnames(data$a))
	stopifnot(!any(is.na(mtc)))
	col.use = is.na(match(colnames(data$a),colnames.remove))
	for(j in which(col.use)) {
		colname = colnames(data$a)[j]
		
		IS.NA = is.na(data$a[,j])
		# Crude mean-value imputation for covariates
		data$a[IS.NA,j] = mean(data$a[!IS.NA,j])
		data$a[,j] = rescale(data$a[,j,drop=FALSE])

		if(any(!is.finite(data$a[,j]))) {
			col.use[j] = FALSE
			# We need one chain to standardise the matrix so all chains can have the same excluded columns otherwise we have trouble in postprocessing
			if(taskid==1) {
				warning("Excluding column ", j, " with zero variance: ", colname)
				write(colname,file=exclude_columns,append=TRUE)
				system(paste0("sort -u ", exclude_columns, " > tmp_exclude.txt" ) )
				system(paste0("cat tmp_exclude.txt > ", exclude_columns))
			}
		}
		if(j%%10==0) print(paste("Done",j,"of",ncol(data$a)))
	}
}
	
print("NA a")
print(colnames(data$a)[which(is.na(data$a))])

# Assert there are no longer any NAs
stopifnot(!any(is.na(data$a)))

###################################################
# Calculate correlations between the columns of a #
###################################################

# Output the correlations so they can be used to create chain initialisation groups using make_starting_groups.py
corr_file = paste0(outname, "_corr.csv")
start_group_file = paste0(outname,"_groups.csv")
# This file will stop all other chains from also creating a correlation and starting group file. We only need one 
temp_file = "hold_file1"

# No longer: Check if the temporary file exists. Give priority to the first chain instead
if(taskid==1) {

	# If it doesn't then make correlation and starting group file
	print("This is the leading chain, will create correlation and chain initialisation file")
	system(paste0("touch ", temp_file))

	# Create correlation matrix
	stopifnot(!any(duplicated(colnames(data$a))))
	# Compute the correlation structure of a (slow)
	print("Calculating correlation matrix")
	data$a.rsq = cor(data$a)^2
	
	print(colnames(data$a.rsq)[which(is.na(data$a.rsq))])

	# Get column names in
	colnames(data$a.rsq) = colnames(data$a)
	rownames(data$a.rsq) = colnames(data$a)

	# Write output
	write.csv(data$a.rsq, corr_file)

	# Now run the chain initialisation with furthest neighbour in python
	start_length = as.integer(1/MU)
	chain_init_cmd = paste0(python_loc, " ", sourcedir, "/make_starting_groups.py -n ", outname, " -m ", chain_numbers," -l ", start_length)
	system(chain_init_cmd)
} else {  # If it does exist wait for all the necessary files to be created

	print("Waiting for the leading chain to create chain initialisation file")
	
	# Wait for the start group file to be created
	while(!file.exists(start_group_file)) {
		  Sys.sleep(5)
	}
	
	print("Done waiting")
	data$a.rsq = read.csv(corr_file, row.names = "X")
}

############################################
# Model selection Markov chain Monte Carlo #
############################################

# Read the columns that are supposed to be excluded from the analysis
x = readLines(exclude_columns)

# No empty lines
x[x != ""]

# Only exclude if they're truly in the dataset
x = x[(x %in% colnames(data$a))]

# Store them in remove
data$colnames.remove = x

# Start the chain
mcmc.result = lr.mcmc(
	pheno = data$pheno$pheno,
	a = data$a,
	colnames.remove = data$colnames.remove,
	MU = MU,
	chainid = taskid,
	outname = outname,
	f_prior = f_prior,
	niter = niter,
	report.move = FALSE,
	report.state = TRUE,
	a.rsq = data$a.rsq,
	nlevmax = data$nlevmax,
)

# Save the results
saveRDS(object = mcmc.result, file = mcmc.outfilename)

# Summarize the MCMC
mcmc.summary = summarize.mcmc(
	mcmc.result = mcmc.result,
	burnin = burnin,
	thinning = thinning,
	version = 2
)

# Save the summary
saveRDS(object = mcmc.summary, file = mcmc.summary.outfilename)

# Trace plot of the log-likelihood or log-posterior
trace.plot(mcmc.result, traceplot.outfilename)

# Display top hits
tp = data.frame(substr(names(mcmc.summary$mvec.mean),1,50),"%"=round(1000*mcmc.summary$mvec.mean)/10,check.names=F); head(tp[order(-mcmc.summary$mvec.mean),],n=pmin(30,sum(mcmc.summary$mvec.mean>0)))

# Done
