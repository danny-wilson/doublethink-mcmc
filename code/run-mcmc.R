#########################
# NICK MODIFIED VERSION #
#########################
# R code to run an MCMC on input data
#
###################

##################
# Read config file #
##################

help = paste(
"Usage: Rscript run-mcmc.R config_file taskid",
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
if(!exists("flat_likelihood")) {
	flat_likelihood = FALSE
	warning('In config, flat_likelihood not specified. Setting to default value of "FALSE".')
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


h_prior = as.numeric(h_prior)
if(is.na(h_prior)) stop("Nonsensical h_prior")
if(h_prior<=0) stop("h_prior must be positive")

mu_prior = as.numeric(mu_prior)
if(is.na(mu_prior)) stop("Nonsensical mu_prior")
if(mu_prior<=0) stop("mu_prior must be positive")

chain_numbers = as.numeric(chain_numbers)
nlevmax = as.numeric(nlevmax)
burnin = as.numeric(burnin)
thinning = as.numeric(thinning)

if(burnin>=niter) stop("iterations must exceed burnin")

if(is.character(flat_likelihood)) flat_likelihood = tolower(flat_likelihood)=="true"

# Obtain the ID column name
if(!exists("idcol")) idcol = "eid"

# Obtain the phenotype column name
if(!exists("phenocol")) phenocol = "pheno"

# Obtain the glm family
if(!exists("glm_family")) glm_family = "binomial"

# Set an intercept?
if(!exists("intercept")) {
	intercept = TRUE
} else {
	intercept = toupper(intercept)=="TRUE"
}

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

	# Summarize column types
	table(col$class)

	# R doesn't have class float, change to double
	col$class[col$class=="float"] = "double"

	# Read the data
	a = read.csv(input_filename,check.names=F,colClasses=as.character(col$class))

	# Get the ids of participants to remove and remove them
	exclude.eids = scan(exclude_filename)
	a=a[!a[,idcol] %in% exclude.eids,]

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
	case.eids = case.eids[!case.eids[,idcol] %in% exclude.eids,]
	pheno = data.frame(case.eids[match(a[,idcol], case.eids[,idcol]),])
	if(!exists("controls_subsample")) controls_subsample = 0
	if(controls_subsample !=0){
		controls = pheno[pheno[,phenocol] == 0,]
		cases = pheno[pheno[,phenocol] == 1,]
		# Temporarily fix the seed at a pre-specified level for consistently
		# sub-sampling the controls across chains. Then restore the state of the RNG
		saved.seed = get(".Random.seed", envir=.GlobalEnv)
		set.seed(seed_subsample)
		subsampled_eids = sample(controls[,idcol], controls_subsample, replace=FALSE)
		write(subsampled_eids, paste0("subsampled_control_eids.", taskid, ".txt"), ncol=1)
		assign(".Random.seed", saved.seed, envir=.GlobalEnv)
		all_eids = c(cases$eid, subsampled_eids)
		a=a[a[,idcol] %in% all_eids,]
		pheno=pheno[pheno[,idcol] %in% all_eids,]
	}

	# Check eids match between pheno and a
	stopifnot(all(a[,idcol]==pheno[,idcol]))
	
	# Now drop them from a
	rownames(a) = a[,idcol]
	a[,idcol] = NULL
	
	# Also drop first (unnamed) column of row numbers
	a[,match("",colnames(a))] = NULL

	ret = list("pheno"=pheno, "a"=a)
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


#######################################################
# Create matrix b, a purely numeric counterpart to a  #
# This matrix is scaled so it can be encoded as dummy #
# genotypes if necessary for analysis in SAIGE etc    #
#######################################################

# First column: eids
temp_file = "hold_file"
done_file = "done_file"
matrix_file = "standardised_matrix.csv"
# We need one chain to standardise the matrix so all chains can have the same excluded columns otherwise we have trouble in postprocessing
if(taskid==1) {
    stopifnot(!file.exists(temp_file))
	system(paste("touch", temp_file))
	print("This is the leading chain, will create standardised matrix")

	# Create standardised matrix
	data$b = create.standardized.matrix(data$a, exclude_columns, impute.NA=TRUE, nlevmax=data$nlevmax)
	write.csv(data$b, matrix_file)

	system(paste("touch", done_file))

} else {  # If it does exist wait for all the necessary files to be created
	while (!file.exists(done_file)) {
		  Sys.sleep(5)
	}
	print("Waiting for the leading chain to create standardised matrix")
	data$b = read.csv(matrix_file, row.names = "X")
}
print("NA a")
print(colnames(data$a)[which(is.na(data$a))])
print("NA b")
print(colnames(data$b)[which(is.na(data$b))])


# Assert there are no longer any NAs
stopifnot(!any(is.na(data$a)))


#######################################################################
# Read the columns that are supposed to be excluded from the analysis #
#######################################################################
x = readLines(exclude_columns)

# No empty lines
x[x != ""]

# Only exclude if they're truly in the dataset
x = x[(x %in% colnames(data$a))]

# Store them in remove
data$colnames.remove = x

#############################################################
# Convert h_prior and mu_prior to internal parameterization #
#############################################################
data_a_nlev = unlist(lapply(data$a,nlevels))
data_a_colnums_remove = match(data$colnames.remove, colnames(data$a))
stopifnot(!any(is.na(data_a_colnums_remove)))
data_a_nmax = length(setdiff(which(data_a_nlev != 1 & data_a_nlev <= data$nlevmax), data_a_colnums_remove))
MU = 1.0 / ( mu_prior/(1.0 + mu_prior) * data_a_nmax )
f_prior = 1.0 / h_prior

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
	data$a.rsq = create.regressor.correlation.matrix(a=data$a, b=data$b, colnames.remove=setdiff(colnames(data$a), colnames(data$b)))
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

# Start the chain
mcmc.result = lr.mcmc(
	pheno = data$pheno[,phenocol],
	a = data$a,
	colnames.remove = data$colnames.remove,
	MU = MU,
	chainid = taskid,
	outname = outname,
	f.prior = f_prior,
	niter = niter,
	report.move = FALSE,
	report.state = TRUE,
	a.rsq = data$a.rsq,
	nlevmax = data$nlevmax,
	flat_likelihood = flat_likelihood,
	glm_family = glm_family,
	intercept = intercept
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
