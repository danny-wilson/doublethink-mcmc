# Packages are loaded from R_lib_loc. The variable is inherited from the parent environment
library(FMStable, lib.loc=R_lib_loc)                           
library(harmonicmeanp, lib.loc=R_lib_loc) 

############################################
#         Error handling function          #
############################################

# Courtesy of https://renkun.me/2020/03/31/a-simple-way-to-show-stack-trace-on-error-in-r/
if(!interactive()) options(error = function() {
  calls <- sys.calls()
  if (length(calls) >= 2L) {
    sink(stderr())
    on.exit(sink(NULL))
    cat("Backtrace:\n")
    calls <- rev(calls[-length(calls)])
    for (i in seq_along(calls)) {
      cat(i, ": ", deparse(calls[[i]], nlines = 1L), "\n", sep = "")
    }
  }
  coredump = paste0("coredump.",Sys.getpid(),".RData")
  if(coredump!="") save.image(coredump)
  if (!interactive()) {
    q(status = 1)
  }
})

############################################
#            Plotting functions            #
############################################

# Plotting function for the tracing plots
pl = function(expr,filename,height=6,width=10,units="in",res=300,...) {
	bitmap(filename,height=height,width=width,units=units,res=res,...)
	tryCatch(expr,error=function(e) cat(e$message,"\n"))
	dev.off()
	invisible()
}

# Trace plot of the log-likelihood or log-posterior
trace.plot = function(mcmc.ret, outfilename) {
	invisible(pl({
		par(mfrow=c(2,1))
		# Temporarily record the (unnormalized) log-posterior probability of each model over time
		tp = sapply(mcmc.ret$state[!is.na(mcmc.ret$state)],function(i)mcmc.ret$loglik.rec[[i]]+mcmc.ret$logPrP[[i]])
		# Plot it on a log scale
		plot(tp,xlab="Iteration",ylab="Unnormalized log-posterior model probability")
		# Normalize the log-posterior model probabilities
		denom = logsum(tp)
		# Plot them on an unlogged scale
		plot(exp(tp-denom),xlab="Iteration",ylab="Posterior model probability")
	},filename = outfilename))
}

############################################
#        Simple stats functions            #
############################################

# Simes' test (analogous to harmonic mean p-value test)
p.Simes = function(p, w = NULL, L = NULL, w.sum.tolerance = 1e-6, multilevel = TRUE) {
  if(is.null(L) & multilevel) {
    warning("L not specified: for multilevel testing set L to the total number of individual p-values")
    L = length(p)
  }
  if(length(p) == 0) return(NA)
  if(length(p) > L) stop("The number of p-values cannot exceed L")
  if(is.null(w)) {
    w = rep(1/L,length(p))
  } else {
    if(any(w<0)) stop("No weights can be negative")
    if(length(w)!=length(p)) stop("When specified, length of w must equal length of p")
  }
  w.sum = sum(w)
  if(w.sum>1+w.sum.tolerance) {
    stop("Weights cannot exceed 1")
  }
  pstar = p/w
  return(c(p.Simes = min(sort(pstar)/(1:length(p)))))
}


# Convert a factor vector to a matrix of dummy variables
factor2dummy = function(ftr) {
  LEV = levels(ftr)
  NLEV = nlevels(ftr)
  x = ftr==LEV[1]
  if(NLEV>2) {
	for(j in 2:(NLEV-1)) {
	  x = cbind(x,ftr==LEV[j])
	}
  }
  return(1*x)
}
dp = function(x, dp=1) round(10^dp*x)/10^dp

# probility to odds function
p2o = function(p) p/(1-p)

# odds to probability function
o2p = function(o) o/(1+o)

# sum of logarithms function 
logsum = function(x,na.rm=FALSE) {
	mx = max(x,na.rm=na.rm)
	log(sum(exp(x-mx),na.rm=na.rm))+mx
}

# mean of logarithms function
logmean = function(x,na.rm=FALSE) {
	mx = max(x,na.rm=na.rm)
	log(mean(exp(x-mx),na.rm=na.rm))+mx
}

# Function to compute the p-value for the nested pair of models excluding/including j
calc.neglogp.pair = function(loglik.0,loglik.A) as.numeric(-pchisq(2*(loglik.A-loglik.0),attributes(loglik.A)$df-attributes(loglik.0)$df,low=FALSE,log=TRUE))

# Rescale a dummy variable matrix to have summed variance 1
# ignoring covariance between the columns of the matrix
rescale = function(mat) {
	n = nrow(mat)
	mat = t(t(mat)-colMeans(mat))
	v = (n-1)/n*sum(apply(mat,2,var))
	return(mat/sqrt(v))
}

# Create a standardized numeric matrix for covariates and factors
create.standardized.matrix = function(a, exclude_columns, colnames.remove=c(), impute.NA=FALSE, nlevmax=50) {
	mtc = match(colnames.remove,colnames(a))
	stopifnot(!any(is.na(mtc)))
	col.use = which(is.na(match(colnames(a),colnames.remove)))
	nlev = unlist(lapply(a,nlevels))
	b = matrix(0.0,nrow=nrow(a),ncol=0)
	dummy_threshold = round(nrow(a) * 0.001, digits=0)
	print(paste0("Dummy threshold ", dummy_threshold))
	for(j in col.use) {
		colname = colnames(a)[j]
		skip = FALSE
		if(nlev[j]==0) {
			b.elem = as.matrix(a[,j,drop=FALSE])
			IS.NA = is.na(b.elem)
			if(any(IS.NA)) stopifnot(impute.NA)
			# Crude mean-value imputation for covariates
			b.elem[IS.NA] = mean(b.elem[!IS.NA])
			b.elem = rescale(b.elem)
			colnames(b.elem) = colname
		} else {
			if(nlev[j]<=nlevmax) {
				b.elem = sapply(levels(a[,j]),function(LEV) a[,j,drop=FALSE]==LEV)
				IS.NA = is.na(b.elem)
				if(any(IS.NA)) stopifnot(impute.NA)
				# For factors, crudely impute as false for all levels
				b.elem[is.na(b.elem)] = FALSE
				b.elem = rescale(1*b.elem)
				colnames(b.elem) = paste0(colname,":",levels(a[,j]))
			} else skip = TRUE
		}
		if(any(!is.finite(b.elem))) {
			warning("Excluding column ", j, " with zero variance: ", colname)
			write(colname,file=exclude_columns,append=TRUE)
			system(paste0("sort -u ", exclude_columns, " > tmp_exclude.txt" ) )
			system(paste0("cat tmp_exclude.txt > ", exclude_columns))
			skip = TRUE
		} else if(startsWith(colname, "dummy")) {
			tab = table(a[,j])
			cases_count = tab[names(tab)==1]
			if(cases_count<dummy_threshold) {
				print(paste0("Excluding dummy column ", j, " because of less than 0.1% positive cases: ", colname, " has only ", cases_count, " positive labels"))
				write(colname,file=exclude_columns,append=TRUE)
				system(paste0("sort -u ", exclude_columns, " > tmp_exclude.txt" ) )
				system(paste0("cat tmp_exclude.txt > ", exclude_columns))
				skip = TRUE
		} }
		if(!skip) {
			b = cbind(b,b.elem)
		}
		if(j%%10==0) print(paste("Done",j,"of",ncol(a)))
	}
	return(b)
}

# Use matrix b to calculate correlations between the columns of a #
create.regressor.correlation.matrix = function(a, b=NULL, colnames.remove=c(), impute.NA=FALSE, nlevmax=50) {
	stopifnot(!any(duplicated(colnames(a))))
	if(is.null(b)) {
		print("Creating standardized matrix")
		b = create.standardized.matrix(a=a, colnames.remove=colnames.remove, impute.NA=impute.NA, nlevmax=nlevmax)
		colnames.removed = setdiff(colnames(a), colnames(b))
		# Make sure create.standardized.matrix removed the specified columns
		stopifnot(length(setdiff(colnames.remove,colnames.removed))==0)
		# Add to the vector any other columns removed by create.standardized.matrix
		colnames.remove = colnames.removed
	}
	# Which columns of a to use
	col.use = which(is.na(match(colnames(a),colnames.remove)))
	# Compute the correlation structure of b (slow)
	print("Calculating correlation matrix")
	b.cor = cor(b)
	# Map the column names of a to the column names of b
	a2b = sapply(colnames(a),grep,colnames(b),fixed=TRUE)
	stopifnot(!any(sapply(col.use,function(i) length(a2b[[i]])==0)))
	# Compute the squared correlation structure of a from b
	a.rsq = matrix(1.0,ncol(a),ncol(a))
	for(i in 2:ncol(a)) {
		for(j in 1:(i-1)) {
			tp = b.cor[a2b[[i]],a2b[[j]],drop=FALSE]
			if(length(tp)>0) {
				# For factors, take the mean squared correlation
				a.rsq[i,j] = mean(tp^2)
				a.rsq[j,i] = mean(tp^2)
			} else {
				a.rsq[i,j] = 0.0
				a.rsq[j,i] = 0.0
			}
		}
		if(i%%50==0) print(paste("Done",i,"of",length(col.use)))
	}
	return(a.rsq)
}

############################################
# Model selection Markov chain Monte Carlo #
############################################

# pheno				vector of phenotypes (binary)
# a					data frame of regressors, # rows equal to length of pheno
# colnames.remove	names of columns of a to ignore
# MU				prior parameter, inverse expected number of regressors
# f.prior			g-prior parameter, e.g. 1 (AIC) or 1/length(pheno) (BIC)
# niter				number of iterations of the MCMC chain
# report.move		flag: output MCMC proposals to screen?
# report.state		flag: output MCMC state to screen?
# a.rsq				pre-computed correlation matrix between regressors

lr.mcmc = function(pheno, a, colnames.remove, MU, chainid, outname, f.prior=1, niter, report.move=FALSE, report.state=TRUE, a.rsq=NULL, nlevmax=50, flat_likelihood=FALSE, glm_family="binomial", intercept=TRUE) {
	stopifnot(nrow(a)==length(pheno))
# Parameters: 1/MU is the expected number of regressors in the model
#             assuming a geometric prior
#	MU = MU
#             samp.size is the sample size
	samp.size = length(pheno)
#             f.prior is the expected magnitude of effect sizes
#             on a scale relative to the per-observation precision
#             defined via the Fisher information matrix
#             Standard magnitudes for f.prior are 1 (implies interest in
#             effect sizes that should be readily identified in large
#             samples and an AIC-like criterion for model selection)
#             and 1/samp.size (implies interest in effect sizes that are
#             weakly identifiable given the sample size and a BIC-like
#             criterion for model selection)
#	f.prior = f.prior
#             xi.prior is a compound parameter depending on samp.size
#             and f.prior
	xi.prior = 1/(1+samp.size*f.prior)
#             number of iterations of the MCMC
#	niter = niter
#             whether to output proposal details to screen
#	report.move = report.move
#             whether to output MCMC state to screen
#	report.state = report.state

# Data: define the regressors allowed in the model
#       Exclude uninformative (single-level) factors, factors with
#       too many levels, column 1, "eid" and the noise columns 522:533
#	pheno = pheno

	# initialise a prediction table where we store all the predictions of the logistic regression
	prediction_table = data.frame(matrix(,nrow=0 ,ncol=nrow(a)))

	# This is where the file will be written into
	prediction_file = paste0(workdir,"/",outname,"/doublethink-prediction.ID", chainid, ".csv") 

	if(is.null(a.rsq)) a.rsq = create.regressor.correlation.matrix(a=a, colnames.remove=colnames.remove, impute.NA=TRUE, nlevmax=nlevmax)
#	colnames.remove = colnames.remove
	nlev = unlist(lapply(a,nlevels))
	colnums.remove = match(colnames.remove,colnames(a))
	stopifnot(!any(is.na(colnums.remove)))
	mvec.gd = setdiff(which(nlev!=1 & nlev<=nlevmax),colnums.remove)
#       Maximum number of regressors in the model, needed for the prior
# Functions: log-prior model probability
	Nmax = length(mvec.gd)
	logpr.prior = function(n,mu=MU,nmax=length(mvec.gd)) dbinom(n,nmax,1/mu/nmax,log=TRUE)-lchoose(nmax,n) # binomial distribution
	base_model = "pheno ~ 1"
	if(!intercept) base_model = "pheno ~ 0"
	gen.calc.loglik = function(flat_likelihood=FALSE) {
	#        log-likelihood of the null model is used as a baseline
		fit0.logLik = logLik(glm(as.formula(base_model), family=glm_family, singular.ok=FALSE))
	#        mvec is a binary encoding of the regressors to include
		function(mvec) {
		#    which columns of a act as regressors
			wh = which(mvec==1)
		#    sanity check
			stopifnot(!any(nlev[wh]>nlevmax))
			if(flat_likelihood==TRUE) return(list("loglik"=fit0.logLik,"coef"=NA,"predict"=NA))
		#    define the model as a string
			if(length(wh)==0) {
				model = base_model
			} else {
				model = paste(base_model, "+" ,paste0("a[,",wh,"]",collapse=" + "))
			}
		#    fit the model
			fit = glm(as.formula(model), family=glm_family, singular.ok=TRUE)
			# Get a prediction of the model
		        if(length(wh)<2) {
			     prediction = rep(1, nrow(a))

			} else {
			     prediction = predict(fit, a[,wh], type="response")			   }
		#    compute the log-likelihood
			fit.logLik = logLik(fit)
		#    obtain the coefficients
			fit.coef = coef(summary(fit))
		#    return value is the difference in log-likelihood to the null
			ret.loglik = fit.logLik-fit0.logLik
		#    has attribute equal to difference in degrees of freedom
			attributes(ret.loglik)$df = attributes(fit.logLik)$df - attributes(fit0.logLik)$df
		#    if any coefficient is NA indicates a singular fit: return -Inf
			if(any(is.na(coef(fit)))) ret.loglik = ret.loglik + log(0)
		#    return
			return(list("loglik"=ret.loglik,"coef"=fit.coef,"predict"=prediction))
		}
	}; calc.loglik = gen.calc.loglik(flat_likelihood=flat_likelihood)
#            calculate a negative log p-value for each model vs the null
	calc.neglogp = function(loglik) as.numeric(-pchisq(2*loglik,attributes(loglik)$df,low=FALSE,log=TRUE))
#            calculate a log Bayes factor for each model vs the null
	calc.logBF = function(loglik) as.numeric(0.5*attributes(loglik)$df*log(xi.prior) + (1-xi.prior)*loglik)
#            log proposal probability for moves that add a regressor
#            needed for the MCMC Hastings ratio
	logpr.addmove = function(n,mu=MU) {
		pr.swap = 0.1
		pr.add = (1.0-pr.swap)*0.5
		if(n==0) pr.add = 1.0
		if(n==length(mvec.gd)) pr.add = 0.0
		log(pr.add) - lchoose(length(mvec.gd)-n,1)
	}
#            log proposal probability for moves that remove a regressor
#            needed for the MCMC Hastings ratio
	logpr.remmove = function(n,mu=MU) {
		# Bug fixed 15 Feb 2022: need 0.95 to account for swap moves
		pr.rem = 0.475
		if(n==0) pr.rem = 0
		if(n==length(mvec.gd)) pr.rem = 1
		log(pr.rem) - lchoose(n,1)
	}
#            convert the binary encoding of regressors to include to a string
	mvec2string = function(mvec) paste0(mvec,collapse="")
	
# Initialize: iteration number
	iter=0
#             storage to record the log-likelihood
	loglik.rec = list()
#             storage for the thinned model parameters
beta = list()
#             storage to record how often the MCMC stays in each state
	mcmc.occupancy = list()
#				storage to record number of variables included in each model
	loglik.rec.n = list()
#				storage to record log prior probability
	logPrP = list()
#				storage to record log Bayes factor vs the grand null
	logBF = list()
#				storage to record degrees of freedom vs the grand null
	df.rec = list()
#             storage to record the MCMC state for each iteration
	state = rep(NA,niter+1)
#             storage to record each proposed MCMC move
	mcmc.move = rep(NA,niter)
#             storage to record the acceptance probability for each proposed MCMC move
	mcmc.logalpha = rep(NA,niter)
#             storage to record the MCMC state for each proposed move
	state.proposed = rep(NA,niter)
#             variable encoding inclusion status of each column of a
#             initialize with 10 randomly drawn columns, excluding
#             columns not allowed to be regressors
	start_group_file = paste0(workdir,"/",outname,"/", outname, "_groups.csv")
	mvec = as.vector(as.matrix(read.csv(start_group_file, nrows=1, skip=as.integer(chainid)-1)))
	#mvec = mvec[1:length(mvec.gd)]
	mvec[is.na(match(1:length(mvec),mvec.gd))] = 0
#             convert mvec to a string
	string.mvec = mvec2string(mvec)
#             record the initial log-likelihood and coefficients
	fit.loglik.coef = calc.loglik(mvec)
	loglik.rec[[string.mvec]] = fit.loglik.coef$loglik
	beta[[string.mvec]] = fit.loglik.coef$coef
	mvec.n = sum(mvec)
	loglik.rec.n[[string.mvec]] = mvec.n
	logPrP[[string.mvec]] = logpr.prior(mvec.n)
	logBF[[string.mvec]] = calc.logBF(loglik.rec[[string.mvec]])
	df.rec[[string.mvec]] = attributes(loglik.rec[[string.mvec]])$df
#             record the initial MCMC state
	state[1] = match(string.mvec,names(loglik.rec))
#             record the initial MCMC occupancy
	mcmc.occupancy[[string.mvec]] = 1

	# Run the MCMC
	for(iter in 1:niter) {
	# while(TRUE)  iter=iter+1  # replace the above line for indeterminate chain length
		# Create a variable to store the proposed MCMC state
		mvec.new = mvec
		# Create a variable storing the previous MCMC state
		mvec.old = mvec
		# Count the number of regressors in the curren
		mvec.n = sum(mvec)

		# Simulate the MCMC proposal:
		#  do.swap  do.add  Proposal
		#     TRUE    TRUE  swap      (i.e. ignore do.add if do.swap==TRUE)
		#     TRUE   FALSE  swap
		#    FALSE    TRUE  add
		#    FALSE   FALSE  remove
		# The move probabilities need to match logpr.addmove etc
		do.swap = runif(1)<=0.1
		do.add = runif(1)<=0.5
		# Disallow some moves in some cases
		if(mvec.n==0) {
			do.swap=FALSE
			do.add=TRUE
		}
		if(mvec.n==length(mvec.gd)) {
			do.swap=FALSE
			do.add=FALSE
		}
		# Implement the move
		if(do.swap) {
			# This move replaces one regressor with another, correlated variable
			mcmc.move[iter] = "swap"
			# List of possible variables to remove
			choosefrom = which(mvec==1)
			# Choose the variable to remove
			ch1 = choosefrom[sample.int(length(choosefrom), 1)]
			# Sanity check
			stopifnot(mvec[ch1]==1)
			# List of possible variables to add
			choosefrom = setdiff(mvec.gd, which(mvec==1))
            # Transform a.rsq to modify the relative swap probabilities
            ftt = function(rsq) rsq^4
			# Choose to add a variable correlated with the removed variable
			ch2 = choosefrom[sample.int(length(choosefrom), 1, prob=pmax(0.001, ftt(a.rsq[ch1,choosefrom])))]
			# Record the proposed move
			mvec.new[ch1] = 0
			mvec.new[ch2] = 1
			# To calculate the Hastings ratio, need to consider the possible
			# backwards move probability
			reverse.choosefrom = setdiff(mvec.gd, which(mvec.new==1))
			# Calculate the log-Hastings ratio
			swap.Hastings = -log(sum(pmax(0.001, ftt(a.rsq[ch2, reverse.choosefrom])))) + log(sum(pmax(0.001, ftt(a.rsq[ch1, choosefrom]))))
		} else if(do.add) {
			# This move adds a regressor to the model
			mcmc.move[iter] = "add"
			# List of possible variables to choose from
			choosefrom = setdiff(mvec.gd,which(mvec==1))
			# Choose the variable to add
			ch = choosefrom[sample.int(length(choosefrom),1)]
			# Sanity check
			stopifnot(mvec[ch]==0)
			# Record the proposed move
			mvec.new[ch] = 1
		} else {
			# This move removes a regressor from the model
			mcmc.move[iter] = "remove"
			# List of possible variables to remove
			choosefrom = which(mvec==1)
			# Choose the variable to remove
			ch = choosefrom[sample.int(length(choosefrom),1)]
			# Sanity check
			stopifnot(mvec[ch]==1)
			# Record the proposed move
			mvec.new[ch] = 0
		}
		# Convert mvec.new to a string
		string.mvec.new = mvec2string(mvec.new)
		# Count the number of regressors in the current MCMC state
		mvec.new.n = sum(mvec.new)
		# Output the proposed move to screen if desired
		if(report.move) cat(iter,"Proposing to",ifelse(do.swap,paste("swap",colnames(a)[ch1],colnames(a)[ch2]),paste(ifelse(mvec.new[ch]==1,"add","remove"),colnames(a))[ch]),"...")
		# Calculate and record the log-likelihood if not previously computed
		if(is.null(loglik.rec[[string.mvec.new]])) {
			fit.loglik.coef = calc.loglik(mvec.new)
			loglik.rec[[string.mvec.new]] = fit.loglik.coef$loglik
			beta[[string.mvec.new]] = fit.loglik.coef$coef
			prediction = fit.loglik.coef$predict
			mcmc.occupancy[[string.mvec.new]] = 0
			loglik.rec.n[[string.mvec.new]] = mvec.new.n
			logPrP[[string.mvec.new]] = logpr.prior(mvec.new.n)
			logBF[[string.mvec.new]] = calc.logBF(loglik.rec[[string.mvec.new]])
			df.rec[[string.mvec.new]] = attributes(loglik.rec[[string.mvec.new]])$df
		}
		# Calculate the log-Hastings ratio
		if(do.swap) {
			hastings = swap.Hastings
		} else if(do.add) {
			hastings = logpr.remmove(mvec.new.n)-logpr.addmove(mvec.n)
		} else {
			hastings = logpr.addmove(mvec.new.n)-logpr.remmove(mvec.n)
		}
		# Calculate the Metropolis-Hastings ratio
		mh.ratio = logPrP[[string.mvec.new]] - logPrP[[string.mvec]] + logBF[[string.mvec.new]] - logBF[[string.mvec]] + hastings
		# Accept or reject
		if(is.na(mh.ratio)) {
			print("doublethink-functions: lr.mcmc(): ERROR MH Ratio is NA. Exiting")
			print("logPrP stin.mvec.new")
			print(logPrP[[string.mvec.new]])
			print("logPrP stin.mvec")
			print(logPrP[[string.mvec]])
			print("logBF string mvec new")
			print(logBF[[string.mvec.new]])
			print("logBF string mvec")
			print(logBF[[string.mvec]])
			print("hastings")
			print(hastings)
			print("mh ratio")
			print(mh.ratio)
			print("string mvec new")
			print(string.mvec.new)
			print("string mvec")
			print(string.mvec)
		}
		if(log(runif(1))<=mh.ratio) { # accept
			# Update the current MCMC state
			mvec = mvec.new
			string.mvec = string.mvec.new
			# Update the MCMC occupancy record
			mcmc.occupancy[[string.mvec.new]] = mcmc.occupancy[[string.mvec.new]] + 1
			# Report to screen if desired
			if(report.move) cat(" Accepted\n")
		} else { # reject
			# Update the MCMC occupancy record
			mcmc.occupancy[[string.mvec]] = mcmc.occupancy[[string.mvec]] + 1
			# Report to screen if desired
			if(report.move) cat(" Rejected\n")
		}
		# Not used unless MCMC chain length is indeterminate
		if(length(state)==iter) {
			state = c(state,rep(NA,niter))
			state.proposed = c(state.proposed,rep(NA,niter))
			mcmc.move = c(mcmc.move,rep(NA,niter))
			mcmc.logalpha = c(mcmc.logalpha,rep(NA,niter))
		}
		# Update the MCMC state, the proposed MCMC state and the acceptance probability
		state[iter+1] = match(string.mvec,names(loglik.rec))
		state.proposed[iter] = match(string.mvec.new,names(loglik.rec))
		mcmc.logalpha[iter] = mh.ratio
		# Report the new MCMC state to screen, if desired
		if((iter%%1000==0)) {
			# Identify regressors whose inclusion status has changed
			tp = colnames(a)[mvec==1 | mvec.old==1]
			# Storage to colour the column names
			col = rep("",length(tp))
			# Colour additions green
			col[is.na(match(tp,colnames(a[mvec.old==1])))] = "\\e[32m"
			# Colour removals light red
			col[is.na(match(tp,colnames(a[mvec==1])))] = "\\e[91m"
			# System command to print to screen (needed for colours)
			out = paste("echo -e '",iter,":\n",paste(col,tp,"\\e[0m",collapse="\n",sep=""),"\n'",sep="")
			# Execute the system command
			system(out)
			prediction_table[nrow(prediction_table) +1,] = prediction
			# After 100 iterations, only report to screen periodically (if at all)
			if(iter>100) report.state=FALSE
		}
	}
	prediction_table = cbind(rownames(a), colMeans(prediction_table))
	write.table(format(prediction_table,digits=5), file = prediction_file, sep = ",", quote = FALSE)

	ret = list("colnames.a"=colnames(a), "MU"=MU, "samp.size"=samp.size, "f.prior"=f.prior, "xi.prior"=xi.prior, "niter"=niter, "mvec.gd"=mvec.gd, "nmax"=length(mvec.gd), "logpr.prior"=logpr.prior, "calc.loglik"=calc.loglik, "calc.logBF"=calc.logBF, "logpr.addmove"=logpr.addmove, "logpr.remmove"=logpr.remmove, "loglik.rec"=loglik.rec, "beta"=beta, "mcmc.occupancy"=mcmc.occupancy, "loglik.rec.n"=loglik.rec.n, "logPrP"=logPrP, "logBF"=logBF, "df.rec"=df.rec, "state"=state, "mcmc.move"=mcmc.move, "state.proposed"=state.proposed, "mcmc.logalpha"=mcmc.logalpha, "nlev"=nlev)
	return(ret)
}


############################################
#              HMP functions               #
############################################


# Calculate posterior odds and harmonic mean p-values per chain
# (importance sampling method using removal moves, assumes high model prob)
hmp.stats.chain = function(rds, mvec.gd.j, burnin=0) {
	j = rds$mvec.gd[mvec.gd.j]
	# Storage indices of models excluding variable j
	j.state.neg = which(substr(names(rds$mcmc.occupancy),j,j)=="0")
	# Models counterparts with variable j included
	j.state.pos.names = sapply(names(rds$mcmc.occupancy)[j.state.neg],function(s) {ret=s; substr(ret,j,j)="1"; return(ret)})
	# Storage indices of models corresponding to j.state.neg
	# but including variable j. NA if not explored
	j.state.pos = match(j.state.pos.names,names(rds$mcmc.occupancy))

	# Iterations in which j.state.neg was proposed from j.state.pos
	# removing burn-in
	gd = (1:rds$niter)>burnin & rds$state[-rds$niter]==j.state.pos[match(rds$state.proposed,j.state.neg)]
	gd = which(!is.na(gd) & gd)
	stopifnot(all(rds$mcmc.move[gd]=="remove"))
	if(length(gd)==0) {
		return(c("log.hmp"=NA,"log.PoO"=NA,"log.PoO.approx"=NA, "M"=length(gd),"log.inv.hmp.num"=NA,"log.inv.hmp.den"=NA, "log.PoO.num"=NA, "log.PoO.den"=NA,"log.PoO.approx.num"=NA,  "log.PoO.approx.den"=NA))
	}
	
	# The proposal weights (unnormalized)
	# Proportional to posterior prob. of current state (j included)
	# times the proposal prob. of removing one variable (j)
	log.proposal = unlist(rds$logPrP[rds$state[gd]]) + unlist(rds$logBF[rds$state[gd]]) + log(ifelse(unlist(rds$loglik.rec.n[rds$state[gd]])==rds$nmax,1,0.475)) - log(unlist(rds$loglik.rec.n[rds$state[gd]]))
	
	# Precompute quantities for the function h
	# Prior odds (of including vs excluding variable j when !is.na(iter.rem))
	logpro = unlist(rds$logPrP[rds$state[gd]]) - unlist(rds$logPrP[rds$state.proposed[gd]])
	# Difference in degrees of freedom (should always be the same when !is.na(iter.rem))
	delta.df = unlist(rds$df.rec[rds$state[gd]]) - unlist(rds$df.rec[rds$state.proposed[gd]])
	# Compute the p-values
	neglogp = sapply(gd,function(iter) calc.neglogp.pair(loglik.0=rds$loglik.rec[[rds$state.proposed[iter]]], loglik.A=rds$loglik.rec[[rds$state[iter]]]))

	# The function h
	# log.h.PoO for the exact PoO:
	log.h.PoO = logpro + unlist(rds$logBF[rds$state[gd]]) - unlist(rds$logBF[rds$state.proposed[gd]])
	# log.h.PoO.approx for the approximate PoO:
	# Define BF ~ BF0*p0/p where BF(p0)==Pr(M-)/Pr(M+),
	# i.e. exact when posterior odds = 1 (Method 4b of setup.R)
	# Specify the point around which to approximate the Bayes factor
	logBF0 = -logpro
	neglogp0 = -pchisq(2/(1-rds$xi.prior)*(logBF0-delta.df/2*log(rds$xi.prior)), delta.df,low=FALSE,log=TRUE)
	# h = PrO * approxBF = PrO * BF0 * p0 / p = PrO * 1/PrO * p0 / p = p0 / p
	log.h.PoO.approx = neglogp - neglogp0
	# log.h.hmp for the (inverse) HMP:
	log.h.inv.hmp = neglogp
	
	# The target weights (unnormalized) - note that log.weight + log.h is (approx) identical for all three
	# For PoO and approx PoO: Proportional to posterior prob. of proposed model (j excluded when !is.na(iter.rem))
	log.weight.PoO = unlist(rds$logPrP[rds$state.proposed[gd]]) + unlist(rds$logBF[rds$state.proposed[gd]])
	# For HMP: includes everything except 1/p itself
	log.weight.hmp = unlist(rds$logPrP[rds$state.proposed[gd]]) + unlist(rds$logBF[rds$state.proposed[gd]]) - neglogp0

	# Sampling estimate of the posterior odds of the model including vs excluding variable j
	log.PoO.num = logsum(log.weight.PoO - log.proposal + log.h.PoO, na.rm=TRUE)
	log.PoO.den = logsum(log.weight.PoO - log.proposal, na.rm=TRUE)
	log.PoO = log.PoO.num - log.PoO.den
	log.PoO.approx.num = logsum(log.weight.PoO - log.proposal + log.h.PoO.approx, na.rm=TRUE)
	log.PoO.approx.den = logsum(log.weight.PoO - log.proposal, na.rm=TRUE)
	log.PoO.approx = log.PoO.approx.num - log.PoO.approx.den
	log.inv.HMP.num = logsum(log.weight.hmp - log.proposal + log.h.inv.hmp, na.rm=TRUE)
	log.inv.HMP.den = logsum(log.weight.hmp - log.proposal, na.rm=TRUE)
	log.inv.HMP = log.inv.HMP.num - log.inv.HMP.den

	return(c("log.hmp"=-log.inv.HMP,"log.PoO"=log.PoO,"log.PoO.approx"=log.PoO.approx,"M"=length(gd),"log.inv.hmp.num"=log.inv.HMP.num,"log.inv.hmp.den"=log.inv.HMP.den,"log.PoO.num"=log.PoO.num, "log.PoO.den"=log.PoO.den,"log.PoO.approx.num"=log.PoO.approx.num, "log.PoO.approx.den"=log.PoO.approx.den))
}
hmp.merge = function(list.hmp.stats, logL=NULL, minM=30, quiet=FALSE) {
	if(length(list.hmp.stats)==0) {
		return(c(
			"log.PoO.pe"=NA,
			"log.PoO.stderr"=NA,
			"log.PoO.approx.pe"=NA,
			"log.PoO.approx.stderr"=NA,
			"log.raw.hmp.pe"=NA,
			"log.raw.hmp.stderr"=NA,
			"log.raw.hmp.bayes.thresh"=NA,
			"log.p.hmp"=NA,
			"mn.M"=NA,
			"sd.M"=NA,
			"nchain"=0
		))
	}
	hmp.j = matrix(list.hmp.stats[[1]],ncol=1)
	rownames(hmp.j) = names(list.hmp.stats[[1]])
	if(length(list.hmp.stats)>1) for(i in 2:length(list.hmp.stats)) hmp.j = cbind(hmp.j, list.hmp.stats[[i]])
	
	gd = hmp.j["M",]>minM
	if(any(!gd) & !quiet) print(paste0("WARNING in hmp.merge(): ",sum(!gd)," out of ",length(gd)," chains had fewer than ",minM," iterations"))
	log.hmp.stat.pe =
	-(logsum(hmp.j["log.inv.hmp.num",gd])-logsum(hmp.j["log.inv.hmp.den",gd]))
	if (sqrt(sum(gd)) >1){
		log.PoO.stderr = sd(hmp.j["log.PoO",gd])/sqrt(sum(gd))
		log.hmp.stat.stderr = sd(hmp.j["log.hmp",gd])/sqrt(sum(gd))
		log.hmp.stat.1.stderr = sd(hmp.j["log.hmp.1",gd])/sqrt(sum(gd))
	} else {
		log.hmp.stat.stderr = 0
		log.PoO.stderr = 0
		log.hmp.stat.1.stderr = 0
	}
	log.PoO.pe = logsum(hmp.j["log.PoO.num",gd])-logsum(hmp.j["log.PoO.den",gd])
	log.PoO.stderr = sd(hmp.j["log.PoO",gd])/sqrt(sum(gd))
	log.PoO.approx.pe = logsum(hmp.j["log.PoO.approx.num",gd])-logsum(hmp.j["log.PoO.approx.den",gd])
	log.PoO.approx.stderr = sd(hmp.j["log.PoO.approx",gd])/sqrt(sum(gd))

	if(is.null(L)) {
		log.p.hmp = NA
	} else {
		log.p.hmp = pharm_workaround(exp(log.hmp.stat.pe),logL=logL,log=TRUE)
	}
	# Implied HMP threshold - not simply the prior odds. It's the p-value
	# needed to make the approx PoO above 1. It's a mean of the p0 values.
	# ApproxBF_i = BF0_i * p0_i / p_i
	# ApproxPoO_i = PrO_i * BF0_i * p0_i / p_i
	#             = p0(BF=1/PrO_i) / p_i
	# So the implied threshold is the difference in denominator of 1/HMP (includes p0) vs ApproxPoO (excludes it):
	log.hmp.bayes.thresh = logsum(hmp.j["log.inv.hmp.den",gd]) - logsum(hmp.j["log.PoO.approx.den",gd])
	
	mn.M = mean(hmp.j["M",gd])
	sd.M = sd(hmp.j["M",gd])
	nchain = sum(gd)
	
	return(c(
		"log.PoO.pe"=log.PoO.pe,
		"log.PoO.stderr"=log.PoO.stderr,
		"log.PoO.approx.pe"=log.PoO.approx.pe,
		"log.PoO.approx.stderr"=log.PoO.approx.stderr,
		"log.raw.hmp.pe"=log.hmp.stat.pe,
		"log.raw.hmp.stderr"=log.hmp.stat.stderr,
		"log.raw.hmp.bayes.thresh"=log.hmp.bayes.thresh,
		"log.p.hmp"=log.p.hmp,
		"mn.M"=mn.M,
		"sd.M"=sd.M,
		"nchain"=nchain
	))
}

# Calculate posterior odds and harmonic mean p-values per chain
# (importance sampling method using removal moves, assumes high model prob - second approach)
hmp.stats.chain2 = function(rds, mvec.gd.j, burnin=0) {
	# Identify the variable of interest
	j = rds$mvec.gd[mvec.gd.j]
	# Odds (based on calc.PoP.chain)
	# NB: the first state is always ignored
	wh = which(substr(names(rds$mcmc.occupancy),j,j)=="1")
	occupancy.pos = sum(!is.na(match(rds$state[-(1:(burnin+1))],wh)))
	occupancy.neg = rds$niter-burnin-occupancy.pos
	# Sampling estimate of the posterior odds of the model including vs excluding variable j
	log.PoO.num = log(occupancy.pos)
	log.PoO.den = log(occupancy.neg)
	log.PoO = log.PoO.num - log.PoO.den

	# Denominator
	# Storage indices of models excluding variable j
	j.state.neg = which(substr(names(rds$mcmc.occupancy),j,j)=="0")
	# Models counterparts with variable j included
	j.state.pos.names = sapply(names(rds$mcmc.occupancy)[j.state.neg],function(s) {ret=s; substr(ret,j,j)="1"; return(ret)})
	# Storage indices of models corresponding to j.state.neg
	# but including variable j. NA if not explored
	j.state.pos = match(j.state.pos.names,names(rds$mcmc.occupancy))

	# Iterations in which j.state.neg was proposed from j.state.pos
	# removing burn-in
	gd = (1:rds$niter)>burnin & rds$state[-rds$niter]==j.state.pos[match(rds$state.proposed,j.state.neg)]
	gd = which(!is.na(gd) & gd)
	stopifnot(all(rds$mcmc.move[gd]=="remove"))
	if(length(gd)==0) {
		return( c("log.hmp"=NA,"log.hmp.1"=NA,"log.PoO"=log.PoO,"log.PoO.approx"=NA,"log.PoO.num"=log.PoO.num,"log.PoO.den"=log.PoO.den,"log.mean.logp0"=NA,"log.mean.logp0.1"=NA,"log.nonlinearity"=NA,"n.num"=NA,"n.num.1"=NA,"n.den"=length(gd),"n.PoO"=rds$niter-burnin,"log.mean.logp0.num"=NA,"log.mean.logp0.den"=NA,"log.nonlinearity.num"=NA,"log.nonlinearity.den"=NA) )
	}
	
	# The proposal weights (unnormalized)
	# Proportional to posterior prob. of current state (j included)
	# times the proposal prob. of removing one variable (j)
	log.proposal = unlist(rds$logPrP[rds$state[gd]]) + unlist(rds$logBF[rds$state[gd]]) + log(ifelse(unlist(rds$loglik.rec.n[rds$state[gd]])==rds$nmax,1,0.475)) - log(unlist(rds$loglik.rec.n[rds$state[gd]]))
	
	# Precompute quantities for the function h
	# Prior odds (of including vs excluding variable j when !is.na(iter.rem))
	logpro = unlist(rds$logPrP[rds$state[gd]]) - unlist(rds$logPrP[rds$state.proposed[gd]])
	# Difference in degrees of freedom (should always be the same when !is.na(iter.rem))
	delta.df = unlist(rds$df.rec[rds$state[gd]]) - unlist(rds$df.rec[rds$state.proposed[gd]])
	# Compute the p-values
	logp = sapply(gd,function(iter) -calc.neglogp.pair(loglik.0=rds$loglik.rec[[rds$state.proposed[iter]]], loglik.A=rds$loglik.rec[[rds$state[iter]]]))

	# Define BF ~ BF0*p0/p where BF(p0)==Pr(M-)/Pr(M+),
	# i.e. exact when posterior odds = 1 (Method 4b of setup.R)
	# Specify the point around which to approximate the Bayes factor
	logBF0 = -logpro
	logp0 = pchisq(2/(1-rds$xi.prior)*(logBF0-delta.df/2*log(rds$xi.prior)), delta.df,low=FALSE,log=TRUE)
	
	# Final term for the denominator
	#den.mean.nonlinearity = logmean( logp0 - logp - (log.proposal - unlist(rds$logPrP[rds$state.proposed[gd]]) - unlist(rds$logBF[rds$state.proposed[gd]]))   , na.rm=TRUE)
	den.nonlinearity.num = logsum(unlist(rds$logPrP[rds$state.proposed[gd]]) + unlist(rds$logBF[rds$state.proposed[gd]]) + logp0 - logp - log.proposal)
	den.nonlinearity.den = logsum(unlist(rds$logPrP[rds$state[gd]]) + unlist(rds$logBF[rds$state[gd]]) - log.proposal)
	den.mean.nonlinearity = den.nonlinearity.num - den.nonlinearity.den

	# Numerator
	# Assert the DF is always the same
	stopifnot(all(delta.df)==delta.df[1])
	# Iterations visiting models excluding variable j
	# Exclude the first iteration for consistency with earlier
	num.gd = which((1:rds$niter)>burnin & !is.na(match(rds$state[-1],j.state.neg)))
	num.n = length(num.gd)
	# Prior odds (of including vs excluding variable j)
	logpr.prior = function(n,mu=rds$MU,nmax=rds$nmax) dbinom(n,nmax,1/mu/nmax,log=TRUE)-lchoose(nmax,n) # binomial distribution
	if(num.n>0) {
		# Number of variables present among models excluding variable j
		num.nvar = unlist(rds$loglik.rec.n[rds$state[num.gd]])
		num.logpro = logpr.prior(num.nvar+delta.df[1])-logpr.prior(num.nvar)
		num.logBF0 = -num.logpro
		num.logp0 = pchisq(2/(1-rds$xi.prior)*(num.logBF0-delta.df[1]/2*log(rds$xi.prior)), delta.df[1],low=FALSE,log=TRUE)
		num.mean.logp0 = logmean(num.logp0)
	} else {
		num.mean.logp0 = NA
	}

	# Numerator version 2
	# Storage indices of models including variable j
	add.j.state.pos = which(substr(names(rds$mcmc.occupancy),j,j)=="1")
	# Models counterparts with variable j excluded
	add.j.state.neg.names = sapply(names(rds$mcmc.occupancy)[add.j.state.pos],function(s) {ret=s; substr(ret,j,j)="0"; return(ret)})
	# Storage indices of models corresponding to j.state.pos
	# but excluding variable j. NA if not explored
	add.j.state.neg = match(add.j.state.neg.names,names(rds$mcmc.occupancy))

	# Iterations in which j.state.pos was proposed from j.state.neg
	# removing burn-in
	add.gd = (1:rds$niter)>burnin & rds$state[-rds$niter]==add.j.state.neg[match(rds$state.proposed,add.j.state.pos)]
	add.gd = which(!is.na(add.gd) & add.gd)
	stopifnot(all(rds$mcmc.move[add.gd]=="add"))
	stopifnot(length(intersect(gd,add.gd))==0)
	# Unnormalized proposal probability: remove then add moves
	num2.logproposal = c(log.proposal,
	#	Proportional to posterior prob. of current state (j excluded)
	#	times the proposal prob. of adding one variable (j)
		unlist(rds$logPrP[rds$state[add.gd]]) + unlist(rds$logBF[rds$state[add.gd]]) + log(ifelse(unlist(rds$loglik.rec.n[rds$state[add.gd]])==0,1,0.475)) - log(rds$nmax-unlist(rds$loglik.rec.n[rds$state[add.gd]]))
	)
	# Number of variables present among models excluding variable j
	num2.nvar = c(unlist(rds$loglik.rec.n[rds$state.proposed[gd]]), unlist(rds$loglik.rec.n[rds$state[add.gd]]))
	# Log prior odds of adding variable j to those base models
	num2.logpro = logpr.prior(num2.nvar+delta.df[1])-logpr.prior(num2.nvar)
	num2.logBF0 = -num2.logpro
	# The p-value thresholds for adding variable j to those base models
	num2.logp0 = pchisq(2/(1-rds$xi.prior)*(num2.logBF0-delta.df[1]/2*log(rds$xi.prior)), delta.df[1],low=FALSE,log=TRUE)
	# The posterior odds of those base models
	num2.PoO.base = c(unlist(rds$logPrP[rds$state.proposed[gd]]) + unlist(rds$logBF[rds$state.proposed[gd]]),
		unlist(rds$logPrP[rds$state[add.gd]]) + unlist(rds$logBF[rds$state[add.gd]]))
	# The numerator excluding the posterior odds of excluding variable j
	num2.mean.logp0.num = logsum(num2.logp0 + num2.PoO.base - num2.logproposal)
	num2.mean.logp0.den = logsum(num2.PoO.base - num2.logproposal)
	num2.mean.logp0 = num2.mean.logp0.num - num2.mean.logp0.den
	num2.n = length(gd)+length(add.gd)
	
	# HMP
	log.HMP = -log.PoO + num.mean.logp0 - den.mean.nonlinearity
	log.HMP.num2 = -log.PoO + num2.mean.logp0 - den.mean.nonlinearity
	# Approximate PoO: only for thresholding because its exact value is biased away from or towards 0, depending on degrees of freedom
	log.PoO.approx = num.mean.logp0 - log.HMP

	return(c("log.hmp"=log.HMP.num2,"log.hmp.1"=log.HMP,"log.PoO"=log.PoO,"log.PoO.approx"=log.PoO.approx,"log.PoO.num"=log.PoO.num,"log.PoO.den"=log.PoO.den,"log.mean.logp0"=num2.mean.logp0,"log.mean.logp0.1"=num.mean.logp0,"log.nonlinearity"=den.mean.nonlinearity,"n.num"=num2.n,"n.num.1"=num.n,"n.den"=length(gd),"n.PoO"=rds$niter-burnin,"log.mean.logp0.num"=num2.mean.logp0.num,"log.mean.logp0.den"=num2.mean.logp0.den,	"log.nonlinearity.num"=den.nonlinearity.num,"log.nonlinearity.den"=den.nonlinearity.den))
}

# Merge HMPs not used atm 
hmp.merge2 = function(list.hmp.stats, logL=NULL, minM=30, quiet=FALSE) {
	if(length(list.hmp.stats)==0) {
		return(c(
			"log.PoO.pe"=NA,
			"log.PoO.stderr"=NA,
			"log.PoO.med"=NA,
			"log.PoO.approx.pe"=NA,
			"log.PoO.approx.stderr"=NA,
			"log.PoO.approx.med"=NA,
			"log.raw.hmp.pe"=NA,
			"log.raw.hmp.stderr"=NA,
			"log.raw.hmp.med"=NA,
			"log.raw.hmp.1.pe"=NA,
			"log.raw.hmp.1.stderr"=NA,
			"log.raw.hmp.1.med"=NA,
			"log.raw.hmp.bayes.thresh"=NA,
			"log.p.hmp"=NA,
			"log.p.hmp.1"=NA,
			"mn.n.PoO"=NA,
			"sd.n.PoO"=NA,
			"mn.n.num"=NA,
			"sd.n.num"=NA,
			"mn.n.num.1"=NA,
			"sd.n.num.1"=NA,
			"mn.n.den"=NA,
			"sd.n.den"=NA,
			"nchain"=0
		))
	}
	hmp.j = matrix(list.hmp.stats[[1]],ncol=1)
	rownames(hmp.j) = names(list.hmp.stats[[1]])
	if(length(list.hmp.stats)>1) for(i in 2:length(list.hmp.stats)) hmp.j = cbind(hmp.j, list.hmp.stats[[i]])
	
	gd = hmp.j["n.den",]>minM
	if(any(!gd) & !quiet) print(paste0("WARNING in hmp.merge(): ",sum(!gd)," out of ",length(gd)," chains had fewer than ",minM," iterations"))
	log.hmp.stat.pe = logsum(hmp.j["log.PoO.den",gd]) - logsum(hmp.j["log.PoO.num",gd]) + logsum(hmp.j["log.mean.logp0.num",gd]) - logsum(hmp.j["log.mean.logp0.den",gd]) - logsum(hmp.j["log.nonlinearity.num",gd]) + logsum(hmp.j["log.nonlinearity.den",gd])
	log.hmp.stat.1.pe = logsum(hmp.j["log.PoO.den",gd]) - logsum(hmp.j["log.PoO.num",gd]) + logsum(hmp.j["log.mean.logp0.1",gd]+log(hmp.j["n.num.1",gd])) - logsum(log(hmp.j["n.num.1",gd])) - logsum(hmp.j["log.nonlinearity.num",gd]) + logsum(hmp.j["log.nonlinearity.den",gd])
	log.hmp.stat.med = median(hmp.j["log.hmp",gd])
	log.hmp.stat.1.med = median(hmp.j["log.hmp.1",gd])
	log.PoO.pe = logsum(hmp.j["log.PoO.num",gd])-logsum(hmp.j["log.PoO.den",gd])
	if (sqrt(sum(gd)) >1){
		log.PoO.stderr = sd(hmp.j["log.PoO",gd])/sqrt(sum(gd))
		log.hmp.stat.stderr = sd(hmp.j["log.hmp",gd])/sqrt(sum(gd))
		log.hmp.stat.1.stderr = sd(hmp.j["log.hmp.1",gd])/sqrt(sum(gd))
	} else {
		log.hmp.stat.stderr = 0
		log.PoO.stderr = 0
		log.hmp.stat.1.stderr = 0
	}
	log.PoO.med = median(hmp.j["log.PoO",gd])
	log.PoO.approx.pe = -( logsum(hmp.j["log.PoO.den",gd])-logsum(hmp.j["log.PoO.num",gd]) - (logsum(hmp.j["log.nonlinearity",gd]+log(hmp.j["n.den",gd])) - logsum(log(hmp.j["n.den",gd]))) )
	log.PoO.approx.stderr = sd(hmp.j["log.PoO.approx",gd])/sqrt(sum(gd))
	log.PoO.approx.med = median(hmp.j["log.PoO.approx",gd])

	if(is.null(L)) {
		log.p.hmp = NA
		log.p.hmp.1 = NA
	} else {
		log.p.hmp = pharm_workaround(exp(log.hmp.stat.pe),logL=logL,log=TRUE)
		log.p.hmp.1 = pharm_workaround(exp(log.hmp.stat.1.pe),logL=logL,log=TRUE)
	}
	# Implied HMP threshold - not simply the prior odds. It's the p-value
	# needed to make the approx PoO above 1. It's a mean of the p0 values.
	# ApproxBF_i = BF0_i * p0_i / p_i
	# ApproxPoO_i = PrO_i * BF0_i * p0_i / p_i
	#             = p0(BF=1/PrO_i) / p_i
	# So the implied threshold is the difference in denominator of 1/HMP (includes p0) vs ApproxPoO (excludes it):
	log.hmp.bayes.thresh = logsum(hmp.j["log.mean.logp0.num",gd]) - logsum(hmp.j["log.mean.logp0.den",gd])
	log.hmp.1.bayes.thresh = logsum(hmp.j["log.mean.logp0.1",gd]+log(hmp.j["n.num.1",gd])) - logsum(log(hmp.j["n.num.1",gd]))
	
	mn.n.PoO = mean(hmp.j["n.PoO",gd])
	sd.n.PoO = sd(hmp.j["n.PoO",gd])
	mn.n.num = mean(hmp.j["n.num",gd])
	sd.n.num = sd(hmp.j["n.num",gd])
	mn.n.num.1 = mean(hmp.j["n.num.1",gd])
	sd.n.num.1 = sd(hmp.j["n.num.1",gd])
	mn.n.den = mean(hmp.j["n.den",gd])
	sd.n.den = sd(hmp.j["n.den",gd])
	
	nchain = sum(gd)
	
	return(c(
		"log.PoO.pe"=log.PoO.pe,
		"log.PoO.stderr"=log.PoO.stderr,
		"log.PoO.med"=log.PoO.med,
		"log.PoO.approx.pe"=log.PoO.approx.pe,
		"log.PoO.approx.stderr"=log.PoO.approx.stderr,
		"log.PoO.approx.med"=log.PoO.approx.med,
		"log.raw.hmp.pe"=log.hmp.stat.pe,
		"log.raw.hmp.stderr"=log.hmp.stat.stderr,
		"log.raw.hmp.med"=log.hmp.stat.med,
		"log.raw.hmp.1.pe"=log.hmp.stat.1.pe,
		"log.raw.hmp.1.stderr"=log.hmp.stat.1.stderr,
		"log.raw.hmp.1.med"=log.hmp.stat.1.med,
		"log.raw.hmp.bayes.thresh"=log.hmp.bayes.thresh,
		"log.p.hmp"=log.p.hmp,
		"log.p.hmp.1"=log.p.hmp.1,
		"mn.n.PoO"=mn.n.PoO,
		"sd.n.PoO"=sd.n.PoO,
		"mn.n.num"=mn.n.num,
		"sd.n.num"=sd.n.num,
		"mn.n.num.1"=mn.n.num.1,
		"sd.n.num.1"=sd.n.num.1,
		"mn.n.den"=mn.n.den,
		"sd.n.den"=sd.n.den,
		"nchain"=nchain
	))
}

############################################
#        MCMC evaluation functions         #
############################################

# Calculate posterior probabilities (direct method) per chain
calc.PoP.chain = function(rds, burnin=0) {
	regressor.occupancy = rep(0,length(rds$mvec.gd))
	for(j in 1:length(rds$mvec.gd)) {
		wh = which(substr(names(rds$mcmc.occupancy),rds$mvec.gd[j],rds$mvec.gd[j])=="1")
		# NB: the first state is always ignored
		regressor.occupancy[j] = sum(!is.na(match(rds$state[-(1:(burnin+1))],wh)))
	}
	total.occupancy = rds$niter-burnin
	ret = regressor.occupancy/total.occupancy
	names(ret) = rds$colnames.a[rds$mvec.gd]
	return(ret)
}

# Posterior covariance between inclusion probabilities
calc.regressor.cooccupancy = function(rds, burnin=0) {
	regressor2rec = sapply(1:length(rds$mvec.gd),function(j) {
			(substr(names(rds$mcmc.occupancy),rds$mvec.gd[j],rds$mvec.gd[j])=="1")
		})
	# NB: the first state is always ignored
	regressor2rec = regressor2rec[rds$state[-(1:(burnin+1))],]
	regressor.cooccupancy = crossprod(regressor2rec)
	total.occupancy = rds$niter-burnin

	ret = regressor.cooccupancy/total.occupancy
	rownames(ret) = rds$colnames.a[rds$mvec.gd]
	colnames(ret) = rds$colnames.a[rds$mvec.gd]
	return(ret)
}

# Posterior covariance between inclusion probabilities
calc.regressor.coexclusion = function(rds, burnin=0) {
	regressor2rec = sapply(1:length(rds$mvec.gd),function(j) {
			(substr(names(rds$mcmc.occupancy),rds$mvec.gd[j],rds$mvec.gd[j])=="1")
		})
	# NB: the first state is always ignored
	regressor2rec = regressor2rec[rds$state[-(1:(burnin+1))],]
	regressor.coexclusion = crossprod(1-regressor2rec)
	total.occupancy = rds$niter-burnin

	ret = regressor.coexclusion/total.occupancy
	rownames(ret) = rds$colnames.a[rds$mvec.gd]
	colnames(ret) = rds$colnames.a[rds$mvec.gd]
	return(ret)
}

calc.PoP.merge = function(list.rds) {
	regressor.occupancy = matrix(0,length(list.rds),length(list.rds[[1]]$mvec.mean))
	colnames(regressor.occupancy) = names(list.rds[[1]]$mvec.mean)
	for(i in 1:length(list.rds)) regressor.occupancy[i,] = list.rds[[i]]$mvec.mean
	mn.regressor.occupancy = colMeans(regressor.occupancy)
	sd.regressor.occupancy = apply(regressor.occupancy,2,sd)
	ret = data.frame("PoP"=mn.regressor.occupancy,"SE"=sd.regressor.occupancy/sqrt(nrow(regressor.occupancy)))
	rownames(ret)=colnames(regressor.occupancy)
	return(ret)
}

calc.regressor.cooccupancy.merge = function(list.rds) {
	E.regressor.occupancy = list.rds[[1]]$mvec.mean
	E.regressor.cooccupancy = list.rds[[1]]$mvec.cooc
	#E2.regressor.cooccupancy = list.rds[[1]]$mvec.cooc^2
	if(length(list.rds)>1) for(i in 2:length(list.rds)) {
		E.regressor.occupancy = E.regressor.occupancy + list.rds[[i]]$mvec.mean
		E.regressor.cooccupancy = E.regressor.cooccupancy + list.rds[[i]]$mvec.cooc
		#E2.regressor.cooccupancy = E2.regressor.cooccupancy + list.rds[[i]]$mvec.cooc^2
	}
	E.regressor.occupancy = E.regressor.occupancy/length(list.rds)
	E.regressor.cooccupancy = E.regressor.cooccupancy/length(list.rds)
	#E2.regressor.cooccupancy = E2.regressor.cooccupancy/length(list.rds)
	#sd.regressor.cooccupancy = sqrt(E2.regressor.cooccupancy - E.regressor.cooccupancy^2)
	cov.regressor.cooccupancy = E.regressor.cooccupancy - outer(E.regressor.occupancy,E.regressor.occupancy)
	cor.regressor.cooccupancy = cov.regressor.cooccupancy/outer(sqrt(E.regressor.occupancy*(1-E.regressor.occupancy)),sqrt(E.regressor.occupancy*(1-E.regressor.occupancy)))
	rownames(cor.regressor.cooccupancy) = rownames(list.rds[[1]]$mvec.cooc)
	colnames(cor.regressor.cooccupancy) = colnames(list.rds[[1]]$mvec.cooc)
	return(cor.regressor.cooccupancy)
}

calc.regressor.coexclusion.merge = function(list.rds) {
	E.regressor.occupancy = list.rds[[1]]$mvec.mean
	E.regressor.coexclusion = list.rds[[1]]$mvec.coex
	if(length(list.rds)>1) for(i in 2:length(list.rds)) {
		E.regressor.occupancy = E.regressor.occupancy + list.rds[[i]]$mvec.mean
		E.regressor.coexclusion = E.regressor.coexclusion + list.rds[[i]]$mvec.coex
	}
	E.regressor.occupancy = E.regressor.occupancy/length(list.rds)
	E.regressor.coexclusion = E.regressor.coexclusion/length(list.rds)
	cov.regressor.coexclusion = E.regressor.coexclusion - outer(1-E.regressor.occupancy,1-E.regressor.occupancy)
	cor.regressor.coexclusion = cov.regressor.coexclusion/outer(sqrt(E.regressor.occupancy*(1-E.regressor.occupancy)),sqrt(E.regressor.occupancy*(1-E.regressor.occupancy)))
	rownames(cor.regressor.coexclusion) = rownames(list.rds[[1]]$mvec.coex)
	colnames(cor.regressor.coexclusion) = colnames(list.rds[[1]]$mvec.coex)
	return(cor.regressor.coexclusion)
}

simulate.beta.chain = function(rds, mvec.i, burnin=0, thinning.interval=1) {
	stopifnot(length(rds$state)>1+burnin)
	col.a = rds$mvec.gd[mvec.i]
	lab = paste0("a[, ",col.a,"]")
	# Expecting one coefficient for covariates and nlev-1 for factors
	n.coef = ifelse(rds$nlev[col.a]==0,1,rds$nlev[col.a]-1)
	zero.ret = rep(0,n.coef)
	q = 1-rds$xi.prior
	# NB: the first state is always ignored
	ret = sapply(rds$state[seq(2+burnin,length(rds$state), by=thinning.interval)], function(STATE) {
		COEF = rds$beta[[STATE]]
		wh = which(!is.na(match(substr(rownames(COEF),1,nchar(lab)),lab)))
		if(length(wh)==0) return(zero.ret)
		stopifnot(length(wh)==n.coef)
		rnorm(length(wh),q*COEF[wh,"Estimate"],sqrt(q)*COEF[wh,"Std. Error"])
	})
	return(matrix(ret,nrow=n.coef))
}
summarize.beta.chain = function(rds, mvec.i, burnin=0, thinning.interval=1) {
	beta.sim = simulate.beta.chain(rds, mvec.i, burnin, thinning.interval)
	is.zero = apply(beta.sim==0,2,all)
	return(list("zero.prob"=mean(is.zero),"non.zero.mean"=rowMeans(beta.sim[,!is.zero,drop=FALSE]),"non.zero.sd"=apply(beta.sim[,!is.zero,drop=FALSE],1,sd)))
}
simulate.beta.merge = function(list.rds, mvec.i, burnin=0, thinning.interval=1) {
	beta.sim = simulate.beta.chain(list.rds[[1]], mvec.i, burnin, thinning.interval)
	if(length(list.rds)>1) for(i in 2:length(list.rds)) {
		beta.sim = cbind(beta.sim,simulate.beta.chain(list.rds[[i]], mvec.i, burnin, thinning.interval))
	}
	return(beta.sim)
}
summarize.beta.merge = function(list.rds, mvec.i, burnin=0, thinning.interval=1) {
	beta.sim = simulate.beta.merge(list.rds, mvec.i, burnin, thinning.interval)
	is.zero = apply(beta.sim==0,2,all)
	return(list("zero.prob"=mean(is.zero),"non.zero.mean"=rowMeans(beta.sim[,!is.zero,drop=FALSE]),"non.zero.sd"=apply(beta.sim[,!is.zero,drop=FALSE],1,sd)))
}

# Merge regression coefficients
summarize.beta.sim.merge = function(list.summary, lab) {
	first = TRUE
	for(i in 1:length(list.summary)) {
		beta.sim.elem = list.summary[[i]]$beta.sim[[lab]]
		if(!is.null(beta.sim.elem)) {
			if(first) {
				beta.sim = beta.sim.elem
				first = FALSE
			}
			else {
				beta.sim = cbind(beta.sim,beta.sim.elem)
			}
		}
	}
	if(first) {
		return(list("zero.prob"=NA, "non.zero.mean"=NA,"non.zero.sd"=NA))
	} else {
		is.zero = apply(beta.sim==0,2,all)
		return(list("zero.prob"=mean(is.zero), "non.zero.mean"=rowMeans(beta.sim[,!is.zero,drop=FALSE]),"non.zero.sd"=apply(beta.sim[,!is.zero,drop=FALSE],1,sd)))
	}
}
summarize.mcmc = function(mcmc.result, burnin, thinning, version=1) {
	# Set the number of iterations to remove as burn-in (requires visual inspection of the MCMC trace)
	mcmc.summary = list()
	mcmc.summary$burnin = burnin
	mcmc.summary$thinning = thinning
	# Output to screen the posterior inclusion probabilities for the most included regressors
	mcmc.summary$mvec.mean = calc.PoP.chain(mcmc.result, burnin=mcmc.summary$burnin)

	mcmc.summary$mvec.cooc = calc.regressor.cooccupancy(mcmc.result, mcmc.summary$burnin)
	mcmc.summary$mvec.coex = calc.regressor.coexclusion(mcmc.result, mcmc.summary$burnin)

	# Summarize parameters for the top hits
	# Anything with posterior odds greater than or equal to (1-MU)/MU
	od = order(-mcmc.summary$mvec.mean)
	# or the top 10, whichever is longer
	top = sort(od)

	for(TOP in top) {
		cTOP = as.character(TOP)
		mcmc.summary$beta.sim[[cTOP]] = simulate.beta.chain(mcmc.result, TOP, burnin=mcmc.summary$burnin, thinning.interval=mcmc.summary$thinning)
		if(version==1) mcmc.summary$hmp.stats[[cTOP]] = hmp.stats.chain(mcmc.result, TOP)
		if(version==2) mcmc.summary$hmp.stats[[cTOP]] = hmp.stats.chain2(mcmc.result, TOP)
	}
	
	# Output full-model occupancy counts
	mcmc.summary$mcmc.occupancy = sort(table(mcmc.result$state[-(1:(burnin+1))]),decreasing=TRUE)
	names(mcmc.summary$mcmc.occupancy) = names(mcmc.result$loglik.rec)[as.numeric(names(mcmc.summary$mcmc.occupancy))]
	
	return(mcmc.summary)
}
calc.model.PoP.merge = function(list.summary) {
	all.models = names(list.summary[[1]]$mcmc.occupancy)
	if(length(list.summary)>1) for(i in 2:length(list.summary)) all.models = c(all.models,names(list.summary[[i]]$mcmc.occupancy))
	all.models = unique(all.models)
	E.occupancy = rep(0, length(all.models))
	E2.occupancy = rep(0, length(all.models))
	M.occupancy = rep(0, length(all.models))
	for(i in 1:length(list.summary)) {
		E.occupancy.i = as.vector(list.summary[[i]]$mcmc.occupancy[all.models])
		E.occupancy.i[is.na(E.occupancy.i)] = 0
		E.occupancy.i = E.occupancy.i/sum(E.occupancy.i)
		E.occupancy = E.occupancy + E.occupancy.i
		E2.occupancy = E2.occupancy + E.occupancy.i^2
		M.occupancy = M.occupancy + (E.occupancy.i>0)
	}
	E.occupancy = E.occupancy/length(list.summary)
	E2.occupancy = E2.occupancy/length(list.summary)
	SE.occupancy = sqrt( (E2.occupancy - E.occupancy^2) / length(list.summary) )
	names(E.occupancy) = all.models
	names(SE.occupancy) = all.models
	names(M.occupancy) = all.models
	return(data.frame("PoP"=E.occupancy,"SE"=SE.occupancy,"nchains"=M.occupancy))
}

##############################
# Exhaustive Model selection #
##############################
# pheno				vector of phenotypes (binary)
# a				data frame of regressors, # rows equal to length of pheno
# colnames.remove		names of columns of a to ignore
# MU				prior parameter, inverse expected number of regressors
# f.prior			g-prior parameter, e.g. 1 (AIC) or 1/length(pheno) (BIC)

lr.exhaustive = function(pheno, a, colnames.remove, MU, f.prior=1, nlevmax=50) {
	stopifnot(nrow(a)==length(pheno))
#             samp.size is the sample size
	samp.size = length(pheno)
#             xi.prior is a compound parameter depending on samp.size
#             and f.prior
	xi.prior = 1/(1+samp.size*f.prior)

	nlev = unlist(lapply(a,nlevels))
	colnums.remove = match(colnames.remove,colnames(a))
	stopifnot(!any(is.na(colnums.remove)))
	mvec.gd = setdiff(which(nlev!=1 & nlev<=nlevmax),colnums.remove)
#       Maximum number of regressors in the model, needed for the prior
	nmax = length(mvec.gd)

	# Matrix of all possible combinations of mvec.gd variables
	mvec.comb = matrix(0:1,ncol=1)
	for(j in 2:length(mvec.gd)) mvec.comb = cbind(mvec.comb[rep(1:nrow(mvec.comb),2),],rep(0:1,each=nrow(mvec.comb)))
	# Vectors of all possible models in vector and string form
	mvec = matrix(0,nrow(mvec.comb),ncol(a))
	string.mvec = rep("",nrow(mvec.comb))
	mvec2string = function(mvec) paste0(mvec,collapse="")
	for(i in 1:nrow(mvec.comb)) {
		mvec[i,mvec.gd[mvec.comb[i,]==1]] = 1
		string.mvec[i] = mvec2string(mvec[i,])
	}
	# Functions to compute log-likelihoods and log-priors
	logpr.prior = function(n,mu=MU,nmax=length(mvec.gd)) dbinom(n,nmax,1/mu/nmax,log=TRUE)-lchoose(nmax,n) # binomial distribution
	gen.calc.loglik = function() {
		fit0.logLik = logLik(glm(pheno ~ 1, family="binomial"))
		function(mvec) {
			wh = which(mvec==1)
			stopifnot(!any(nlev[wh]>nlevmax))
			if(length(wh)==0) {
				model = "pheno ~ 1"
			} else {
				model = paste("pheno ~",paste0("a[,",wh,"]",collapse=" + "))
			}
			fit = glm(as.formula(model), family="binomial")
			fit.logLik = logLik(fit)
			ret = fit.logLik-fit0.logLik
			attributes(ret)$df = attributes(fit.logLik)$df - attributes(fit0.logLik)$df
			# Do not allow singular model frames
			if(any(is.na(coef(fit)))) ret = ret + log(0)
			return(ret)
		}
	}; calc.loglik = gen.calc.loglik()
	calc.logBF = function(loglik) as.numeric(0.5*attributes(loglik)$df*log(xi.prior) + (1-xi.prior)*loglik)
	# Compute log-likelihoods and log-priors
	loglik = rep(NA,nrow(mvec.comb))
	bayes.loglik = rep(NA,nrow(mvec.comb))
	df = rep(NA,nrow(mvec.comb))
	logpr = rep(NA,nrow(mvec.comb))
	for(i in 1:nrow(mvec.comb)) {
		if(i==1) {
			loglik[i] = 0
			df[i] = 0
			bayes.loglik[i] = 0
		} else {
			tp = calc.loglik(mvec[i,])
			loglik[i] = tp
			df[i] = attributes(tp)$df
			bayes.loglik[i] = calc.logBF(tp)
		}
		logpr[i] = logpr.prior(sum(mvec[i,]),mu=MU)
		cat("Done",i,"of",nrow(mvec),"\n")
	}
	logpo = bayes.loglik + logpr; logpo = logpo - logsum(logpo)
	return(list(
		"MU"=MU,
		"f.prior"=f.prior,
		"samp.size"=samp.size,
		"xi.prior"=xi.prior,
		"nlevmax"=nlevmax,
		"nlev"=nlev,
		"nmax"=nmax,
		"mvec.gd"=mvec.gd,
		"mvec.comb"=mvec.comb,
		"mvec"=mvec,
		"string.mvec"=string.mvec,
		"logpr.prior"=logpr.prior,
		"calc.loglik"=calc.loglik,
		"calc.logBF"=calc.logBF,
		"loglik"=loglik,
		"bayes.loglik"=bayes.loglik,
		"df"=df,
		"logpr"=logpr,
		"logpo"=logpo
	))
}

# Get more HMP stats
hmp.stats.exhaustive = function(res, j, method="4b", BF0=NA, p0=NA) {
	j.incl.neg = which(res$mvec.comb[,j]==0)
	j.incl.pos = sapply(j.incl.neg,function(k) {
		tp = ifelse((1:ncol(res$mvec.comb))==j,1,res$mvec.comb[k,])
		wh = which(apply(t(res$mvec.comb)==tp,2,all))
		if(length(wh)!=1) return(NA)
		return(wh)
	})
	# exp(logpro) = log prior odds (pos vs neg)
	logpro = logsum(res$logpr[j.incl.pos]) - logsum(res$logpr[j.incl.neg])
	# exp(logmu) = Base model prior probabilities, appropriately normalized
	# Appears on the numerator and denominator
	logmu = res$logpr[j.incl.neg]+res$bayes.loglik[j.incl.neg]
	# exp(logp) = p-value
	# Appears on the numerator only
	logp = pchisq(2*(res$loglik[j.incl.pos]-res$loglik[j.incl.neg]),res$df[j.incl.pos]-res$df[j.incl.neg], low=FALSE,log=TRUE)
	# exp(logxi) = Multipliers to calculate the posterior odds from the p-value
	# Appears on the numerator only
	#   1. Standard HMP definition
	if(method=="1") {
		logxi = res$logpr[j.incl.pos] - res$logpr[j.incl.neg] + log(res$xi.prior)
	}
	#   2. Second HMP definition
	if(method=="2") {
		logxi = res$logpr[j.incl.pos] - res$logpr[j.incl.neg] + res$df[j.incl.neg]/2*log(res$xi.prior)
	}
	#   3. AIC like definition
	if(method=="3") {
		logxi = res$logpr[j.incl.pos] - res$logpr[j.incl.neg] - res$df[j.incl.neg]
	}
	#   4. Series-like expansion around exact point
	#      BF(p) ~ BF(p0) (p/p0)^(-alpha),            alpha ~ 1
	#   4a. 'Optimal' definition occurs when BF(p0)==BF(p) - in practice would induce
	#       an undesirable correlation between the weights and the inverse p-values
	if(method=="4a") {
		logxi = res$logpr[j.incl.pos] - res$logpr[j.incl.neg] + # prior model odds (pos vs neg)
		res$bayes.loglik[j.incl.pos] - res$bayes.loglik[j.incl.neg] + logp  # log[BF_0 * p_0]
	}
	#   4b. Define BF(p0)==Pr(M-)/Pr(M+), i.e. exact when posterior odds = 1
	#       This method aims to get the sign of the posterior odds correct
	if(method=="4b") {
		# Need p(BF0), the exact p-value corresponding to a particular Bayes factor
		# An intermediate step is R(BF0), the correponding LRT statistic
		# BF = xi^(df/2) * R^(1-xi)
		# (BF * xi^(-df/2))^(1/(1-xi)) = R
		delta.df = res$df[j.incl.pos]-res$df[j.incl.neg]
		logBF0 = -logpro
		logp0 = pchisq(2/(1-res$xi.prior)*(logBF0-delta.df/2*log(res$xi.prior)), delta.df,low=FALSE,log=TRUE)
		logxi = res$logpr[j.incl.pos] - res$logpr[j.incl.neg] + # prior model odds (pos vs neg)
		  logBF0 + logp0                                # log[BF_0 * p_0]
	}
	#   4c. Define BF(p0)==1
	#       This method aims to get the sign of the Bayes factor correct
	if(method=="4c") {
		# Need p(BF0), the exact p-value corresponding to a particular Bayes factor
		# An intermediate step is R(BF0), the correponding LRT statistic
		# BF = xi^(df/2) * R^(1-xi)
		# (BF * xi^(-df/2))^(1/(1-xi)) = R
		delta.df = res$df[j.incl.pos]-res$df[j.incl.neg]
		logBF0 = 0
		logp0 = pchisq(-delta.df/(1-res$xi.prior)*log(res$xi.prior),delta.df,low=FALSE,log=TRUE)
		logxi = res$logpr[j.incl.pos] - res$logpr[j.incl.neg] + # prior model odds (pos vs neg)
		  logBF0 + logp0                                # log[BF_0 * p_0]
	}
	#   4d. This method specifies BF0 as an argument
	#       Could warn if the actual gradient around p0 is not near -1
	if(method=="4d") {
		stopifnot(!is.na(BF0))
		# (BF * xi^(-df/2))^(1/(1-xi)) = R
		delta.df = res$df[j.incl.pos]-res$df[j.incl.neg]
		logBF0 = log(BF0)
		logp0 = pchisq(-delta.df/(1-res$xi.prior)*log(res$xi.prior),delta.df,low=FALSE,log=TRUE)
		logxi = res$logpr[j.incl.pos] - res$logpr[j.incl.neg] + # prior model odds (pos vs neg)
		  logBF0 + logp0                                # log[BF_0 * p_0]
	}
	#   4e. This method specifies p0 as an argument
	#       Could warn if the actual gradient around p0 is not near -1
	if(method=="4e") {
		stopifnot(!is.na(p0))
		# BF = xi^(df/2) * R^(1-xi)
		delta.df = res$df[j.incl.pos]-res$df[j.incl.neg]
		logp0 = log(p0)
		R = qchisq(logp0,delta.df,low=FALSE,log.p=TRUE)
		logBF0 = (delta.df)/2*log(res$xi.prior) + (1-res$xi.prior)*log(R)
		logxi = res$logpr[j.incl.pos] - res$logpr[j.incl.neg] + # prior model odds (pos vs neg)
		  logBF0 + logp0                                # log[BF_0 * p_0]
	}
	# For the HMP computation, normalize everything except logp
	logw = logmu + logxi - logsum(logmu + logxi)
	# Compute the HMP
	loghmp = -logsum(logw - logp)
	# Compute the approximate BF
	approx.logbf = logsum(logmu + logxi - logp) - logsum(logmu) - logpro
	# Factor to convert the HMP into the approximate BF
	logxisum = approx.logbf + loghmp
	#        = logsum(logmu + logxi - logp) - logsum(logmu) - logpro -logsum(logmu + logxi - logsum(logmu + logxi) - logp)
	#        = logsum(logmu + logxi) - logsum(logmu) - logpro
	if(FALSE) approx.logbf = logsum(
	  logpr[j.incl.neg] + bayes.loglik[j.incl.neg] - # Prior model weight (numerator)
	  logsum(logpr[j.incl.neg] + bayes.loglik[j.incl.neg]) +  # (denominator)
	  logpr[j.incl.pos] - logpr[j.incl.neg] +        # Prior model adjustment
	  bayes.loglik[j.incl.pos] - bayes.loglik[j.incl.neg] -  # marginal likelihood ratio
	  logsum(logpr[j.incl.pos]) + logsum(logpr[j.incl.neg]))  # Remove prior odds
	  # Exact BF and posterior odds
	  exact.logPoO = logsum(res$logpr[j.incl.pos] + res$bayes.loglik[j.incl.pos]) - logsum(res$logpr[j.incl.neg] + res$bayes.loglik[j.incl.neg])
	  exact.logBF = exact.logPoO - logsum(res$logpr[j.incl.pos]) + logsum(res$logpr[j.incl.neg])
	  c("HMP"=exp(loghmp),"thresh.approx.PoO"=exp(logxisum+logpro),"approx.PoO"=exp(logpro+approx.logbf),
	    #"alpha.PoP"=pharmonicmeanp(exp(logxisum+logpro),nrow(mvec.comb)/2),
		"thresh.approx.BF"=exp(logxisum), "approx.BF"=exp(approx.logbf), "exact.PoO"=exp(exact.logPoO), "exact.BF"=exp(exact.logBF))
}

# Workarounds
pharm_workaround = Vectorize(function(p,logL,log=FALSE) {
	if(is.na(p)) return(NA)
	tryCatch(pharmonicmeanp_workaround(p,logL,log=log), error=function(e) NA)
})

pharmonicmeanp_workaround = Vectorize(function(x, logL, log=FALSE, lower.tail=TRUE) {
    return(pLandau(1/x, mu=logL+1+psigamma(1)-log(2/pi), sigma=pi/2, log=log, lower.tail=!lower.tail))
})

