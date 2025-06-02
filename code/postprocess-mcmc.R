################
# Post-process #
################
# R code to merge MCMC chains and produce final results
#
###################

help = paste(
"Usage: Rscript postprocess-covid.R config_file",
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
results.posterior.inclusion.probs = "results.posterior-inclusion-probs.tsv"
results.group.posteriors.inclusion.probs = "results.group-posterior-inclusion-probs.tsv"
results.matrix.hmp.stats = "results.matrix-hmp-stats.tsv"
results.top.log.posteriors = "results.top-log-posteriors.tsv"
results.top.log.hmp = "results.top-log-hmp.tsv"
results.posterior.inclusion.cor = "results.posterior-inclusion-cov.tsv"
results.top.coefficients = "results.top-coefficients.tsv"

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

#################
# Merge results #
#################

# Merge posterior probabilities
mcmc.PoP = calc.PoP.merge(rds)

# Write to table
write.table(mcmc.PoP,results.posterior.inclusion.probs,row=T,col=T,quote=F,sep="\t")

# Define the top hits (or the regressors of interest)

# I changed this so all covariates get an output instead of only top 10. We need this for getting the correlations of groups later
top = which(p2o(mcmc.PoP[,"PoP"])>=0.00)
all =  which(p2o(mcmc.PoP[,"PoP"])>=0.00)

# Reprocess missing summaries
is.summarized = sapply(1:length(rds),function(i) !is.na(match(as.character(top),names(rds[[i]]$beta.sim))))
rownames(is.summarized) = as.character(top)
if(any(!is.summarized)) print("Reprocessing chains to summarize missing regressors")

# Iterate over chains
for(i in 1:length(rds)) {
	wh = which(!is.summarized[,i])
	if(length(wh)>0) {
		tbeg = Sys.time()
		chain.i = readRDS(chainfiles[i])
		for(TOP in top[wh]) {
			# See summarize.mcmc()
			cTOP = as.character(TOP)
			rds[[i]]$beta.sim[[cTOP]] = simulate.beta.chain(chain.i, TOP, burnin=rds[[i]]$burnin, thinning.interval=rds[[i]]$thinning)
			if(version==1) rds[[i]]$hmp.stats[[cTOP]] = hmp.stats.chain(chain.i, TOP)
			if(version==2) rds[[i]]$hmp.stats[[cTOP]] = hmp.stats.chain2(chain.i, TOP)
		}
		tdelta = Sys.time()-tbeg
		print(paste("Reprocessed chain",i,"in",dp(tdelta,2),attributes(tdelta)$units))
	}
}

# Merge HMP statistics for the top hits
hmp.stats = list()
for(j in top) {
	list.hmp.stats = list()
	for(i in 1:length(rds)) {
		tp = rds[[i]]$hmp.stats[[as.character(j)]]
		if(!is.null(tp)) list.hmp.stats[[as.character(i)]] = tp
	}
	print(paste("Variable",j,regressor.names[j]))
	if(version==1) hmp.stats[[regressor.names[j]]] = hmp.merge(list.hmp.stats,logL=L*log(2))
	if(version==2) hmp.stats[[regressor.names[j]]] = hmp.merge2(list.hmp.stats,logL=L*log(2))
}
mat.hmp.stats = t(sapply(1:length(top),function(i) hmp.stats[[i]]))
rownames(mat.hmp.stats) = regressor.names[top]

# Write to table
write.table(mat.hmp.stats,results.matrix.hmp.stats,row=T,col=T,quote=F,sep="\t")

# Prior probability of including a specific variable in any model
stopifnot(L==chain$nmax)
PrP = sum(exp(chain$logpr.prior(0:L)+lchoose(L,0:L))*(0:L)/L)/sum(exp(chain$logpr.prior(0:L)+lchoose(L,0:L)))

# Top hits table
top.od = order(-mat.hmp.stats[,"log.PoO.pe"])
if(version==1) {
	tophits = cbind("(prior)"=log10(p2o(PrP)), "mcmc.lo"=log10(p2o(mcmc.PoP[top,"PoP"]-2*mcmc.PoP[top,"SE"])),"mcmc"=log10(p2o(mcmc.PoP[top,"PoP"])),"mcmc.hi"=log10(p2o(mcmc.PoP[top,"PoP"]+2*mcmc.PoP[top,"SE"])),  "imsa.lo"=log10(exp(1))*(mat.hmp.stats[,"log.PoO.pe"]-2*mat.hmp.stats[,"log.PoO.stderr"]), "imsa"=log10(exp(1))*mat.hmp.stats[,"log.PoO.pe"], "imsa.hi"=log10(exp(1))*(mat.hmp.stats[,"log.PoO.pe"]+2*mat.hmp.stats[,"log.PoO.stderr"]))
	consistent.PoO =
		(tophits[,"mcmc.lo"]<=tophits[,"imsa"] & tophits[,"imsa"]<=tophits[,"mcmc.hi"]) |
		k
		(tophits[,"imsa.lo"]<=tophits[,"mcmc"] & tophits[,"mcmc"]<=tophits[,"imsa.hi"])
	rownames(tophits) = rownames(mat.hmp.stats)
	print("Table of log10 Posterior Odds:")
	data.frame(dp(tophits,2),"warn"=c("!","")[1+consistent.PoO],check.names=F)[top.od,]
} else if(version==2) {
	tophits = cbind("(prior)"=log10(p2o(PrP)), "mcmc"=log10(exp(1))*mat.hmp.stats[,"log.PoO.pe"], "stderr"=log10(exp(1))*mat.hmp.stats[,"log.PoO.stderr"])
	consistent.PoO =
	abs(mat.hmp.stats[,"log.PoO.pe"]-mat.hmp.stats[,"log.PoO.med"])<=2*mat.hmp.stats[,"log.PoO.stderr"] & abs(mat.hmp.stats[,"log.raw.hmp.pe"]-mat.hmp.stats[,"log.raw.hmp.med"])<=2*mat.hmp.stats[,"log.raw.hmp.stderr"]
	rownames(tophits) = rownames(mat.hmp.stats)
	print("Table of log10 Posterior Odds:")
	data.frame(dp(tophits,2),"warn"=c("!","")[1+consistent.PoO],check.names=F)[top.od,]
}

# Write table
rownames(tophits) = rownames(mat.hmp.stats)
write.table(tophits[top.od,],results.top.log.posteriors,row=T,col=T,quote=F,sep="\t")

# Give MCMC estimates of the HMP p-value
#pharm = Vectorize(function(p,l) {
#	if(is.na(p)) return(NA)
#	tryCatch(pharmonicmeanp(p,l), error=function(e) NA)
#})
tp.hmp = dp(cbind("log10.hmp"=log10(pharm_workaround(exp(mat.hmp.stats[,"log.raw.hmp.pe"]),L*log(2))),"stderr"=log10(exp(1))*(mat.hmp.stats[,"log.raw.hmp.stderr"]),"Bayes.thresh"=log10(pharm_workaround(exp(mat.hmp.stats[,"log.raw.hmp.bayes.thresh"]),L*log(2))),"Bonf.5%thresh"=log10(0.05/L)),2)
nlev2df = function(df) ifelse(df==0,1,df-1)
tp.hmp = data.frame("df"=nlev2df(chain$nlev[chain$mvec.gd][top]),tp.hmp[,1:3]," "=c("","*")[1+(tp.hmp[,1]<tp.hmp[,3])], "5%Bonf.thresh"=tp.hmp[,4]," "=c("","*")[1+(tp.hmp[,1]<tp.hmp[,4])],check.names=F)
rownames(tp.hmp) = rownames(mat.hmp.stats)
print("Table of Harmonic Mean p-Values:")
tp.hmp[top.od,]

# Write table
rownames(tp.hmp) = rownames(mat.hmp.stats)
write.table(tp.hmp[top.od,],results.top.log.hmp,row=T,col=T,quote=F,sep="\t")

# Merge posterior correlations in inclusion probabilities between regressors
mcmc.cooc = calc.regressor.cooccupancy.merge(rds)

# Output highest correlations to screen
wh = which(abs(mcmc.cooc*lower.tri(mcmc.cooc))>0.3,arr=TRUE)
print(paste("Found",length(wh),"pairs of regressors with |cor|>0.3"))
data.frame("rgr1"=substr(regressor.names[wh[,1]],1,30),"rgr2"=substr(regressor.names[wh[,2]],1,30),"cor"=dp(mcmc.cooc[wh],2))

# Write table
write.table(mcmc.cooc,results.posterior.inclusion.cor,row=T,col=T,quote=F,sep="\t")

# Posterior parameter estimates
list.mcmc.beta = lapply(all,function(j) summarize.beta.sim.merge(rds, as.character(j)))
mcmc.beta = matrix(0,0,3)
colnames(mcmc.beta) = names(list.mcmc.beta[[1]])
for(i in top.od) {
	mcmc.beta.elem = cbind(list.mcmc.beta[[i]][[1]],list.mcmc.beta[[i]][[2]],list.mcmc.beta[[i]][[3]])
	if(nrow(mcmc.beta.elem)==1) {
		rownames(mcmc.beta.elem) = regressor.names[all][i]
	} else {
		rownames(mcmc.beta.elem) = paste0(regressor.names[all][i],":",1:nrow(mcmc.beta.elem))
	}
	mcmc.beta = rbind(mcmc.beta,mcmc.beta.elem)
}

tp = mcmc.beta
rownames(tp) = ifelse(nchar(rownames(mcmc.beta))>32, paste0(substr(rownames(mcmc.beta),1,25),"...",substr(rownames(mcmc.beta),nchar(rownames(mcmc.beta))-4,nchar(rownames(mcmc.beta)))), rownames(mcmc.beta))
print("Estimates of Coefficients:")
dp(tp,2)

# Write table
write.table(mcmc.beta,results.top.coefficients,row=T,col=T,quote=F,sep="\t")

#####################
# Define the groups #
#####################

chainfiles = sort(dir(".",glob2rx("doublethink-mcmc-out.ID*.RDS")))

E = mcmc.PoP$PoP; names(E) <- rownames(mcmc.PoP)
mcmc.cooc = read.delim(results.posterior.inclusion.cor,as.is=TRUE,check.names=FALSE)
COR = as.matrix(mcmc.cooc); # To save memory: rm(mcmc.cooc)

### HACK: fill in NAs assuming zero correlation
diag.is.na = is.na(diag(COR))
diag(COR)[diag.is.na] <- 1
offdiag.is.na = is.na(COR[lower.tri(COR)])
COR[lower.tri(COR)][offdiag.is.na] <- 0
offdiag.is.na = is.na(COR[upper.tri(COR)])
COR[upper.tri(COR)][offdiag.is.na] <- 0
stopifnot(!any(is.na(COR)))

# Sanity checks
stopifnot(nrow(COR)==ncol(COR))
stopifnot(rownames(COR)==names(E))

# Define the posterior variance in inclusion probabilities
V = E*(1-E)
# Define the posterior correlation in inclusion probabilities
COR.DENOM = outer(sqrt(V),sqrt(V))
COV = COR * COR.DENOM
# Define a positive number close to zero (not used)
# eps = 1e-6
# Compute the probability both variables are included (not used)
# EE = COV + outer(E,E); EE[EE<eps] = 0


# Group everything with strong negative correlation
#gp = 1:nrow(COR)
#for(i in 1:nrow(COR)) {
#  # Change the rule to depend on COV not COR
#  #wh = !is.na(COR[i,]) & COR[i,]< -0.25
#  #wh = !is.na(COV[i,]) & COV[i,]< -0.05
#  # Or both
#  wh = (!is.na(COR[i,]) & COR[i,]< -0.25) &
#	   (!is.na(COV[i,]) & COV[i,]< -0.01)
#  gp[wh] = gp[i]
#}
#ct.gp = table(factor(gp,levels=1:nrow(COR)))
# Create factors for groups
system(paste0(python_loc, " ", sourcedir,"/cluster_posterior_probabilities.py"))
gp = as.integer(readLines("results.clusters.out"))
fgp = factor(gp, levels=0:max(gp))
# Ensure gp is a vector of length ncol(mvec) with integer group
# assignments between 1 and ngp, where ngp is the number of groups
gp = as.integer(fgp)
ngp = max(gp, na.rm=TRUE)
stopifnot(nlevels(factor(gp))==ngp)

# For each group thus defined, calculate maximum sum across PIPs
maxE.group = sapply(1:ngp,function(GP) sum(E[!is.na(gp) & gp==GP]))
# Use delimiter \037 (ASCII 31) known as Unit [i.e. field] Separator (US)
who.group = sapply(1:ngp,function(GP) paste0(rownames(COR)[!is.na(gp) & gp==GP],collapse="\037"))
names(maxE.group) <- who.group
head(sort(maxE.group,decreasing=TRUE))

# Create a binary matrix where mgp[i,j] indicates whether variable i
# is a member of group j (1) or not (0).
mgp = 1*outer(gp, 1:ngp, function(x, y) ifelse(is.na(x==y), FALSE, x==y))

##############################################
# Reprocess the chains to compute group PoPs #
##############################################
string2mvec = function(s) as.integer(substring(s,1:nchar(s),1:nchar(s)))

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

# Repeat for the remaining chains
for(i in 2:length(chainfiles)) {
	chain.i = readRDS(chainfiles[i])
	string.mvec = with(chain.i, names(loglik.rec)[state])[-(1:(burnin+1))]
	mvec = t(sapply(string.mvec, string2mvec))
	mvec.gp = mvec[,mvec.gd] %*% mgp
	# mcmc.PoP.gp[i,j] is the posterior inclusion probability of group j in chain i
	# This enables calculation of a standard error for mcmc.PoP.gp across chains
	mcmc.PoP.gp = rbind(mcmc.PoP.gp, colMeans(mvec.gp>0))
	print(paste("Reprocessed chain",i))
}

colnames(mcmc.PoP.gp) <- who.group
out_means = sort(colMeans(mcmc.PoP.gp), decreasing=TRUE)
#write.table(out_means, file=results.group.posteriors.inclusion.probs, quote=F,sep="\t", col.names=c("group_name\tcombined_PoP"))
if(length(out_means)==1) write.table(out_means, file=results.group.posteriors.inclusion.probs, quote=F,sep="\t", col.names=c("group_name\tcombined_PoP"))
if(length(out_means)>1) write.table(out_means[2:length(out_means)], file=results.group.posteriors.inclusion.probs, quote=F,sep="\t", col.names=c("group_name\tcombined_PoP"))
system(paste0(python_loc, " ", sourcedir, "/postpostprocessing.py"))

# Compile main table
system(paste0(sourcedir, "/compile-main-table.R ", args[1]))
