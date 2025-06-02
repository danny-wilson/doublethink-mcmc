#!/usr/bin/env Rscript

######################
# Compile main table #
######################
# Merge output files for
# main table
######################

help = paste(
"Usage: Rscript compile-main-table.R config_file",
sep="\n")

# Functions
null2NA = function(x) ifelse(is.null(x), NA, x)
NA2dash = function(x) ifelse(is.na(x), "-", x)
NA20 = function(x) ifelse(is.na(x), 0, x)

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
infile.top.coefficients = "results.top-coefficients.tsv"
infile.group.posterior.inclusion.probs = "results.group-posterior-inclusion-probs.tsv"
#infile.group.posterior.inclusion.probs.prior.groups = "results.group-posterior-inclusion-probs.prior-groups.tsv"

results.table.variables = "results.table.variables.tsv"
results.table.prior.groups = "results.table.prior.groups.tsv"
results.table.post.groups = "results.table.post.groups.tsv"

# Read files

posterior.inclusion.probs = read.delim(infile.posterior.inclusion.probs, check.names=FALSE)
top.coefficients = read.delim(infile.top.coefficients, check.names=FALSE); top.coefficients = top.coefficients[match(rownames(posterior.inclusion.probs), rownames(top.coefficients)),]
group.posterior.inclusion.probs = read.delim(infile.group.posterior.inclusion.probs, check.names=FALSE)
#group.posterior.inclusion.probs.prior.groups = read.delim(infile.group.posterior.inclusion.probs.prior.groups, check.names=FALSE)

# Work out sample size
# See postpostprocessing.py
pheno = read.csv(cases_filename)
number_cases = sum(pheno$pheno==1)
if(!exists("controls_subsample")) controls_subsample = 0
if(controls_subsample!=0) {
	number_ctrls = min(sum(pheno$pheno==0), controls_subsample)
} else {
	number_ctrls = sum(pheno$pheno==0)
}

# Paper parameterizations from software parameterizations:
# See "Convert doublethink parameterizations.txt"
paper.nu = nrow(posterior.inclusion.probs)
stopifnot(paper.nu == nrow(infile.top.coefficients))
paper.mu = mu_prior
paper.h = h_prior
paper.n = number_cases + number_ctrls
paper.xi = paper.h/(paper.n + paper.h)

# Create a table for individual variables
format.pop = function(PoP) round(1000*PoP)/10
#format.poo = function(PoP, maxexp=6) round(10*log10(pmax(10^-maxexp, pmin(10^maxexp, PoP/(1-PoP)))))/10
format.poo = function(PoP, maxexp=6) round(100*log10(PoP/(1-PoP)))/100
cutoff0.02 = function(x) ifelse(x> -log10(0.02), x, "-")
calc.neglog10p = function(PoP) cutoff0.02(round(100*-log10(pchisq(2*log(PoP/(1-PoP)/paper.nu/paper.mu/sqrt(paper.xi)), 1, low=FALSE)))/100)
format.beta = function(x) round(100*x)/100
table.variables = data.frame("Variable"=rownames(posterior.inclusion.probs), "Posterior inclusion probability (%)"=format.pop(posterior.inclusion.probs$PoP), "Monte Carlo SE (%)"=format.pop(10*posterior.inclusion.probs$SE)/10, "log10 Posterior odds"=format.poo(posterior.inclusion.probs$PoP), "-log10 Adjusted p-value"=calc.neglog10p(posterior.inclusion.probs$PoP), "Effect size when included"=format.beta(top.coefficients$non.zero.mean), "Standard error when included"=format.beta(top.coefficients$non.zero.sd), check.names=FALSE)
table.variables.od = order(posterior.inclusion.probs$PoP, decreasing=TRUE) #[1:sum(table.variables$"-log10 Adjusted p-value"!="-")]
out.table.variables = table.variables[table.variables.od,]

if(FALSE) {
	# Create a table for groups defined a priori by variable correlation
	group.posterior.inclusion.probs.prior.groups$"Group #" = 1:nrow(group.posterior.inclusion.probs.prior.groups)
	table.variables.prior.group = sapply(table.variables$Variable, function(pattern) null2NA(grep(pattern, group.posterior.inclusion.probs.prior.groups$Group, fixed=TRUE)))
	table.variables.prior.group.pop = posterior.inclusion.probs$PoP; table.variables.prior.group.pop[!is.na(table.variables.prior.group)] = group.posterior.inclusion.probs.prior.groups$PoP[table.variables.prior.group[!is.na(table.variables.prior.group)]]
	table.prior.groups = data.frame("Group"=NA2dash(table.variables.prior.group), "Group posterior inclusion probability (%)"=format.pop(table.variables.prior.group.pop), "Group log10 Posterior odds"=format.poo(table.variables.prior.group.pop), "Group -log10 Adjusted p-value"=calc.neglog10p(table.variables.prior.group.pop), table.variables, check.names=FALSE)
	table.prior.groups.od = order(table.variables.prior.group.pop - 1e-6*NA20(table.variables.prior.group), posterior.inclusion.probs$PoP, decreasing=TRUE) #[1:sum(table.prior.groups$"Group -log10 Adjusted p-value"!="-")]
	out.table.prior.groups = table.prior.groups[table.prior.groups.od,]
}

# Create a table for groups defined a posteriori by variable inclusion correlation
group.posterior.inclusion.probs$"Group #" = 1:nrow(group.posterior.inclusion.probs)
table.variables.post.group = sapply(table.variables$Variable, function(pattern) null2NA(grep(pattern, group.posterior.inclusion.probs$group_name, fixed=TRUE)))
table.variables.post.group.pop = posterior.inclusion.probs$PoP; table.variables.post.group.pop[!is.na(table.variables.post.group)] = group.posterior.inclusion.probs$combined_PoP[table.variables.post.group[!is.na(table.variables.post.group)]]
table.post.groups = data.frame("Group"=NA2dash(table.variables.post.group), "Group posterior inclusion probability (%)"=format.pop(table.variables.post.group.pop), "Group log10 Posterior odds"=format.poo(table.variables.post.group.pop), "Group -log10 Adjusted p-value"=calc.neglog10p(table.variables.post.group.pop), table.variables, check.names=FALSE)
table.post.groups.od = order(table.variables.post.group.pop - 1e-6*NA20(table.variables.post.group), posterior.inclusion.probs$PoP, decreasing=TRUE) #[1:sum(table.post.groups$"Group -log10 Adjusted p-value"!="-")]
out.table.post.groups = table.post.groups[table.post.groups.od,]

# Output results
write.table(out.table.variables, file=results.table.variables, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
#write.table(out.table.prior.groups, file=results.table.prior.groups, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
write.table(out.table.post.groups, file=results.table.post.groups, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
