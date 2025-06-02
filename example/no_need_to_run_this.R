# No need to run this - it simply documents the provenance and the
# preparation of the example data for doublethink-mcmc

# Load R object 'amd_example' from remote .RData file at https://github.com/verena-zuber/demo_AMD
load(url("https://github.com/verena-zuber/demo_AMD/raw/refs/heads/master/amd_example"))

# Process the data
# Estimates of the effects of 148 markers on levels of 49 biomarkers
betaX = amd_example$betaX
# Estimates of the effects of 148 markers on age-related macular degeneration
amd_beta = amd_example$amd_beta
# Standard errors of the effects of 148 markers on age-related macular degeneration
amd_se = amd_example$amd_se
# Reference SNP cluster IDs of the 148 markers
rs = amd_example$annotate[,1]
# Genes in which the SNPs are situated
genes = amd_example$annotate[,7]

# Inverse variance weighting (IVW) based on the standard error of the amd_beta effects prior to subsequent analysis
# Scaled estimates of the effects of 148 markers on levels of 49 biomarkers
betaX_ivw = betaX / amd_se
# Scaled estimates of the effects of 148 markers on age-related macular degeneration
amd_beta_ivw = matrix(amd_beta / amd_se, ncol=1)
# Names of the 49 biomarkers
biomarkers = gsub("^beta_", "", colnames(betaX_ivw))
# Append the rs numbers to both matrices
betaX_ivw = cbind("rs" = as.numeric(gsub("rs", "", rs)), betaX_ivw)
amd_beta_ivw = cbind("rs" = as.numeric(gsub("rs", "", rs)), amd_beta_ivw)
# Name the columns of the features (betaX_ivw) and outcomes (amd_beta_ivw)
colnames(amd_beta_ivw) <- c("eid", "pheno")
colnames(betaX_ivw) <- c("eid", biomarkers)

# There was one influential variant in the LIPC gene region and two outliers in the APOE and FUT2 gene region
# Remove these 3 data points for the analysis
LIPC = which(genes == "LIPC")
FUT2 = which(genes == "FUT2")
APOE = which(genes == "APOE")
exclude_vec = c(LIPC,FUT2,APOE)

# Create the example data directory
dir.create("~/Documents/GitHub/doublethink-mcmc/example", recursive=TRUE)
# Save the features (betaX_ivw)
write.csv(betaX_ivw, "~/Documents/GitHub/doublethink-mcmc/example/features.csv")
# Save the outcomes (amd_beta_ivw)
write.csv(amd_beta_ivw, "~/Documents/GitHub/doublethink-mcmc/example/outcomes.csv", row.names=FALSE)
# Save a binary version of the outcomes
write.csv(cbind("eid"=amd_beta_ivw[,"eid"], "pheno"=(1+sign(amd_beta_ivw[,"pheno"]))/2), "~/Documents/GitHub/doublethink-mcmc/example/binary_outcomes.csv", row.names=FALSE)
# Save the column names and R data types (set all equal to double)
write.table(cbind(colnames(betaX_ivw), "double"), "~/Documents/GitHub/doublethink-mcmc/example/columns.csv", row.names=FALSE, col.names=FALSE, quote=TRUE, sep=",")
# Columns to exclude
keep = c("Ace", "Ala", "ApoA1", "ApoB", "Gln", "His", "IDL.TG", "L.HDL.C", "LDL.D", "M.HDL.C", "S.HDL.TG", "S.VLDL.TG", "VLDL.D", "XL.HDL.C", "XL.HDL.TG")
write(setdiff(biomarkers, keep), "~/Documents/GitHub/doublethink-mcmc/example/exclude_columns.csv", ncol=1)
# Rows to exclude
write(amd_beta_ivw[exclude_vec, "eid"], "~/Documents/GitHub/doublethink-mcmc/example/exclude_rows.csv", ncol=1)
