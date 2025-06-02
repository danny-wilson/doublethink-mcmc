# Function to convert NAs to zeros
NAto0 = function(x) ifelse(x<0,0,x)
# Read the processed phenotype data
pheno20201231 = read.delim("/fullpath/ukb41482.bd.gwasdata.DOUBLETHINK-2020-12-31.txt", check.names=FALSE, stringsAsFactors=FALSE)

# Output phenotype A2
write.table(data.frame("eid" = pheno20201231$eid, "pheno" = pheno20201231$A2), row=FALSE, col=TRUE, sep=",", quote=FALSE, file="/fullpath/A2-2020-12-31.csv")
write(pheno20201231$eid[pheno20201231$A2<0], "/fullpath/exclude_eids.A2-2020-12-31.csv", ncol=1)

# Output phenotype B2
write.table(data.frame("eid" = pheno20201231$eid, "pheno" = pheno20201231$B2), row=FALSE, col=TRUE, sep=",", quote=FALSE, file="/fullpath/B2-2020-12-31.csv")
write(pheno20201231$eid[pheno20201231$B2<0], "/fullpath/exclude_eids.B2-2020-12-31.csv", ncol=1)

# Output phenotype C2
write.table(data.frame("eid" = pheno20201231$eid, "pheno" = pheno20201231$C2), row=FALSE, col=TRUE, sep=",", quote=FALSE, file="/fullpath/C2-2020-12-31.csv")
write(pheno20201231$eid[pheno20201231$C2<0], "/fullpath/exclude_eids.C2-2020-12-31.csv", ncol=1)
