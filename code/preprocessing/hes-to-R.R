# Create a data frame whose row names are EIDs, columns ICD10 codes and
# entries are the first admission date associated with diagnosis (primary
# or secondary) of the ICD10 code in the participant with that EID, or NA.

# screen -x 15689
# module add R/4.2.1-foss-2022a
# R
# Based on hgi-phenotypes-v2.R and reprocess-raw-hes.R

#############
# Functions #
#############
anymatch = function(x, y, ...) !is.na(match(x,c(y,...)))

#############
# Constants #
#############
# Output file name prefix
out.prefix = "/fullpath/ukb41482.hes"

# Hospital episode statistics
# November 2021 update
hesin.filename = "/fullpath/hesin.17November2021.txt.gz"
hesin_diag.filename = "/fullpath/hesin_diag.19November2021.txt.gz"

# EIDs and Sample size
system.time((eid = as.vector(as.matrix(read.delim(pipe("cut -f1 /fullpath/ukb41482.tab"), colClasses="integer")))))
#  user  system elapsed
#73.336  13.053  86.393

(n = length(eid))
#[1] 502505

# Filter out ICD10 codes at a frequency less than this
# Deliberately lenient: allows more stringency in downstream analysis
(min.icd10.freq = 0)

###############################
# Hospital episode statistics #
###############################
system.time((hesin = read.delim(hesin.filename,as.is=TRUE)))
#  user  system elapsed
#29.126   0.863  30.019

system.time((hesin_diag = read.delim(hesin_diag.filename,as.is=TRUE)))
#  user  system elapsed
#15.452   0.521  16.039

# Format dates correctly
hesin$epistart = as.Date(hesin$epistart,tryFormats=c("%d/%m/%Y","%Y-%m-%d"))
hesin$epiend = as.Date(hesin$epiend,tryFormats=c("%d/%m/%Y","%Y-%m-%d"))
hesin$elecdate = as.Date(hesin$elecdate,tryFormats=c("%d/%m/%Y","%Y-%m-%d"))
hesin$admidate = as.Date(hesin$admidate,tryFormats=c("%d/%m/%Y","%Y-%m-%d"))
hesin$disdate = as.Date(hesin$disdate,tryFormats=c("%d/%m/%Y","%Y-%m-%d"))

# Merge, filtering on non-NA
hesin.gd = !is.na(hesin$admidate)
hesin_diag.gd = !is.na(hesin_diag$diag_icd10)
# NB: need to make all=FALSE to impose filtering criteria above
system.time((hesin_merge = merge(hesin[hesin.gd, ], hesin_diag[hesin_diag.gd, ], by=c("eid", "ins_index"), all=FALSE, stringsAsFactors=FALSE)))
#  user  system elapsed
#57.702   3.180  60.881

dim(hesin_merge)
# [1] 13786843       48

# All observed ICD10 codes
lev.icd10 = setdiff(sort(unique(hesin_merge$diag_icd10)),"")
length(lev.icd10)
# [1] 12097

# De-duplicate by eid (participant ID)
hesin_merge$eid_diag_icd10 = paste0(hesin_merge$eid, "_", hesin_merge$diag_icd10)
hesin_merge$eid_diag_icd10_dup = duplicated(hesin_merge$eid_diag_icd10)

# De-duplicated counts by ICD10 code
tb.lev.icd10 = table(factor(hesin_merge$diag_icd10[!hesin_merge$eid_diag_icd10_dup], levels=lev.icd10))
# Number of codes remaining after minimum frequency filter
sum(tb.lev.icd10 >= n*min.icd10.freq)
#[1] 12097

# Read in descriptions
icd10code = read.delim("/fullpath/ukb-icd10-data-coding-19.tsv", stringsAsFactors=FALSE)
desc.icd10 = icd10code$meaning[match(lev.icd10, icd10code$coding)]
# Check there are no unrecognised codes
lev.icd10[is.na(desc.icd10)]
#character(0)

# Apply filter on number of observed participants with code
lev.icd10.gd = tb.lev.icd10/n >= min.icd10.freq
table(lev.icd10.gd)
# TRUE
#12097

# Build the data frame containing the first admission date for each ICD10 code
system.time((bd.hes = data.frame(row.names=eid)))
for(LEV in lev.icd10[lev.icd10.gd]) {
	# Subset the entries corresponding to the specified ICD10 code
	gd = hesin_merge$diag_icd10==LEV
	# Order the subset in ascending date
	od = order(hesin_merge$admidate[gd], decreasing=FALSE)
	# Match the first instance of the date-ordered ICD10 codes (or NA if no instance)
	wh = hesin_merge$admidate[gd][od][match(eid, hesin_merge$eid[gd][od])]
	# Add the result to the data frame
	bd.hes[LEV] <- as.Date(wh)
	cat("Done",LEV,"\n")
}

# Find the maximum date of any entry: use this to name the file
(maxdate = max(hesin_merge$admidate,na.rm=TRUE))
#[1] "2021-03-31"
out.filename = paste0(out.prefix,".maxdate.",maxdate,".rds")
system.time((saveRDS(bd.hes, out.filename)))
#   user  system elapsed
#260.331  18.810 279.772
