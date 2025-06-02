#############
# bd-to-R.R #
#############
# This code creates ukb41482.rds from ukb41482.tab, a dataframe with correctly type-cast variables
# according to the UK Biobank Data Dictionary and Data Codings

# Run on dingo 5th May 2023
# screen -x 25415
# module add R/4.2.1-foss-2022a
# R

#############
# Functions #
#############
anymatch = function(x, y, ...) !is.na(match(x,c(y,...)))

########
# Data #
########
# Read the publicly available Data Dictionary, and manually add missing entries
ukbdesc = rbind(read.delim("/fullpath/Data_Dictionary_Showcase.tsv", as.is=TRUE, stringsAsFactors=FALSE, quote = ""),
c("Assessment centre > Recruitment > Reception", 100024, 20074, "Home location at assessment - east co-ordinate (1km resolution)", 497834, 550740, "Complete", "Integer", "metre-grid", "Data", "Derived", "Unisex", 4, 1, NA, "This location information is based on the postcode for the address to which the participants invitation was sent. If the participant updated their home address during their assessment centre visit, then that postcode was used as the source instead.  The location co-ordinates use the Ordnance Survey (OSGB) reference, rounded to the nearest Kilometre. This East co-ordinate needs to be used in conjunction with the North co-ordinate in Field 20075.  For most purposes, grid co-ordinates rounded to a 1Km distance will provide sufficient information on location of residence at assessment. Research proposals that require location details to a finer level of granularity (i.e. 100m) should request Field 20033 and Field 20034, however special approval will be required to release the more detailed information.", "https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20074"),
c("Assessment centre > Recruitment > Reception", 100024, 20075, "Home location at assessment - north co-ordinate (1km resolution)", 497834, 550740, "Complete", "Integer", "metre-grid", "Data", "Derived", "Unisex", 4, 1, NA, "This location information is based on the postcode for the address to which the participants invitation was sent. If the participant updated their home address during their assessment centre visit, then that postcode was used as the source instead.  The location co-ordinates use the Ordnance Survey (OSGB) reference, rounded to the nearest Kilometre. This North co-ordinate needs to be used in conjunction with the East co-ordinate in Field 20074.  For most purposes, grid co-ordinates rounded to a 1Km distance will provide sufficient information on location of residence at assessment. Research proposals that require location details to a finer level of granularity (i.e. 100m) should request Field 20033 and Field 20034, however special approval will be required to release the more detailed information.","https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20075")
)
# Add partial entries for the following
# eid:         Coded participant ID, unique per project
# 22100-22125: Coded links to genotype data
# 22400-22402: Imaging data
# 22411-22414: Imaging data
# 22800-22823: Coded links to imputed genotype data
ukbdesc.blank = ukbdesc[1,]; ukbdesc.blank[1,1:ncol(ukbdesc.blank)] <- NA; for(i in 1:56) ukbdesc.blank = rbind(ukbdesc.blank,ukbdesc.blank[1,])
ukbdesc.blank$FieldID = c('eid', '22100', '22101', '22102', '22103', '22104', '22105', '22106', '22107', '22108', '22109', '22110', '22111', '22112', '22113', '22114', '22115', '22116', '22117', '22118', '22119', '22120', '22121', '22122', '22123', '22124', '22125', '22400', '22401', '22402', '22411', '22412', '22414', '22800', '22801', '22802', '22803', '22804', '22805', '22806', '22807', '22808', '22809', '22810', '22811', '22812', '22813', '22814', '22815', '22816', '22817', '22818', '22819', '22820', '22821', '22822', '22823')
ukbdesc.blank$ValueType = c('Text', rep('Text', 26), rep('Continuous', 3), 'Categorical single', 'Categorical multiple', 'Categorical single', rep('Text', 24))
ukbdesc = rbind(ukbdesc, ukbdesc.blank)
ukbdesc$Coding[ukbdesc$FieldID=="22411"] = "403"
ukbdesc$Coding[ukbdesc$FieldID=="22412"] = "401"
ukbdesc$Coding[ukbdesc$FieldID=="22414"] = "400"

# Read the publicly available Data Codings, and manually add missing entries
ukbcode = rbind(read.delim("/fullpath/Codings.tsv", as.is=TRUE, stringsAsFactors=FALSE, quote=""),
c(439,"1900-01-01 00:00:00","Not performed"),
c(586,"1910-01-01","Prefer not to answer"),
c(586,"1920-01-01","Do not know"),
c(586,"1930-01-01","Not applicable"),
c(1313,"1904-04-04","Date in or before calendar year of birth"),
c(272,"1900-01-01","Date is unknown")
)

bd.colnames = colnames(read.table("/fullpath/ukb41482.tab", header=TRUE, sep="\t", as.is=TRUE, stringsAsFactors=FALSE, nrows=1))
bd.nrows = scan((PIPE=pipe("wc -l /fullpath/ukb41482.tab | cut -f1 -d' '")),what=integer(0))-1; close(PIPE)
system.time((bd = read.table("/fullpath/ukb41482.tab", header=TRUE, sep="\t", as.is=TRUE, stringsAsFactors=FALSE, colClasses=rep("character",length(bd.colnames)), nrows=bd.nrows)))
#    user   system  elapsed
#1314.301  124.801 1439.317

dim(bd)
#[1] 502505  18432

##############
# Processing #
##############
cn = sapply(colnames(bd),function(s) unlist(strsplit(s,".",fixed=TRUE))[2:4])
#table(cn[3,])

cn1.n = table(factor(cn[1,],levels=unique(cn[1,])))
length(cn1.n)
#[1] 3543

table(cn1.n)
#   1    2    3    4    5    6    7    8   10   11   12   13   14   15   16   17 
# 930 1369  166  461  303    9    5   28   12    1   68    2    2    9   10    8 
#  18   19   20   22   23   24   25   28   30   32   33   35   40   44   47   48 
#   1    3   17    3    1   12    9   10    4   10   10    2   19    1    2    1 
#  49   60   65   66   85   86  105  117  120  128  136  184  192  213  224  228 
#   2   20    1    2    1    1    1    6    1    5    5    1    1    2    1    5 

# Sanity check the field classes
bd2ukbdesc = match(cn[1,], ukbdesc$FieldID)

# No missing fields wrt ValueType
table(ukbdesc$ValueType[bd2ukbdesc], !is.na(ukbdesc$Coding[bd2ukbdesc]), useNA="a")
#					   FALSE TRUE <NA>
#  Categorical multiple     0 3111    0
#  Categorical single       0 4854    0
#  Compound                32    0    0
#  Continuous            3317  775    0
#  Date                   658   35    0
#  Integer               2988 1084    0
#  Text                  1372    0    0
#  Time                   178   28    0
#  <NA>                     0    0    0

# Explicitly type cast every column of bd.

# For debugging purposes, because re-reading the data is slow
options(warn=2)
dry.run = FALSE
# For every unique field ID in bd, cast to the correct class
for(CN1 in unique(cn[1,])) {
    # Find all columns matching this fieldID
    COLS = which(cn[1,]==CN1)
	# Obtain the field from showcase
	showcase.field = unique(ukbdesc$Field[bd2ukbdesc][COLS])
    # Obtain the class from showcase
    showcase.class = unique(ukbdesc$ValueType[bd2ukbdesc][COLS])
	# Obtain coding (if any) from showcase
	showcase.code = unique(ukbdesc$Coding[bd2ukbdesc][COLS])
    # Sanity checks
    stopifnot(length(showcase.class)==1)
	stopifnot(length(showcase.code)==1)
	# Extract the codings, if any
	field.code = ukbcode[0, c("Value", "Meaning")]
	if(!is.na(showcase.code)) field.code = ukbcode[ukbcode$Coding==showcase.code, c("Value", "Meaning")]
    # Manually type cast the columns as appropriate
    if(is.na(showcase.class)) {
		stop("ValueType missing for some FieldIDs")
    } else if(showcase.class=='Categorical multiple' | showcase.class=='Categorical single') {
		# Levels observed across all instances of the fieldID
		LEVS = sort(setdiff(unique(unlist(apply(bd[,COLS,drop=FALSE],2,function(x) as.character(unique(x))))),NA))
		# Warn about fields with many levels
		if(length(LEVS)>1000) {
			cat(paste(CN1,showcase.field,"has",length(LEVS),"levels\n"))
		}
		# Identify the expected levels and labels
		field.levs = as.character(field.code$Value)
		field.labs = paste(field.levs,":",as.character(field.code$Meaning))
		# Sanity check: observed values match expected coding
		if(!all(anymatch(LEVS, field.levs))) {
			cat(paste(CN1,showcase.field,"has unmatched levels:",paste(setdiff(LEVS,field.levs),collapse=" "),"\n"))
		}
		for(COL in COLS) {
			# Unmatched levels are converted to NA (even if sometimes these are obvious data entry typos)
			if(!dry.run) bd[,COL] <- factor(bd[,COL], levels=field.levs, labels=field.labs)
		}
    } else if(showcase.class=='Compound') {
		for(COL in COLS) {
			is.special = anymatch(bd[,COL], field.code$Value)
			if(!dry.run) {
				# Keep Compound variables as text
				bd[,COL] <- as.character(bd[,COL])
				# Replace special values with NA
				bd[is.special, COL] <- NA
			}
		}
	} else if(showcase.class=='Continuous') {
		for(COL in COLS) {
			is.special = anymatch(bd[,COL], field.code$Value)
			if(!dry.run) {
				# Transform Continuous variables to numeric
				bd[,COL] <- as.numeric(bd[,COL])
				# Replace special values with NA
				bd[is.special, COL] <- NA
			}
		}
	} else if(showcase.class=='Date') {
		for(COL in COLS) {
			is.special = anymatch(bd[,COL], field.code$Value)
			if(!dry.run) {
				# Transform Date variables to Date
				bd[,COL] <- as.Date(bd[,COL], "%Y-%m-%d")
				# Replace special values with NA
				bd[is.special, COL] <- NA
			}
		}
	} else if(showcase.class=='Time') {
		for(COL in COLS) {
			is.special = anymatch(bd[,COL], field.code$Value)
			if(!dry.run) {
				# Transform Time variables to POSIXct
				bd[,COL] <- as.POSIXct(bd[,COL], tz="GMT", format="%Y-%m-%d %H:%M:%OS")
				# Replace special values with NA
				bd[is.special, COL] <- NA
			}
		}
	} else if(showcase.class=='Integer') {
		for(COL in COLS) {
			is.special = anymatch(bd[,COL], field.code$Value)
			if(!dry.run) {
				# Transform Integer variables to integer
				bd[,COL] <- as.integer(bd[,COL])
				# Replace special values with NA
				bd[is.special, COL] <- NA
			}
		}
	} else if(showcase.class=='Text') {
		for(COL in COLS) {
			is.special = anymatch(bd[,COL], field.code$Value)
			if(!dry.run) {
				# Transform Text variables to character
				bd[,COL] <- as.character(bd[,COL])
				# Replace special values with NA
				bd[is.special, COL] <- NA
			}
		}
	} else {
		stop("Unexpected ValueType")
	}
	#cat("Done",CN1,"\n")
}
#393 Program (tactus) version ID (compiler timestamp) has unmatched levels: <TMar 7 2006 08:34:01> [TAug 4 2009 16:07:02] [TAug 8 2007 10:26:36] [TJul 5 2007 15:03:00]
#680 Own or rent accommodation lived in has unmatched levels: -1
#1190 Nap during day has unmatched levels: -1
#2129 Answered sexual history questions has unmatched levels: -1
#2207 Wears glasses or contact lenses has unmatched levels: -1
#2296 Falls in the last year has unmatched levels: -1
#5090 Refractometry result unreliable (left) has unmatched levels: 2 E6 E8 E9 e4
#5091 Refractometry result unreliable (right) has unmatched levels: 2 E6 E8 E9 e4
#20003 Treatment/medication code has 3735 levels
#40001 Underlying (primary) cause of death: ICD10 has 1018 levels
#40002 Contributory (secondary) causes of death: ICD10 has 1291 levels
#41200 Operative procedures - main OPCS4 has 6394 levels
#41201 External causes - ICD10 has 1493 levels
#41202 Diagnoses - main ICD10 has 8360 levels
#41203 Diagnoses - main ICD9 has 2447 levels
#41204 Diagnoses - secondary ICD10 has 10161 levels
#41205 Diagnoses - secondary ICD9 has 2275 levels
#41210 Operative procedures - secondary OPCS4 has 7070 levels
#41270 Diagnoses - ICD10 has 11726 levels
#41271 Diagnoses - ICD9 has 3337 levels
#41272 Operative procedures - OPCS4 has 8124 levels

# Check
table(ukbdesc$ValueType[bd2ukbdesc], unlist(lapply(bd,function(x) class(x)[1])), useNA="a")
#					   Date POSIXct character factor integer numeric <NA>
#  Categorical multiple    0       0         0   3111       0       0    0
#  Categorical single      0       0         0   4854       0       0    0
#  Compound                0       0        32      0       0       0    0
#  Continuous              0       0         0      0       0    4092    0
#  Date                  693       0         0      0       0       0    0
#  Integer                 0       0         0      0    4072       0    0
#  Text                    0       0      1372      0       0       0    0
#  Time                    0     206         0      0       0       0    0
#  <NA>                    0       0         0      0       0       0    0

##########
# Output #
##########
# Save the main data frame
system.time((saveRDS(bd, "/fullpath/ukb41482.rds")))
#    user   system  elapsed
#1174.798   24.598 1201.095

# Save the helper data frames
saveRDS(ukbdesc, "/fullpath/ukb41482.ukbdesc.rds")
saveRDS(ukbcode, "/fullpath/ukb41482.ukbcode.rds")
