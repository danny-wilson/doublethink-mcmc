# Based on covars-self-report.R
# Add HES records too

# Run interactively
# screen -x 1363
# module add R/4.2.1-foss-2022a
# R

##############
# Parameters #
##############

# datestamp for output files
datestamp = "2023-09-18"

# Maximum date for HES codes as risk factors
bd.hes.maxdate = as.Date("2018-12-31")

# Maximum proportion of NAs in any column
missingness.threshold = 0.15

# Withdrawn individuals to remove
withdrawn.participants.filename = "/fullpath/w53100_2023-04-25.csv"

# Participant eids: English and not lost to follow up by 8 April 2020
#EnglishNotLostToFollowUp.filename = "/fullpath/ukb41482.English-not-lost-to-followup-8-April-2020.txt"

# Field classes to remove
remove.fieldclass = c("Compound", "Date", "Text", "Time")

# Maximum number of levels for factors, i.e. Categorical (single or multiple)
factors.max.nlevels = 50

# Override factors.max.nlevels for specific fields
override.factors.max.nlevels = unique(c(
20001, # Cancer code, self-reported, nlevels = 89
20002  # Non-cancer illness code, self-reported, nlevels = 474
#22601, # Job coding, nlevels = 995 [SPECIAL HANDLING NEEDED]
#41270, # Diagnoses - ICD10, nlevels = 19156 [FILTER NEEDED ON DATES]
))

# Minimum frequency for inclusion of a factor level
factors.min.freq = 0.002

# Remove negative factor levels EXCEPT when all factor levels are negative
remove.negative.levels = TRUE

# Specific fields to remove: genotyping and process-related fields
remove.field = unique(sort(paste0("f.", c(
3, 4, 5, 6, 21, 23, 35, 630, 5198, 6022, 10241, 20005, 21611, 21621, 21622, 21623, 21625, 21631, 21632, 21633, 21634, 21636, 21638, 21641, 21642, 21651, 21661, 21662, 21663, 21664, 21665, 21666, 21671, 21711, 21721, 21722, 21723, 21725, 21731, 21732, 21733, 21734, 21736, 21738, 21741, 21742, 21751, 21761, 21762, 21763, 21764, 21765, 21766, 21771, 21811, 21821, 21822, 21823, 21825, 21831, 21832, 21833, 21834, 21836, 21838, 21841, 21842, 21851, 21861, 21862, 21863, 21864, 21865, 21866, 21871, 22000, 22001, 22002, 22003, 22004, 22005, 22006, 22007, 22008, 22009, 22010, 22011, 22012, 22013, 22018, 22019, 22020, 22021, 22022, 22023, 22024, 22025, 22026, 22027, 22028, 22029, 22030, 22050, 22051, 22052, 22182, 23080, 23081, 23082, 23083, 23084, 23085, 23086, 23087, 23088, 23089, 23090, 23091, 23092, 23093, 23094, 23649, 23650, 23651, 23652, 23653, 23654, 23655, 23658, 23659, 23660, 24030, 24031, 24032, 24033, 24034, 24035, 24036, 24037, 24038, 24039, 24040, 24041, 24042, 24043, 24044, 24045, 24046, 24047, 24048, 24049, 24050, 24051, 24052, 24053, 24054, 24055, 24056, 24057, 24058, 24059, 24060, 24061, 24062, 24063, 24064, 24065, 24066, 24067, 24068, 24069, 28001, 28003, 28004, 28005, 28006, 28009, 28141, 28142, 28143, 30001, 30002, 30003, 30004, 30011, 30012, 30013, 30014, 30021, 30022, 30023, 30024, 30031, 30032, 30033, 30034, 30041, 30042, 30043, 30044, 30051, 30052, 30053, 30054, 30061, 30062, 30063, 30064, 30071, 30072, 30073, 30074, 30081, 30082, 30083, 30084, 30091, 30092, 30093, 30094, 30101, 30102, 30103, 30104, 30111, 30112, 30113, 30114, 30121, 30122, 30123, 30124, 30131, 30132, 30133, 30134, 30141, 30142, 30143, 30144, 30151, 30152, 30153, 30154, 30161, 30162, 30163, 30164, 30171, 30172, 30173, 30174, 30181, 30182, 30183, 30184, 30191, 30192, 30193, 30194, 30201, 30202, 30203, 30204, 30211, 30212, 30213, 30214, 30221, 30222, 30223, 30224, 30231, 30232, 30233, 30234, 30241, 30242, 30243, 30244, 30251, 30252, 30253, 30254, 30261, 30262, 30263, 30264, 30271, 30272, 30273, 30274, 30281, 30282, 30283, 30284, 30291, 30292, 30293, 30294, 30301, 30302, 30303, 30304, 30314, 30324, 30334, 30344, 30354, 30364, 30374, 30384, 30394, 30404, 30414, 30424, 30502, 30503, 30512, 30513, 30522, 30523, 30532, 30533, 30601, 30602, 30603, 30604, 30605, 30606, 30611, 30612, 30613, 30614, 30615, 30616, 30621, 30622, 30623, 30624, 30625, 30626, 30631, 30632, 30633, 30634, 30635, 30636, 30641, 30642, 30643, 30644, 30645, 30646, 30651, 30652, 30653, 30654, 30655, 30656, 30661, 30662, 30663, 30664, 30665, 30666, 30671, 30672, 30673, 30674, 30675, 30676, 30681, 30682, 30683, 30684, 30685, 30686, 30691, 30692, 30693, 30694, 30695, 30696, 30701, 30702, 30703, 30704, 30705, 30706, 30711, 30712, 30713, 30714, 30715, 30716, 30721, 30722, 30723, 30724, 30725, 30726, 30731, 30732, 30733, 30734, 30735, 30736, 30741, 30742, 30743, 30744, 30745, 30746, 30751, 30753, 30754, 30755, 30756, 30761, 30762, 30763, 30764, 30765, 30766, 30771, 30772, 30773, 30774, 30775, 30776, 30781, 30782, 30783, 30784, 30785, 30786, 30791, 30792, 30793, 30794, 30795, 30796, 30801, 30802, 30803, 30804, 30805, 30806, 30811, 30812, 30813, 30814, 30815, 30816, 30821, 30822, 30823, 30824, 30825, 30826, 30831, 30832, 30833, 30834, 30835, 30836, 30841, 30842, 30843, 30844, 30845, 30846, 30851, 30852, 30853, 30854, 30855, 30856, 30861, 30862, 30863, 30864, 30865, 30866, 30871, 30872, 30873, 30874, 30875, 30876, 30881, 30882, 30883, 30884, 30885, 30886, 30891, 30892, 30893, 30894, 30895, 30896, 30897, 40022, 40425, 41206, 41207, 41208, 41209, 41211, 41212, 41213, 41229, 41230, 41231, 41232, 41233, 41235, 41244, 41245, 41246, 41247, 41248, 41249, 41250, 41251, 41253, 42000, 42001, 42002, 42003, 42004, 42005, 42006, 42007, 42008, 42009, 42010, 42011, 42012, 42013, 42014, 42015, 42016, 42017, 42018, 42019, 42020, 42021, 42022, 42023, 42024, 42025, 42026, 42027, 42028, 42029, 42030, 42031, 42032, 42033, 42034, 42035, 42036, 42037
))))

# From thin-fields.ukb.R: identify English participants
f.assesscentre = "f.54.0.0"
assess.centre.England = c(
"11012 : Barts",
"11021 : Birmingham",
"11011 : Bristol",
"11008 : Bury",
#11003 : Cardiff
"11024 : Cheadle (revisit)",
"11020 : Croydon",
#11005 : Edinburgh
#11004 : Glasgow
"11018 : Hounslow",
"11010 : Leeds",
"11016 : Liverpool",
"11001 : Manchester",
"11017 : Middlesborough",
"11009 : Newcastle",
"11013 : Nottingham",
"11002 : Oxford",
"11007 : Reading",
"11014 : Sheffield",
"10003 : Stockport (pilot)",
"11006 : Stoke",
#11022 : Swansea
#11023 : Wrexham
"11025 : Cheadle (imaging)",
"11026 : Reading (imaging)",
"11027 : Newcastle (imaging)",
"11028 : Bristol (imaging)")

# Output file names
doublethink.input_filename = paste0("/fullpath/data.ukb41482.bd.selfreport.hes.", datestamp, ".csv")
doublethink.columns_filename = paste0("/fullpath/columns.ukb41482.bd.selfreport.hes.", datestamp, ".csv")

#############
# Functions #
#############
anymatch = function(x, y, ...) !is.na(match(x,c(y,...)))
remove.NA = function(x) ifelse(is.na(x),"",x)
remove.empty.brackets = function(s) ifelse(s=="()","",s)

# Expand the *first instance* of a Categorical (single) variable
# Does not create a separate column for NA values
expand.categorical.single = function(fieldID, bd, row.gd, factors.min.freq, check.levs=TRUE, remove.negative.levels=TRUE) {
	stopifnot(length(fieldID)==1)
	wh = grep(paste0("^f.",fieldID,".0.0$"), colnames(bd))
	stopifnot(length(wh)==1)
	levs = setdiff(unique(as.character(bd[row.gd, wh])), NA)
	if(check.levs) {
		fieldCoding = ukbdesc$Coding[ukbdesc$FieldID==fieldID]
		stopifnot(length(fieldCoding)==1)
		fieldLevs = paste(ukbcode$Value[ukbcode$Coding==fieldCoding], ":", ukbcode$Meaning[ukbcode$Coding==fieldCoding])
		stopifnot(all(anymatch(levs, fieldLevs)))
		if(!all(anymatch(levs, fieldLevs))) {
			stop(paste(fieldID,"has unmatched levels:",paste(setdiff(levs,fieldLevs),collapse=" ")))
		}
	}
	if(remove.negative.levels) {
		levs.neg = sapply(levs, function(s) {
			ilev = suppressWarnings(as.integer(unlist(strsplit(s, ":"))[1]))
			!is.na(ilev) & (ilev < 0)
		})
		if(!all(levs.neg)) levs = levs[!levs.neg]
	}
	bd.fieldID = bd[row.gd,0]
	for(LEV in levs) {
		lab = paste0("f.",fieldID,".0.0 : ", LEV)
		x = 1*(!is.na(bd[row.gd,wh]) & bd[row.gd,wh]==LEV)
		if(mean(x) >= factors.min.freq) bd.fieldID[,lab] <- x
	}
	return(bd.fieldID)
}

# Expand the *first instance* of a Categorical (multiple) variable
# Does not create a separate column for NA values
expand.categorical.multiple = function(fieldID, bd, row.gd, factors.min.freq, check.levs=TRUE, remove.negative.levels=TRUE) {
	stopifnot(length(fieldID)==1)
	wh = grep(paste0("^f.",fieldID,".0."), colnames(bd))
	levs = unique(as.character(bd[row.gd, wh[1]]))
	for(WH in wh[-1]) levs = union(levs, as.character(bd[row.gd, WH]))
	levs = setdiff(sort(levs), NA)
	if(check.levs) {
		fieldCoding = ukbdesc$Coding[ukbdesc$FieldID==fieldID]
		stopifnot(length(fieldCoding)==1)
		fieldLevs = paste(ukbcode$Value[ukbcode$Coding==fieldCoding], ":", ukbcode$Meaning[ukbcode$Coding==fieldCoding])
		stopifnot(all(anymatch(levs, fieldLevs)))
		if(!all(anymatch(levs, fieldLevs))) {
			stop(paste(fieldID,"has unmatched levels:",paste(setdiff(levs,fieldLevs),collapse=" ")))
		}
	}
	if(remove.negative.levels) {
		levs.neg = sapply(levs, function(s) {
			ilev = suppressWarnings(as.integer(unlist(strsplit(s, ":"))[1]))
			!is.na(ilev) & (ilev < 0)
		})
		if(!all(levs.neg)) levs = levs[!levs.neg]
	}
	bd.fieldID = bd[row.gd,0]
	for(LEV in levs) {
		lab = paste0("f.",fieldID,".0.* : ", LEV)
		x = 1*(rowSums(bd[row.gd,wh]==LEV, na.rm=TRUE)>0)
		if(mean(x) >= factors.min.freq) bd.fieldID[,lab] <- x
	}
	return(bd.fieldID)
}

# Manually review fields
more = function(txt, advance=45, width=100) {
	stopifnot(advance>0)
	nline = 0
	while(TRUE) {
		if(nline>length(txt)) break
		cat(paste0(substr(txt[pmin(nline+1:advance, length(txt))], 1, width),"\n"))
		nline = nline+advance
		readline()
	}
}

########
# Data #
########

# Read the UK Biobank main data for project 53100
system.time((bd = readRDS("/fullpath/ukb41482.rds")))
#   user  system elapsed
#291.128  24.890 316.031
dim(bd)
#[1] 502505  18432

# Read the Data Dictionary
ukbdesc = readRDS("/fullpath/ukb41482.ukbdesc.rds")

# Read the Data Codings
ukbcode = readRDS("/fullpath/ukb41482.ukbcode.rds")

bd.Rclass = unlist(lapply(bd,function(BD) class(BD)[1]))
table(bd.Rclass)
#	 Date   POSIXct character    factor   integer   numeric
#	  693       206      1404      7965      4072      4092

# Number of non-missing values per field
bd.ct = colSums(!is.na(bd))
summary(bd.ct)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   0     118    3272   44246   37439  502505

bd.field = gsub("^f\\.", "", gsub("\\.[0-9]+\\.[0-9]+$", "", colnames(bd)))
bd.class = ukbdesc$ValueType[match(bd.field, ukbdesc$FieldID)]
table(bd.class, useNA="a")
#Categorical multiple   Categorical single             Compound
#				3111                 4854                   32
#		  Continuous                 Date              Integer
#				4092                  693                 4072
#				Text                 Time                 <NA>
#				1372                  206                    0

# Number of levels per field
bd.nlev = unlist(lapply(bd,nlevels))

# Garbage collect
gc()

############################
# Remove repeated measures #
############################

# Take the first instance/array element of every field
# (This will be partially undone for Categorical (multiple) variables
col.gd = anymatch(1:ncol(bd), c(match("f.eid",colnames(bd)), grep(".0.0$",colnames(bd))))
sum(col.gd)
#[1] 2392

###############################
# Apply missingness threshold #
###############################

# But exclude 'Categorical multiple' from consideration: valid to have no entries, in
# which case the first column may have many NAs
col.gd = col.gd & ((bd.ct/nrow(bd)) >= (1-missingness.threshold) | bd.class=="Categorical multiple")
sum(col.gd)
#[1] 795

#################################
# Remove specific field classes #
#################################

col.gd = col.gd & !anymatch(bd.class, remove.fieldclass)
sum(col.gd)
#[1] 619

table(bd.class[col.gd])
# Categorical multiple   Categorical single           Continuous			 Integer
#                  114                  266                  144				  95

###########################
# Remove specific columns #
###########################

col.gd = col.gd & !anymatch(colnames(bd), paste0(remove.field,".0.0"))
sum(col.gd)
#[1] 392

#######################################
# Remove factors with too many levels #
#######################################

# List fields removed by this filter (unless overridden)
bd.field.nlev50 = unique(bd.field[bd.nlev>50])
cat(paste(bd.field.nlev50, ":", ukbdesc$Field[match(bd.field.nlev50, ukbdesc$FieldID)]), sep="\n# ")
# 20001 : Cancer code, self-reported								***
# 20002 : Non-cancer illness code, self-reported													***
# 20003 : Treatment/medication code
# 20004 : Operation code
# 20115 : Country of Birth (non-UK origin)
# 20199 : Antibiotic codes for last 3 months
# 21711 : Reception authorisation
# 21721 : Consent authorisation
# 21722 : Touchscreen authorisation
# 21723 : Touchscreen cognitive authorisation (pilot)
# 21725 : Touchscreen cognitive authorisation
# 21731 : Verbal interview authorisation
# 21732 : Pilot Spirometry authorisation
# 21733 : Measurement/Impedance authorisation (pilot)
# 21734 : Biometrics authorisation
# 21736 : Eye measures authorisation
# 21738 : ECG during exercise authorisation
# 21741 : Urine collection authorisation
# 21742 : Sample collection authorisation
# 21751 : Conclusion authorisation
# 21761 : Imaging screening authorisation
# 21762 : Brain MRI authorisation
# 21763 : Chest MRI authorisation
# 21764 : DXA assessment authorisation
# 21765 : Carotid ultrasound authorisation
# 21766 : ECG at rest authorisation
# 21771 : Cardiac monitor authorisation
# 22000 : Genotype measurement batch
# 22601 : Job coding												***
# 22617 : Job code - historical
# 40001 : Underlying (primary) cause of death: ICD10
# 40002 : Contributory (secondary) causes of death: ICD10
# 40006 : Type of cancer: ICD10
# 40011 : Histology of cancer tumour
# 40013 : Type of cancer: ICD9
# 41200 : Operative procedures - main OPCS4
# 41201 : External causes - ICD10
# 41202 : Diagnoses - main ICD10
# 41203 : Diagnoses - main ICD9
# 41204 : Diagnoses - secondary ICD10
# 41205 : Diagnoses - secondary ICD9
# 41210 : Operative procedures - secondary OPCS4
# 41229 : PCT responsible for patient data
# 41230 : PCT where patients GP was registered
# 41245 : Main speciality of consultant (recoded)
# 41246 : Treatment speciality of consultant (recoded)
# 41248 : Destinations on discharge from hospital (recoded)
# 41256 : Operative procedures - main OPCS3
# 41258 : Operative procedures - secondary OPCS3
# 41270 : Diagnoses - ICD10											*** (dates 41280)
# 41271 : Diagnoses - ICD9
# 41272 : Operative procedures - OPCS4
# 41273 : Operative procedures - OPCS3

# Implement the filter
col.gd = col.gd & ((bd.nlev <= factors.max.nlevels) | anymatch(bd.field, override.factors.max.nlevels))
sum(col.gd)
#[1] 377

#################################
# Remove withdrawn participants #
#################################

withdrawn.participants = scan(withdrawn.participants.filename, what=character(0))
row.gd = !anymatch(bd$f.eid, withdrawn.participants)
sum(row.gd)
#[1] 502367

############################################
# Remove non-English resident participants #
############################################

EnglishResident = anymatch(bd[,f.assesscentre], assess.centre.England)
sum(EnglishResident)
#[1] 445855
row.gd = row.gd & EnglishResident
sum(row.gd)
#[1] 445725

# NB some individuals are included that are known to have died (up to 2018-02-14)
# Need to exclude these downstream, depending on the outcome variable:
table("recorded_dead"=!is.na(bd$f.40001.0.0), "included"=row.gd)
#			 included
#recorded_dead  FALSE   TRUE
#		FALSE  54090 427983
#		TRUE    2690  17742

###########################
# Build the output matrix #
#  ! HIGH MEMORY USAGE !  #
###########################

bd.out = bd[row.gd,0]
bd.out$eid = as.integer(bd$f.eid[row.gd])
col.gd[match("f.eid", colnames(bd))] = FALSE

# Summarize final numbers of each field class type
table(bd.class[col.gd])
# Categorical multiple   Categorical single           Continuous			 Integer
#                   85                  108                  135				  49

# Construct the output matrix, taking each class in turn (some types not implemented)
for(BD.CLASS in unique(bd.class)) {
	wh = which(bd.class==BD.CLASS & col.gd)
	if(length(wh)>0) {
		if(BD.CLASS=="Integer" | BD.CLASS=="Continuous") {
			bd.out = cbind(bd.out, bd[row.gd, wh])
		} else if(BD.CLASS=="Categorical single") {
			for(WH in wh) {
				bd.out = cbind(bd.out, expand.categorical.single(bd.field[WH], bd, row.gd, factors.min.freq, TRUE, remove.negative.levels))
			}
		} else if(BD.CLASS=="Categorical multiple") {
			for(WH in wh) {
				bd.out = cbind(bd.out, expand.categorical.multiple(bd.field[WH], bd, row.gd, factors.min.freq, TRUE, remove.negative.levels))
			}
		} else {
			stop(paste("Processing",BD.CLASS,"not yet implemented"))
		}
		cat("Done",BD.CLASS,"\n")
	}
}

dim(bd.out)
#[1] 445725   1050

######################################
# Crude imputation - not implemented #
######################################

# Frequency range of missing values
bd.out.ct = colSums(!is.na(bd.out))
range(bd.out.ct/nrow(bd.out))
#[1] 0.847545 1.0000000

# Retain missing values, for downstream analysis to handle
if(FALSE) {
	wh = which(bd.out.ct < nrow(bd.out))
	bd.out.mn = rep(NA, ncol(bd.out))
	bd.out.mn[wh] = colMeans(bd.out[,wh], na.rm=TRUE)
	for(j in wh) bd.out[is.na(bd.out[, j]), j] <- bd.out.mn[j]
}

##############################
# User-friendly column names #
##############################

# Extract the fieldID for each field in bd.out
bd.out.field = c("eid", sapply(colnames(bd.out)[-1], function(s) regmatches(s, regexpr("([0-9]+)", s))))

# Extract additional information in the column names
bd.out.desc = trimws(sapply(colnames(bd.out), function(s) {
	paste0(unlist(strsplit(s, ":"))[-1], collapse=":")
}))

# Construct user-friendly column names
bd.out.new.colnames = trimws(paste(
	bd.out.field,
	remove.NA(ukbdesc$Field[match(bd.out.field, ukbdesc$FieldID)]),
	remove.empty.brackets(paste0("(", remove.NA(ukbdesc$Units[match(bd.out.field, ukbdesc$FieldID)]), ")"))
))
bd.out.new.colnames[bd.out.desc!=""] = paste(bd.out.new.colnames, ":", bd.out.desc)[bd.out.desc!=""]

if(!exists("bd.out.old.colnames")) bd.out.old.colnames = colnames(bd.out)
colnames(bd.out) = bd.out.new.colnames

##########################
# Check included columns #
##########################

# Print column names
more(colnames(bd.out))

########################################
# Add HES records until bd.hes.maxdate #
########################################

# Field 41270 contains combined primary and secondary ICD-10 diagnosis codes
# Field 41280 contains the first date associated with each diagnosis
# Cannot use them as they are out of date:
# range(apply(bd[,bd.field==41280],2,range,na.rm=TRUE),na.rm=TRUE)
# "1992-03-31" "2017-03-31"

# Instead read data directly from HES up to 2021-03-31
system.time((bd.hes = readRDS("/fullpath/ukb41482.hes.maxdate.2021-03-31.rds")))
#  user  system elapsed
#94.103   1.565  95.683
stopifnot(all(as.character(rownames(bd.hes)[row.gd])==bd.out$eid))

# Apply user-friendly column names
icd10code = read.delim("/fullpath/ukb-icd10-data-coding-19.tsv", stringsAsFactors=FALSE)
bd.hes.desc = icd10code$meaning[match(colnames(bd.hes), icd10code$coding)]
stopifnot(!any(is.na(bd.hes.desc)))

# Number of non-missing values per field
#bd.hes.ct = colSums(!is.na(bd.hes))
#summary(bd.hes.ct)
for(wh in 1:ncol(bd.hes)) {
	lab = bd.hes.desc[wh]
	x = 1*(!is.na(bd.hes[row.gd,wh]) & (bd.hes[row.gd,wh] <= bd.hes.maxdate))
	if(mean(x) >= factors.min.freq) {
		bd.out[,lab] <- x
		cat("Done",lab,"\n")
	}
}

dim(bd.out)
#[1] 445725   1913
# 2023-09-18: after moving bd.hes.maxdate forward to "2018-12-31"
#[1] 445725   1840

###############
# Output data #
###############

# Write the data
rownames(bd.out) <- as.integer((1:nrow(bd.out))-1)
system.time((write.table(bd.out, doublethink.input_filename, sep=",", quote=TRUE, row.names=TRUE, col.names=NA)))
#   user  system elapsed
#625.950   2.356 630.386

# Write the column types (all double for simplicity)
doublethink.columns = data.frame("colname"=colnames(bd.out), "type"="double")
write.table(doublethink.columns, doublethink.columns_filename, sep=",", quote=TRUE, row.names=FALSE, col.names=FALSE)
