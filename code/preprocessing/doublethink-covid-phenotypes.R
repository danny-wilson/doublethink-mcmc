#####################################################################
# COVID-19 phenotypes for doublethink analyses                      #
# Based on COVID-19 HGI definitions in hgi-phenotypes-v2.R          #
# With critical care criteria adapted from Steven Lin's version     #
# And a new facility to set a right-censoring date                  #
#####################################################################

	GENOFILTER = FALSE					# Use genotype and ancestry-based filters

# HGI Definition (17 August 2020):
# https://docs.google.com/document/d/1okamrqYmJfa35ClLvCt_vEe4PkvrTwggHq7T3jbeyCI/edit
# See also HGI UKB working group:
# https://docs.google.com/document/d/12XREkcRy8y8wW3Gr7H6PKhsaBVVw0yWVbWcPS-S_QDM/edit#
#
# Do not stratify by ancestry, age and sex
#
# Other criteria:
#          Do not exclude individuals with first degree relatives
#          Assessment centre used to define English participants
#          Not lost to follow-up for any other reason
#############
# Functions #
#############
anymatch = function(x, y, ...) !is.na(match(x,c(y,...)))
resolve.symlink = function(x) {
	y = Sys.readlink(x)
	if(any(is.na(y))) {
		stop("Could not resolve symlink ",x)
	}
	ifelse(y=="",x,file.path(dirname(x),y))
}

################
# Data logging #
################
# Every object in lg will be saved
lg = list()

##########################
# Command-line arguments #
##########################
help = paste(
"Usage: Rscript doublethink-covid-phenotypes.R YYYY-MM-DD",
sep="\n")
# Argument
lg$args = commandArgs(trailingOnly = TRUE)
print("Arguments:")
print(lg$args)
if(length(lg$args)!=1) {
	cat(help,sep="\n")
	stop("\nIncorrect usage\n")
}

lg$MINDATE = as.Date("2020-01-01")		# Hard code the left censoring date
#lg$MAXDATE = as.Date("2020-12-31")		# Specify the right censoring date [NOT YET IMPLEMENTED]
lg$MAXDATE = as.Date(lg$args[1])
if(is.na(lg$MAXDATE) | lg$MAXDATE < as.Date("2019-01-01")) stop("Nonsensical censoring date")
# outfile prefix
lg$stem = paste0("DOUBLETHINK-",lg$MAXDATE)
if(is.na(lg$stem) | lg$stem=="") stop("Nonsensical stem")

################
# Code version #
################
# Record GitLab version
print("User:")
print((lg$username = Sys.info()["user"]))
print("Source directory:")
print((lg$srcdir = paste0("/fullpath/",lg$username,"/GitLab/covid19gwas/")))
print("Git repository version")
system(paste0("(cd ",lg$srcdir," && git show --oneline -s)"))

##########################
# Input and output files #
##########################
# covid results file (use the latest)
lg$covid19_result.filename = resolve.symlink ("/fullpath/covid19_result.latest.txt")
# Mortality data (use the latest)
lg$death.filename = resolve.symlink("/fullpath/death.latest.txt")
lg$death_cause.filename = resolve.symlink("/fullpath/death_cause.latest.txt")
# Hospital episode statistics (use the latest)
lg$hesin.filename = resolve.symlink("/fullpath/hesin.latest.txt.gz")
lg$hesin_diag.filename = resolve.symlink("/fullpath/hesin_diag.latest.txt.gz")
lg$hesin_critical.filename = resolve.symlink("/fullpath/hesin_critical.latest.txt.gz")
# ukb (bd) data fields file
# bd.filename = "/fullpath/ukb41482.r"
# In RData format
lg$bd.RData.filename = "/fullpath/ukb41482.ukb41376.fields.RData"
# Remove individuals with first degree relatives
lg$remrels.filename = "/fullpath/ukb41482.English-remove-first-degree-relatives.eids.txt"
# Inferred MSOAs per individual
lg$homeloc.filename = "/fullpath/ukb41482.homelocation.msoa.txt.gz"
# Genotype QC data
lg$bed.sample.qc.filename = "/fullpath/ukb_sqc_v2.txt"
# Pre-computed eids for the bed-format genotypes
lg$bed.eid.filename = "/fullpath/analysis.bed.eids.txt"
# Ancestries as assigned by pan-UKB
lg$pan_ukb.filename = "/fullpath/all_pops_non_eur_pruned_within_pop_pc_covs.tsv"
lg$pan_ukb_bridge.filename = "/fullpath/ukb53100bridge31063.txt"
# individuals who issued data withdraw
lg$withdraw.filename = resolve.symlink("/fullpath/w53100_latest.csv")

# Output files
lg$log.outfilename = paste0("/fullpath/log.ukb41482.bd.gwasdata.",lg$stem,".txt")
lg$log.rds.outfilename = paste0("/fullpath/log.ukb41482.bd.gwasdata.",lg$stem,".rds")
# phenotype: UK Biobank data field reference
lg$data.outfilename = paste0("/fullpath/ukb41482.bd.gwasdata.",lg$stem,".txt")

# Enclose what follows in a tryCatch statement so the log is output
# even if R quits with an error
tryCatch( {
#############################
# Load epidemiological data #
#############################

# Load core epidemiological data to object 'bd'
system.time((load(lg$bd.RData.filename)))
#   user  system elapsed
# 17.031   2.182  19.272
all.eids = bd[,"f.eid"]

# Home locations and MSOA
homeloc = read.delim(lg$homeloc.filename,h=T,stringsAsFactors=F)

# Sample QC
bed.sample.qc = read.delim(lg$bed.sample.qc.filename,stringsAsFactors=FALSE,sep=' ')
# The corresponding eids
bed.eid = scan(lg$bed.eid.filename)
# Convert to bd.eid order
sample.qc = bed.sample.qc[match(all.eids,bed.eid),]

# Close (first degree) relatives
remrels = scan(lg$remrels.filename)

# Withdraw individuals
withdraw.eids = scan(lg$withdraw.filename)

# Ancestry as identified by pan-UKB
panukb.ancestry = read.delim(lg$pan_ukb.filename, sep = '\t', header = T, stringsAsFactors=F)[,c('s', 'pop')]
panukb.bridge = read.table(lg$pan_ukb_bridge.filename, sep = ' ', stringsAsFactors = F, header = F)
bridge.matched = panukb.bridge[match(panukb.ancestry$s, panukb.bridge[,2]),]
panukb.ancestry$eid = bridge.matched[,1]
panukb.ancestry.matched = panukb.ancestry[match(all.eids, panukb.ancestry$eid),]

##############################
# COVID-19 results wrangling #
##############################
covid19_result = read.delim(lg$covid19_result.filename,stringsAsFactors=FALSE)
# Fix date format
covid19_result$specdate = as.Date(covid19_result$specdate, tryFormats=c("%d/%m/%Y","%Y-%m-%d"))
# Filter by lg$MINDATE
if(!is.null(lg$MINDATE)) {
	covid19_result.keep = covid19_result$specdate >= lg$MINDATE
	covid19_result.keep[is.na(covid19_result.keep)] = FALSE
	covid19_result = covid19_result[covid19_result.keep,]
}
# Filter by lg$MAXDATE
if(!is.null(lg$MAXDATE)) {
	covid19_result.keep = covid19_result$specdate <= lg$MAXDATE
	covid19_result.keep[is.na(covid19_result.keep)] = FALSE
	covid19_result = covid19_result[covid19_result.keep,]
}
# All eids with a test result
covid19_result.eids = unique(covid19_result$eid)
# For each, are there any positive tests
covid19_result.eids.anypos = sapply(covid19_result.eids,function(EID) any(covid19_result$result[!is.na(match(covid19_result$eid,EID))]==1))
# For each, are there any positive tests *while* they were an inpatient
covid19_result.eids.anypos.inpatient = sapply(covid19_result.eids,function(EID) any(covid19_result$result[!is.na(match(covid19_result$eid,EID))]==1 & covid19_result$origin[!is.na(match(covid19_result$eid,EID))]==1))
# Were they *ever* flagged an inpatient when tested, irrespective of result?
covid19_result.eids.inpatient = sapply(covid19_result.eids,function(EID) any(covid19_result$origin[!is.na(match(covid19_result$eid,EID))]==1))
# Compare
print("UKB covid19_result table: any positive vs any positive while inpatient")
print((lg$table.covid19_result.eids.anypos.vs.anypos.inpatient =
table(covid19_result.eids.anypos,covid19_result.eids.anypos.inpatient)))
print("UKB covid19_result table: any positive vs any test while inpatient")
print((lg$table.covid19_result.eids.anypos.vs.inpatient = table(covid19_result.eids.anypos,covid19_result.eids.inpatient)))

# For each eid in the bd file, create anypos and inpatient phenotypes
bd.covid19_result.anypos = rep(NA,nrow(bd))
bd.covid19_result.anypos[match(covid19_result.eids[covid19_result.eids.anypos],all.eids)] = 1
bd.covid19_result.anypos[match(covid19_result.eids[!covid19_result.eids.anypos],all.eids)] = 0
bd.covid19_result.inpat = rep(NA,nrow(bd))
bd.covid19_result.inpat[match(covid19_result.eids[covid19_result.eids.inpatient],all.eids)] = 1
bd.covid19_result.inpat[match(covid19_result.eids[!covid19_result.eids.inpatient],all.eids)] = 0
print("UKB covid19_result table: any positive vs inpatient, sanity check")
print((lg$table.bd.covid19_result.anypos.vs.inpat = table(bd.covid19_result.anypos,bd.covid19_result.inpat,useNA="a")))
stopifnot(all(lg$table.bd.covid19_result.anypos.vs.inpat[1:2,1:2]==lg$table.covid19_result.eids.anypos.vs.inpatient))

# Define four groups and baseline (no.test)
# NB: anypos.in requires only that there was *ever* a positive test and
# they were *ever* an inpatint while tested, not necessarily at the same time
bdc.group = ifelse(bd.covid19_result.anypos,ifelse(bd.covid19_result.inpat,"anypos.in","anypos.nin"),ifelse(bd.covid19_result.inpat,"neg.in","neg.nin"))
bdc.group[is.na(bdc.group)] = "no.test"
bdc.group = factor(bdc.group,levels=c("anypos.in","anypos.nin","neg.in","neg.nin","no.test"))
print("UKB covid19_result table: group assignments")
print((lg$table.bdc.group = table(bdc.group,useNA="a")))

##################
# Mortality data #
##################
death = read.delim(lg$death.filename,as.is=T,stringsAsFactors=F)
death$date_of_death = as.Date(death$date_of_death,tryFormats=c("%d/%m/%Y","%Y-%m-%d"))
# Record eids of any participant that died before mindate
anymort.leftcensor = c()
# Filter by lg$MINDATE
if(!is.null(lg$MINDATE)) {
	anymort.leftcensor = sort(unique(death$eid[death$date_of_death<lg$MINDATE]))
	death.keep = death$date_of_death >= lg$MINDATE
	death.keep[is.na(death.keep)] = FALSE
	death = death[death.keep,]
}
# Filter by lg$MAXDATE
if(!is.null(lg$MAXDATE)) {
	death.keep = death$date_of_death <= lg$MAXDATE
	death.keep[is.na(death.keep)] = FALSE
	death = death[death.keep,]
}
death_cause = read.delim(lg$death_cause.filename,as.is=T,stringsAsFactors=F)
# Only include primary cases to follow definition of phenotype A1 (below)
print("COVID-19 causes of death in UKB death_cause table")
print((lg$table.death_cause.U07 = table(factor(death_cause$cause_icd10,levels=c("U071","U072")))))
# Unique EIDs
mort.primary.eids = unique(death_cause$eid[death_cause$cause_icd10=="U071"])
mort.secondary.eids = unique(death_cause$eid[death_cause$cause_icd10=="U072"])

# Any participant who died in the censored study period
anymort.eids = sort(unique(death$eid))

###############################
# Hospital episode statistics #
###############################
hesin = read.delim(lg$hesin.filename,as.is=TRUE)
hesin_diag = read.delim(lg$hesin_diag.filename,as.is=TRUE)
hesin_critical = read.delim(lg$hesin_critical.filename, as.is=TRUE)
# Format dates correctly
hesin$epistart = as.Date(hesin$epistart,tryFormats=c("%d/%m/%Y","%Y-%m-%d"))
hesin$epiend = as.Date(hesin$epiend,tryFormats=c("%d/%m/%Y","%Y-%m-%d"))
hesin$elecdate = as.Date(hesin$elecdate,tryFormats=c("%d/%m/%Y","%Y-%m-%d"))
hesin$admidate = as.Date(hesin$admidate,tryFormats=c("%d/%m/%Y","%Y-%m-%d"))
hesin$disdate = as.Date(hesin$disdate,tryFormats=c("%d/%m/%Y","%Y-%m-%d"))
hesin_critical$ccstartdate = as.Date(hesin_critical$ccstartdate,tryFormats=c("%d/%m/%Y","%Y-%m-%d"))
hesin_critical$ccdisdate = as.Date(hesin_critical$ccdisdate,tryFormats=c("%d/%m/%Y","%Y-%m-%d"))

# Filter by lg$MINDATE
if(!is.null(lg$MINDATE)) {
	hesin.keep = !is.na(hesin$epistart) & (hesin$epistart >= lg$MINDATE)
	hesin_critical.keep = !is.na(hesin_critical$ccstartdate) & hesin_critical$ccstartdate >= lg$MINDATE
	hesin = hesin[hesin.keep,]
	hesin_critical = hesin_critical[hesin_critical.keep,]
}
# Filter by lg$MAXDATE
if(!is.null(lg$MAXDATE)) {
	hesin.keep = !is.na(hesin$epiend) & (hesin$epiend <= lg$MAXDATE)
	hesin_critical.keep = !is.na(hesin_critical$ccdisdate) & hesin_critical$ccdisdate <= lg$MAXDATE
	hesin = hesin[hesin.keep,]
	hesin_critical = hesin_critical[hesin_critical.keep,]
}

# Merge, filtering on U07? diagnosis code
hesin_merge = merge(hesin[!is.na(hesin$admidate),]
	,hesin_diag[anymatch(hesin_diag$diag_icd10,c("U071","U072")),],by=c("eid","ins_index"),all=TRUE,stringsAsFactors=F)
hesin_merge = merge(hesin_merge, hesin_critical, by=c("eid", "ins_index"), all=TRUE, stringsAsFactors=F)

hesin.u071.eids = setdiff(unique(hesin_merge$eid[
	hesin_merge$classpat_uni==1000 &
	!is.na(hesin_merge$epiend) &
	hesin_merge$diag_icd10=="U071"
]), NA)
hesin.u07.eids = setdiff(unique(hesin_merge$eid[
	hesin_merge$classpat_uni==1000 &
	!is.na(hesin_merge$epiend) &
	(hesin_merge$diag_icd10=="U071" | hesin_merge$diag_icd10=="U072")
]), NA)

# Identify patients with U07 code and receive respiratory support within 30 days
hesin.ressprt.eids = setdiff(unique(hesin_merge$eid[
	hesin_merge$classpat_uni==1000 &
	!is.na(hesin_merge$ccstartdate) &
	( ( hesin_merge$bressupdays > 0 ) | # Basic respiratory support
	  ( hesin_merge$aressupdays > 0 ) ) # Advanced respiratory support
]),NA)
hesin.u07.ressprt.eids.check = intersect(hesin.ressprt.eids, hesin.u07.eids)
hesin.u07.ressprt.eids = hesin.u07.ressprt.eids.check[sapply(hesin.u07.ressprt.eids.check, function(EID) {
	ressprt.episodes = setdiff(which(hesin_merge$eid==EID &
		hesin_merge$classpat_uni==1000 &
		!is.na(hesin_merge$ccstartdate) &
		( ( hesin_merge$bressupdays > 0 ) | ( hesin_merge$aressupdays > 0 ) )), NA)
	ressprt.dates = as.vector(unlist(sapply(ressprt.episodes, function(wh) {
		hesin_merge$ccstartdate[wh]:hesin_merge$ccdisdate[wh]
	})))
	diag.episodes = setdiff(which(hesin_merge$eid==EID &
		hesin_merge$classpat_uni==1000 &
		!is.na(hesin_merge$epiend) &
		(hesin_merge$diag_icd10=="U071" | hesin_merge$diag_icd10=="U072")), NA)
	diag.dates = as.vector(unlist(sapply(diag.episodes, function(wh) {
		hesin_merge$epistart[wh]:hesin_merge$epiend[wh]
	})))
	min(outer(ressprt.dates, diag.dates, function(x,y) abs(x-y))) <= 30
})]

print("Number of participants with diagnosis code U071")
print((lg$nparticipant.u071 = length(hesin.u071.eids)))
print("Number of participants with diagnosis code U071 or U072")
print((lg$nparticipant.u07 = length(hesin.u07.eids)))
print("Number of participants with respiratory support")
print((lg$nparticipant.ressport = length(hesin.ressprt.eids)))
print("Number of participants with respiratory support and U071 or U072")
print((lg$nparticipant.u07.ressport = length(hesin.u07.ressprt.eids)))

# Participants with evidence of hospitalization since 1st January 2020
# query.inpat.2020.eids = union(covid19_result.eids.inpatient, hesin_merge$eid)
query.inpat.2020.eids = union(covid19_result.eids[covid19_result.eids.inpatient], hesin_merge$eid)


###########################################
# Analysis A1 phenotype                   #
# VERY SEVERE RESPIRATORY CONFIRMED COVID #
###########################################
# Cases: Hospitalized laboratory confirmed SARS-CoV-2 infection AND death
#        AND hospitalization with COVID19 as ++ primary reason ++ for admission
# Ctrls: Laboratory confirmed SARS-CoV-2 infection
#        AND not hospitalised ** 21 days after the test **
# ++ NB: No information on primary vs secondary reasons, but use HES ICD10 code U071
# ** NB: Substitute this condition for anypos.nin or any HES record in 2020

# Not done: Add Z99.1 (Dependence on respirator) Z99.11 (Dependence on respirator, status)

# Counts
print("covid19_result group vs U071 death")
print((lg$table.bdc.group.vs.u071.died = table('test_status'=bdc.group,'u071_death'=anymatch(all.eids,mort.primary.eids))))

# What about folk who died (in study period),
# not registered as COVID-19-related, but tested positive?
bdc.mort = ifelse(anymatch(all.eids,mort.primary.eids),"diedU071",
           ifelse(anymatch(all.eids,mort.secondary.eids),"diedU072",
           ifelse(anymatch(all.eids,anymort.eids),"diedothr",
           "other")))
print("covid19_result group vs death ICD10 codes")
print((lg$table.bdc.group.vs.mort.icd10 = table(bdc.group,bdc.mort)))

# Unfiltered, unstratified phenotype definition
pheno.A1 =
#	Case if ever PCR positive AND (died or received respiratory support) with U071 or U072 code AND HES inpatient with diagnosis code U071 or U072
	ifelse(anymatch(all.eids,covid19_result.eids[covid19_result.eids.anypos]) &
	      (anymatch(all.eids,mort.primary.eids) |
		   anymatch(all.eids,mort.secondary.eids) |
		   anymatch(all.eids,hesin.u07.ressprt.eids)) &
		   anymatch(all.eids,hesin.u07.eids),1,
#	Ctrl if ever positive AND did not (die or receive respiratory support) with U071 nor U072 AND no evidence of hospitalization in study period
	ifelse(anymatch(all.eids,covid19_result.eids[covid19_result.eids.anypos]) &
		  !anymatch(all.eids,mort.primary.eids) &
		  !anymatch(all.eids,mort.secondary.eids) &
		  !anymatch(all.eids,hesin.u07.ressprt.eids) &
		  !anymatch(all.eids,query.inpat.2020.eids),0,
#	Excluded otherwise
		   -1))
# Sanity check: should be no NAs
stopifnot(!any(is.na(pheno.A1)))
print("Unfiltered, unstratified phenotype A1 counts:")
print((lg$table.raw.pheno.A1 = table(pheno.A1,useNA="a")))

lenient.pheno.A1 =
#	Case if ever PCR positive AND [(died or respiratory support) with U071 or U072 code OR ((died or respiratory support) AND HES inpatient with diagnosis code U071 or U072)]
	ifelse(anymatch(all.eids,covid19_result.eids[covid19_result.eids.anypos]) &
		  (anymatch(all.eids,mort.primary.eids) |
	       anymatch(all.eids,mort.secondary.eids) |
		   anymatch(all.eids,hesin.u07.ressprt.eids) |
		   (anymatch(all.eids,anymort.eids) &
		    anymatch(all.eids,hesin.u07.eids)) |
		   (anymatch(all.eids,hesin.ressprt.eids) &
			anymatch(all.eids,hesin.u07.eids))),1,
#	Ctrl if ever positive AND did not die with U071 nor U072 AND no evidence of hospitalization in 2020
	ifelse(anymatch(all.eids,covid19_result.eids[covid19_result.eids.anypos]) &
		  !anymatch(all.eids,mort.primary.eids) &
		  !anymatch(all.eids,mort.secondary.eids) &
		  !anymatch(all.eids,query.inpat.2020.eids),0,
#	Excluded otherwise
		   -1))
print("Unfiltered, unstratified phenotype A1 counts (lenient definition):")
print((lg$table.raw.lenient.pheno.A1 = table(lenient.pheno.A1,useNA="a")))

###########################################
# Analysis A2 phenotype                   #
# VERY SEVERE RESPIRATORY CONFIRMED COVID #
###########################################
# Cases: Same as A1
# Ctrls: Everyone that is not a case, i.e. the population

pheno.A2 = c("-1"=0,"0"=0,"1"=1)[as.character(pheno.A1)]
lenient.pheno.A2 = c("-1"=0,"0"=0,"1"=1)[as.character(lenient.pheno.A1)]
print("Unfiltered, unstratified phenotype A2 counts:")
print((lg$table.raw.pheno.A2 = table(pheno.A2,useNA="a")))
print("Unfiltered, unstratified phenotype A2 counts (lenient definition):")
print((lg$table.raw.lenient.pheno.A2 = table(lenient.pheno.A2,useNA="a")))

####################################
# Analysis B1 phenotype             #
# HOSPITALIZED LAB CONFIRMED COVID #
####################################
# Cases: Hospitalized laboratory confirmed SARS-CoV-2 infection AND
#        hospitalization due to corona-related symptoms.
# Ctrls: Laboratory confirmed SARS-CoV-2 infection AND
#        AND not hospitalised ** 21 days after the test **
# ** NB: Substitute this condition for anypos.nin or any HES record in 2020
# !! Also exclude from controls anyone that died with a COVID-19 code

# Unfiltered, unstratified phenotype definition
pheno.B1 =
#	Case if ever PCR positive AND HES inpatient with diagnosis code U071 or U072
	ifelse(anymatch(all.eids,covid19_result.eids[covid19_result.eids.anypos]) &
		   anymatch(all.eids,hesin.u07.eids),1,
#	Ctrl if ever positive AND did not die with U071 nor U072 AND no evidence of hospitalization in 2020
	ifelse(anymatch(all.eids,covid19_result.eids[covid19_result.eids.anypos]) &
		  !anymatch(all.eids,mort.primary.eids) &
		  !anymatch(all.eids,mort.secondary.eids) &
		  !anymatch(all.eids,query.inpat.2020.eids),0,
#	Excluded otherwise
		   -1))
# Sanity check: should be no NAs
stopifnot(!any(is.na(pheno.B1)))
print("Unfiltered, unstratified phenotype B1 counts:")
print((lg$table.raw.pheno.B1 = table(pheno.B1,useNA="a")))

####################################
# Analysis B2 phenotype            #
# HOSPITALIZED LAB CONFIRMED COVID #
####################################
# Cases: Same as B1
# Ctrls: Everyone that is not a case, i.e. the population

pheno.B2 = c("-1"=0,"0"=0,"1"=1)[as.character(pheno.B1)]
print("Unfiltered, unstratified phenotype B2 counts:")
print((lg$table.raw.pheno.B2 = table(pheno.B2,useNA="a")))


##########################
# Analysis C1 phenotype  #
# PARTIAL-SUSCEPTIBILITY #
##########################
# Cases: Laboratory confirmation of SARS-CoV-2 infection OR
#        EHR/ICD coding/ Physician Confirmed COVID-19 ** OR
#        self-reported COVID-19 positive (e.g. by questionnaire)
# Ctrls: (Laboratory tested for SARS-CoV-2 infection AND
#         all tests (if multiple tests) negative*) OR
#        self-reported tested negative for SARS-CoV-2 infection (e.g. by questionnaire)
# ** (See Appendix 1 for suggestive codes: U07.1, U07.2, Z99.1, Z99.11)
# !! Also exclude from controls anyone that died with a COVID-19 code or had a HES diagnosis code


# Unfiltered, unstratified phenotype definition
pheno.C1 =
#	Case if ever PCR positive OR died with U071 or U072 code OR HES inpatient with diagnosis code U071 or U072
	ifelse(anymatch(all.eids,covid19_result.eids[covid19_result.eids.anypos]) |
	       anymatch(all.eids,mort.primary.eids) |
		   anymatch(all.eids,mort.secondary.eids) |
		   anymatch(all.eids,hesin.u07.eids),1,
#	Ctrl if never positive AND did not die with U071 nor U072 AND no HES inpatient with diagnosis code U071 or U072
	ifelse(anymatch(all.eids,covid19_result.eids[!covid19_result.eids.anypos]) &
		  !anymatch(all.eids,mort.primary.eids) &
		  !anymatch(all.eids,mort.secondary.eids) &
		  !anymatch(all.eids,hesin.u07.eids),0,
#	Excluded otherwise
		   -1))
# Sanity check: should be no NAs
stopifnot(!any(is.na(pheno.C1)))
print("Unfiltered, unstratified phenotype C1 counts:")
print((lg$table.raw.pheno.C1 = table(pheno.C1,useNA="a")))

##########################
# Analysis C2 phenotype  #
# PARTIAL-SUSCEPTIBILITY #
##########################
# Cases: Same as C1
# Ctrls: Everyone that is not a case, i.e. the population

pheno.C2 = c("-1"=0,"0"=0,"1"=1)[as.character(pheno.C1)]
print("Unfiltered, unstratified phenotype C2 counts:")
print((lg$table.raw.pheno.C2 = table(pheno.C2,useNA="a")))

#######################
# Join the phenotypes #
#######################

pheno = data.frame("eid"=all.eids, "A1"=pheno.A1, "A2"=pheno.A2, "B1"=pheno.B1, "B2"=pheno.B2, "C1"=pheno.C1, "C2"=pheno.C2, "lenient.A1"=lenient.pheno.A1)

print("Pre any filters pre stratification: phenotype counts")
print((lg$table.raw.pheno = apply(pheno[,-1],2,function(y) {
	table(factor(y,levels=-1:1),useNA="a")
})))

###########
# Filters #
###########
# From thin-fields.ukb.R: identify English participants
f.assesscentre = "f.54.0.0"
assess.centre.England = c(
11012, #	Barts
11021, #	Birmingham
11011, #	Bristol
11008, #	Bury
#11003	Cardiff
11024, #	Cheadle (revisit)
11020, #	Croydon
#11005	Edinburgh
#11004	Glasgow
11018, #	Hounslow
11010, #	Leeds
11016, #	Liverpool
11001, #	Manchester
11017, #	Middlesborough
11009, #	Newcastle
11013, #	Nottingham
11002, #	Oxford
11007, #	Reading
11014, #	Sheffield
10003, #	Stockport (pilot)
11006, #	Stoke
#11022	Swansea
#11023	Wrexham
11025, #	Cheadle (imaging)
11026, #	Reading (imaging)
11027, #	Newcastle (imaging)
11028) #	Bristol (imaging)

# Sanity check: test results must be in England
sample.qc.English.participant = !is.na(match(bd[,f.assesscentre],assess.centre.England))

# Total number of 'non-English' participants with results, according to assessment centre
print("English recruitment centre vs phenotype:")
print((lg$table.English.recruitment.centre.vs.pheno = table(sample.qc.English.participant,pheno.A1,useNA="a")))

# Now use last known address to define English participants
homeloc.English.participant = !is.na(homeloc$ela[match(all.eids,homeloc$eid)])
print("English home location vs phenotype:")
print((lg$table.English.homeloc.vs.pheno = table(homeloc.English.participant,pheno.A1,useNA="a")))

# Stick with assessment centre for now to define English/non-English participants

# Identify those lost to follow-up (withdrawn or died before study period)
sample.qc.lost2followup.eids = sort(unique(c(withdraw.eids, anymort.leftcensor)))
sample.qc.lost2followup = anymatch(bd$f.eid, sample.qc.lost2followup.eids)
# Sanity check: should rule out non-English tests
print("Sanity check: non-English or lost to follow-up versus phenotype")
print((lg$table.sample.qc.lost2followup.vs.pheno = table(sample.qc.lost2followup,pheno.A1,useNA="a")))

# To be sure, wipe out these results
pheno[sample.qc.lost2followup | !sample.qc.English.participant,-1] = -1
# Reproduce phenotype counts after non-English resident
# and lost to follow-up filters
print("Residence and follow-up filters: number included")
print((lg$table.filtered.resfollowup.counts = table(!apply(pheno[,-1]==-1,1,all),useNA="a")))

if(GENOFILTER) {
	# Exclude folk on various other grounds
	# Apply quality control filters
	gd.qc.pre1 = sample.qc.lost2followup==FALSE &  # Repeated for completeness
		sample.qc.English.participant==TRUE &      # Repeated for completeness
		sample.qc$het.missing.outliers==0 &
		sample.qc$putative.sex.chromosome.aneuploidy==0 &
		sample.qc$Submitted.Gender==sample.qc$Inferred.Gender &
		sample.qc$excluded.from.kinship.inference==0 &
		sample.qc$excess.relatives==0 &
		sample.qc$in.Phasing.Input.chr1_22==1 &
		sample.qc$in.Phasing.Input.chrX==1 &
		sample.qc$in.Phasing.Input.chrXY==1
	gd.qc.pre1[is.na(gd.qc.pre1)] = FALSE
	print("Filtering counts post main filters:")
	print((lg$table.filter.pre1.counts = table(gd.qc.pre1,useNA="a")))
	# If one were to also exclude based on used.in.pca.calculation==1
	gd.qc.pre2.old = gd.qc.pre1 & sample.qc$used.in.pca.calculation==1
	print("Filtering counts post old PCA filter (not used):")
	print((lg$table.filter.pre2.old.counts = table(gd.qc.pre2.old,useNA="a")))
	# New: exclude only first degree relatives
	gd.qc.pre2 = gd.qc.pre1 & is.na(match(all.eids,remrels))
	print("Phenotype counts post close relatives filter:")
	print((lg$table.filter.pre2.new.counts = table(gd.qc.pre2,useNA="a")))
	# Post filter summary
	# Note that stratification (ancestry, age, sex) is done downstream
	gd.qc = gd.qc.pre2
	print("Filtering counts post filters:")
	print((lg$table.filter.gd.qc.counts = table(gd.qc,useNA="a")))
	# Apply the filters
	pheno[is.na(gd.qc) | !gd.qc,-1] = -1
}

print("Post all filters pre stratification: phenotype counts")
print((lg$table.filtered.pheno = apply(pheno[,-1],2,table,useNA="a")))


#########################
# Construct a data file #
#########################
f.sex = "f.31.0.0"
f.birthyear = "f.34.0.0"
ancestry = c("oth","eur")[1+sample.qc$in.white.British.ancestry]
# Create data matrix in the same order as the bd matrix
data = pheno
data$sex = as.character(bd[,f.sex])
data$age = (2020-bd[,f.birthyear])
data$ancestry = ancestry
data$ancestry_panukb = panukb.ancestry.matched$pop
# Count NAs per column
print("NA counts per column in data:")
print((lg$na.counts.columns.data = colSums(is.na(data))))

# There is one individual with NA sex and NA age
# Looks like this individual also has all -1 (exclude) codes:
print("Phenotypes for individuals NA for sex or age:")
print((lg$data.na.for.sex.or.age = data[is.na(data$age) | is.na(data$sex),]))

# Summary statistics
print("# cases")
print((lg$n.cases = colSums(pheno[,-1]==1)))
print("# controls")
print((lg$n.controls = colSums(pheno[,-1]==0)))
print("# cases+controls")
print((lg$n.cases.plus.controls = colSums(pheno[,-1]!=-1)))
print("% female")
print((lg$percent.female = round(1000*(mean(data$sex=="Female",na.rm=TRUE)))/10))
print("Mean age")
print((lg$mean.age = round(10*mean(data$age,na.rm=TRUE))/10))
print("SD age")
print((lg$sd.age = round(10*sd(data$age,na.rm=TRUE))/10))
print("Min date of PCR test")
print((lg$min.specdate = min(covid19_result$specdate)))
print("Max date of PCR test")
print((lg$max.specdate = max(covid19_result$specdate)))
print("Min date of hospital episode")
print((lg$min.hes = min(hesin_merge$epistart,na.rm=TRUE)))
print("Max date of hospital episode")
print((lg$max.hes = max(hesin_merge$epiend,na.rm=TRUE)))
print("Min date of critical care")
print((lg$min.hcc = min(hesin_merge$ccstartdate,na.rm=TRUE)))
print("Max date of critical care")
print((lg$max.hcc = max(hesin_merge$ccdisdate,na.rm=TRUE)))
print("Min date of death")
print((lg$min.dod = min(death$date_of_death)))
print("Max date of death")
print((lg$max.dod = max(death$date_of_death)))
print("Max admission date of U07 diagnosis")
print((lg$max.u07.admidate = max(hesin$admidate[match(hesin.u07.eids,hesin$eid)])))
print("Max discharge date of U07 diagnosis (excluding NAs)")
print((lg$max.u07.disdate = max(hesin$disdate[match(hesin.u07.eids,hesin$eid)],na.rm=TRUE)))

# Simple log odds ratios
pheno.COLs = c("A1","A2","B1","B2","C1","C2")
system.time((fit = lapply(pheno.COLs, function(COL) {
	pheno.COL = ifelse(pheno[,COL]==-1,NA,pheno[,COL])
	glm(pheno.COL ~ sex + ancestry + age, data=data, family="binomial")
})))
names(fit) = pheno.COLs
print("log odds ratios:")
print((lg$coef = lapply(fit,function(FIT) signif(summary(FIT)$coef,3))))

# Output the data construct
system.time((write.table(data,file=lg$data.outfilename,row=FALSE,col=TRUE,quote=FALSE,sep='\t')))

}, finally = {
	# On error or clean exit
	saveRDS(lg, file=lg$log.rds.outfilename)
})
