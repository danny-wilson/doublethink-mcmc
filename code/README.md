# Doublethink Readme 

> NB: this repository is published to accompany the paper. It is released with an example dataset (see README.md in the parent directory), which is intended to be functional and representative, but does not replicate the analysis in the paper for several reasons: (I) Expertise implementing the Markov Chain Monte Carlo algorithm in a high performance computing cluster environment is required. Where necessary we have redacted full file paths for security reasons. (II) To access UK Biobank data, researchers must register and submit a research application (https://www.ukbiobank.ac.uk/register-apply). Registration is open to all bona fide researchers for all types of health-related research that is in the public interest. The registration and application process ensures researchers and projects meet UK Biobank's obligations to its participants and funders. (III) Going forward, applications of Doublethink to UK Biobank will need to implement the code within the UK Biobank Research Analysis Platform, due to changes in UK Biobank data access policy. To install the code on another system, please base the software requirements on the Docker file in the parent directory.

> NB: The following README was originally written by N.A. for internal use during development. For the GitHub repo, see README.md in the parent directory.

Doublethink is a program that can calculate risk factor for any given disease based on Bayesian model averaging. The idea and initial code are by Daniel J. Wilson, the theory was developed by Helen R. Fryer and Daniel J. Wilson, while code development, maintenance and chain initialisation are by Nicolas Arning. In the following we'll first explain what the different files in this directory are, how to run the code and how to interpret the results, to enable further use of doublethink. Lastly, I will go through the data pre-processing we performed on UK Biobank data, so the pre-processing can be repeated on other datasets.

## Quick F.A.Q.

This is a very brief summary of the most urgent questions you might have. More details will follow later
**Why doublethink?**
From cliff notes

> In George Orwell's dystopian classic 1984, doublethink is the act of holding, simultaneously, two opposite, individually exclusive ideas or opinions and believing in both simultaneously and absolutely.

Our approach allows you to compute both Bayesian posterior probability and frequentist p-values, and therefore is interpretable from both perspectives. 

**How do I install doublethink?**

You will also need to install the necessary python libraries from `Python_packages.txt` and R libraries from `R_packages.txt`. Doublethink was developed using Python 3.8.3 and R 3.6.2. You will need to alias the python command to run the relevant python3 version.

**How do I run doublethink?**

```
./qsub_wrapper.sh ./configs/<config of your choice>.cfg 
```
There are a few example configs in the config folders. The variables are briefly explained in there

**Where do I run doublethink?**

Wherever you want really. Just change the workdir variable in the config. This is where all the log and output files will be.

**Where are the config files?**

in the configs subfolder

**What's in the config files?**

Check any of the config files, the parameters are explained in comments. Otherwise check in the readme to follow there is a more detailed explanation.

**Where are the pheno files?**

You will have to supply those yourself. It consists of an eid column and a phenotype column with 0s and 1s. You can change the name of the identifier column in the config using idcol. You can change the name of the phenotype column in the config using phenocol.

**Where are the log files?**

Within the workdir folder there are `<name variable>_chain_<chain id>.out` which is the standard out and `<name variable>_chain_<chain id>.err` which is the error log. Troubleshooting starts here. 

**Where are the output files?**

In whatever folder you defined the workdir to be in the config file. The most interesting and human readable will be the parsed table under `results.<name variable>.summary.tsv` and the pdf with the plots `results.<name variable>.plots.pdf`

With that out of the way let's get to a more detailed explanation of the program.

## Description of files
The main body of the bayesian model averaging (and all the theory that goes with its) is contained in the `doublethink-functions.R` which has all the code to run the MCMC chains. The `run-mcmc.R` preprocesses the data and starts the chains calling the `make_starting_groups.py` from inside the script to initialise the chains using the furthest neighbour algorithm. The chain files denoted by `.RDS` in the output directory are consequently processed with the `postprocess-mcmc.R` script.

As the whole procedure is computationally costly, there are three wrapper scripts to run the analysis on a SLURM scheduler as an array job (with Univa grid engine versions commented within the script). The `qsub.sh` creates the working directory if it doesn't already exist and starts the `run-mcmc.R` script. The `qsub_wrapper.sh` script is the array job that wraps around `qsub.sh`, launching one instance for every MCMC chain, whilst also submitting  the `cleanup.sh` which waits for all instances of `qsub.sh` to finish, to `postprocess-mcmc.R`. Within the postprocessing script the `postpostprocessing.py` script gets called to make easily readable tables and make some pretty figures.

The two textfiles `R_packages.txt` and `Python_packages.txt` contain all packages required for R and python respectively. The folder `example _files` contains example formatting of all the files required to run doublethink. The `configs` folder contains all config files containing all the information required to run the analysis. `scripts_not_used_by_nick` contains files containend in previous versions of doublethink that are not currently being used.

To start a doublethink analysis one has to simply supply the `qsub_wrapper.sh` script with a config file on the command line, like so:
```
./qsub_wrapper.sh configs/test.cfg
```

All the information necessary to start the analysis is contained in the config file, which is read by all bash and R scripts (which is why it is important to not change the formatting especially concerning white space around the equal signs). At the time of starting the analysis The config file will be copied into the output/working directory so you can backtrace which variables were used to create the output file even if you subsequently change the config file. When starting a new analysis with novel parameters it is important to change the `name` variable.

We have also bundled a few preprocessing scripts to get your data ready for the analysis. `data_preprocessing.py` takes a table in various input formats, encodes the data in the correct form, removes columns with too many missing values, dummifies categorical columns, removes dummies with too little presence columns. For help on how to use the script use the `--help` flag. You also have the option to impute missing data with using the [MissForest](https://rpubs.com/lmorgan95/MissForest#:~:text=MissForest%20is%20a%20random%20forest,then%20predicts%20the%20missing%20part) algorithm, which uses random forests to impute missing data. This however will take forever and there is no way of knowing how long. So use with care. The `ICD2pheno.py` file can take any ICD (or any dummy column for that matter) or list thereof and turn it into a pheno file that can be used for doublethink. If multiple columns are given, the phenotype will be 1 if any of the columns are 1 for that participants. Multiple columns should be supplied in a text file, with every column seperated by a newline.

Before explaining the config file here's a list quickly summarising all files:

- `doublethink-functions.R`: R script containing all functions necessary to run doublethink
- `run-mcmc.R`: R script that preprocesses data and starts the MCMC chains
- `make_starting_groups.py`: Python script that uses the furthest neighbour algorithm to generate starting points for the MCMC chains. Gets called from within `run-mcmc.R`
- `postprocess-mcmc.R`: R script that does postprocessing to summarise the results of the MCMC chains
- `qsub.sh`: Bash script that submits `run-mcmc.R` to the SLURM scheduler
- `qsub_wrapper.sh`: Bash script that turns qsub.sh into an array job, also launches the `cleanup.sh` script.
- `cleanup.sh`: Bash script that waits for all MCMC chains to finish to run `postprocess-mcmc.R`
- `Python_packages.txt`: All python packages required to run the doublethink.
- `R_packages.txt`: All R packages required to run the doublethink.
- `postpostprocessing.py`: Post-postprocessing script that makes tables and visuals from the results
- `cluster_posterior_probabilities.py`: Python script that get's called from within `postprocess-mcmc.R`. This uses OPTICS clustering to group 
- `ICD2pheno.py`: Turns list of ICDs into phenotype file for the analysis.
covariates together from their negative inclusion probability.
- `preprocessing`: Folder containing scripts used for data preprocessing in the COVID-19 analysis
- - `data_preprocessing.py`: Script for pre-processing data for the analysis.

## Config file

In the following I will go through every variable defined in the  config file that is required to run coublethink. I will also show an example of how to define the variable, to show the required format. All variables are seperated by newlines Brief explanations can also be found in the `test.cfg` file in the config folder. First we need to name the analysis.

The `name="test"` variable denotes the name of the output folder that will be created. For every new run you that requires seperate outputs please specify a new name.

In the next section we will define the input and output files. All file locations will need to be global paths.

`inpute_filename=<global path to doublethink folder>/example_files/input_filename.example` defines the input data. The data should be formatted as a csv file with a header row.

`columns_filename=<global path to doublethink folder>/example_files/columns_filename.example` is a csv containing the name of the columns and their respective types (float, integer or factor). The columns names here and in the input file must match.

`cases_filename=<global path to doublethink folder>/example_files/cases_filename.example` is a csv containing IDs and phenotypes (0 or 1). Columns names are specified by idcol (default "eid") and phenocol (default "pheno"). Changes to the defaults are not fully implemented.

`exclude_filename=<global path to doublethink folder>/example_files/cases_filename.example` is a text file wherein you can passids to exclude from the analysis, in case for some runs certain participants might have to be excluded. The ids are seperated by a new line

`exclude_columns=<global path to doublethink folder>/example_files/exclude_columns.example` denotes columns that should be excluded for the analysis. For example if your phenotype is made up of one or several ICD codes you may wish to exclude them from the input data. The columns are seperated by newlines.

`sourcedir=<global path to doublethink folder>` is the folder containing all the doublethink scripts (the folder you git cloned into)

`workdir=<global path to doublethink folder>/output` is where the output files will be stored and the data processing takes place. Within this folder a new folder will be created that is named after the `name` variable. If the folder already exist it will be written into.

`showcase_file=<global path to UKB showcase file>` is the path to the Data_showcase.csv file from UK Biobank this will help de-encoding categorical columns by noting which columns have which encoding. If you don't have this file you can just get rid of this parameter it won't return an error.

`encoding_file=<global path to encoding file>` is the counterpart to the `showcase_file` and contains the information about precise encodings. Again if you don't have this file simply delete the line it won't run into errors.

`R_lib_loc="/path/to/R_libs"` denotes the location of the R libraries necessary for the analysis

`python_loc="/path/to/python"` denotes the location of the python 3 binary (and by implication its libraries) necessary for the analysis

Now to the parameters that control the MCMC chains

`chain_numbers=50` controls how many seperate MCMC chains are run.

`niter=100000` controls the length of the MCMC chains.

`mu_prior=0.1` controls the prior odds of variable inclusion

`h_prior=1` controls the expected magnitude of effect sizes on a scale relative to the per-observation precision defined via the Fisher information matrix. Standard magnitudes for h_prior are 1 (implies interest in effect sizes that should be readily identified in large samples and an AIC-like criterion for model selection) and samp.size (implies interest in effect sizes that are weakly identifiable given the sample size and a BIC-like criterion for model selection)

`nlevmax=50` controls the maximum number of levels allowed for factor variables

`burnin=1000` defines how many iterations at the beginning of the MCMC chain are discarded as burn in. In order to obtain this number, one can judge an appropriate burn in from the tracing pngs created by initial runs. This is how we arrived at 1000.

`thinning=1` controls the thinning parameter of the MCMC chain.

This is all the information necessary to start the analysis. 

## Files created during the analysis
You will find all files that are being created in the `<workdir>/<name>` folder, including a copy of the config file that was used to run the analysis. 

As the analysis runs, tracing files are created that are named `<name>_chain_<chain_id>.out`. SLURM will create seperate error and out files that are similarily dumped in `<name>_chain_<chain_id>_slurm.out` and `<name>_chain_<chain_id>_slurm.err` respectively. There should therefore be three tracing files for every chain defined through the `chain_numbers` parameter. 

Alongside the human readbale tracing files there are also `doublethink-mcmc-out.ID<chain_id>.RDS` and `doublethink-mcmc-summary-out.ID<chain_id>.RDS` that trace the MCMC sampling. However these are not humanly readable and are only used for further processing. `doublethink-mcmc-summary-out.ID<chain_id>.png` is a graphic tracing of the MCMC sampling. These files can be viewed to assure that sufficient mixing has ocurred.

## Output files and their interpretation
Post-processing will take a while, as a lot of different measures are calculated. The progress of the post-processing will be recorderd in `results.cleanup.out`, which will also give an overview of all outputs being calculated. During the post-processing the following files will be created in order:

- `results.summary.tsv`: This is a table that combines all of the results tables listed below for easy readability. It picks out the top most important covariates and their most important associated stats. If you are looking at a new analysis this is where you begin making sense of the results.

- `results.plots.pdf`: This is a pdf that combines all the plots created by `postpostprocessing.py`. For easy inclusion into papers, all plots also get saved individually as .png files.

- `results.top-coefficients.tsv`: This file shows the probability of the regression coefficient beta being zero (zero.prob), which is only the case if the covariate is not included in the model. 1-zero.prob would give the posterior inclusion probablity. The file also shows the mean probability of a non-zero beta (non.zero.mean), i.e. when the covariate is included in the model and the standard deviation of the value (non.zero.sd)

- `results.posterior-inclusion-probs.tsv`: This file gives the Posterior inclusion probability (PoP) and standard error of the estimate (SE). The posterior inclusion probability is the main result of the doublethink analysis. It will give a ranking of all possible covariates in their contribution to a given phenotype (disease).

- `results.matrix-hmp-stats.tsv`: This file gives multiple values extracted from the MCMC that facilitate claculating the harmonic mean p-value.

- `results.top-log-posteriors.tsv`: This file gives the prior probability of inclusion for the covariates (prior), the log posterior odds for unclusion (mcmc) and the standard error of the value (stderr). Naturally if a flat prior is chosen the first value (prior) will be the same for all covarariates

- `results.top-log-hmp.tsv`: Here the results of the computation of the harmonic mean p-value(hmp) is displayed. The degrees of freedom (df) alongside the the logarithm to the base 10 of the hmp (log10.hmp) with its standard error (sterr). To avoid having to define interesting covariates by hand two thresholds are given, one Bayesian (Bayes.thresh) and one frequentist Bonferroni (5%Bonf.thresh) and whether the covariates surpass the threshold of significance. This gives an automated cut-off as to what covariates to interpret as risk factors for the phenotype (disease)

- `results.posterior-inclusion-cor.tsv`: This is a square matrix which will give the correlation of the posterior inclusion probability of each covariate with each other covariate. To exemplify, if we have to covariates A and B, a strong negative correlation between A and B implies, that if A is already included in the MCMC, the chain is very unlikely to also include B. This would indicate that B does not add much information above and beyond A to predict the given phenotype. The interpretation would be that A and B are inside a group encoding the same important information for the phenotype. For example BMI and weight would both encode the mass of the person and would therefore have a strong negative correlation. A strong positive correlation on the other hand would show that B add information above and beyond what A contains.

- `results.group-posterior-inclusion-probs`: This is a very simple table that lists groups of covariates (as defined through clustering their negative inclusion probabilities with OPTICS), and lists the groupwise posterior probability. The individual covariates of a group are separated by non-printable ASCII 31, coded as backslash 037.

This concludes the description of the output files. Next I will go through the pre-processing steps we undertook for the UK Biobank data so the analysis can be repeated on other dataset

## Data pre-processing
This is an ordered list of pre-processing we performed on the UK Biobank:
1. Removed patients that died or were otherwise lost to follow up

2. Removed non-English participants (as there were no COVID testing available for these)

3. Covariates with more than 10% missingness were removed

4. Covariates were typecast according to their UKB data specifications

5. Factors were encoded as dummy values (making the data very sparse)

6. ICDs present in less than 0.1% of the participants were removed to reduce the sparseness. As the ICD presence is a very long-tailed distribution.

7. Integer and continuous missing values were imputed using the mean of the column. Missing factors were cast as a new level as the integer -999, except for ICDs were missing values were interpreted as absence

8. Columns with 0 variance were removed

9. Factors with more than 50 levels were removed (you can control this with `nlevmax` in the config

10. Biomarkers which have multiple measurements for each participant (result of multiple visits to the hospital after initial registration) were reduced to one value, using mode of all values for factors and mean for continuous values.

Steps 3-10 can be performed with the `data_preprocessing.py` script. That concludes the pre-processing steps we performed on the UK Biobank data. We perform only crude imputation, so there is room for improvement there.
