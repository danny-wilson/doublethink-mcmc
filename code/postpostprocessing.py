"""
showThis script uses follows the postprocess-mcmc.R script that runs on the .RDS
files. This script parses the output further to generate a table called
results.summary.tsv which can be viewed with excel. It also creates a pdf
containing all plots that are created. The script is simply run withing the
results folder without any command line arguments
"""
import pandas as pd
import numpy as np
import glob
import re
import sys
import glob
import pickle
import decimal
import os

import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf

from matplotlib.offsetbox import AnchoredText
from matplotlib.collections import PolyCollection
from collections import defaultdict
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import chi2
from sklearn.preprocessing import MinMaxScaler
from sklearn import metrics
from decimal import Decimal

import seaborn.categorical
seaborn.categorical._Old_Violin = seaborn.categorical._ViolinPlotter

import warnings

# This is a class so we can have our own hue colours on the violin plot


class _My_ViolinPlotter(seaborn.categorical._Old_Violin):

    def __init__(self, *args, **kwargs):
        super(_My_ViolinPlotter, self).__init__(*args, **kwargs)
        self.gray = 'white'


# initialise the violinplotter
seaborn.categorical._ViolinPlotter = _My_ViolinPlotter

import seaborn as sns

# Get own version of defaultdict where the key is returned if the key is not
# in the dictionary. Needed for the factor encoding in UK Biobank


class MyDict(dict):
    def __missing__(self, key): return key.replace("_", " ")


# ______________________________________________________________TABLE FUNCTIONS

def return_default():
    """
    Function for defaultdict, that returns a tuple of 0,0 if key is not present
    """
    return((0, 0))


class IdentityDict(defaultdict):
    def __missing__(self, key):
        return key


def PoP_table(config):
    """
    Function that reads in top posterior probabilities from
    results.top-coefficients.tsv which gets created by postprocess-mcmc.R
    Only outputs the highest 20 probabilities and outputs into the table
    """

    # Read in table
    global samp_size
    input_table, pheno_table, showcase_file, encoding_file, number_subsample, exclude_eids = read_in_config()
    beta_df = pd.read_table("results.top-coefficients.tsv")

    # Get a dictionary of beta values
    beta_zip = list(zip(
        beta_df.index,
        beta_df["non.zero.mean"].astype(float),
        beta_df["non.zero.sd"].astype(float),
        ))

    # Make into dataframe
    beta_df = pd.DataFrame(
        beta_zip,
        columns=[
            "Covariate",
            "Beta",
            "Standard Deviation Beta"])

    # Initialise defaultdict, where if key is not present in dictionary return
    # the key. This is used for encoding of factor levels. If there is no
    # factor encoding then just keep the column name as it is.
    coding_dict = MyDict()
    covariate2fieldID = defaultdict(str)
    covariate2unit = defaultdict(str)

    # Check if the Data_dictionary_showcase.tsv file is present in the config
    # Default
    coding_covariates = IdentityDict()
    
    # This file comes shipped with UKB
    if showcase_file and encoding_file:
        # Read in factor encoding that translate the integer levels to cleartext
        # description
        coding_covariates, covariate2fieldID, covariate2unit = read_in_data_coding_ukb(
            showcase_file, encoding_file, coding_dict)

    # Get a dictionary which returns beta values given the covariates
    beta_dict = defaultdict(return_default)
    for cov, beta, sd in beta_zip:
        beta_dict[coding_covariates[cov]] = (round(beta, 3), round(sd, 3))

    # Write into output string
    output_string = "___________POSTERIOR INCLUSION PROBABLITIES_____________\n\n"

    # Read in the posterior probabilities
    df = pd.read_table("results.posterior-inclusion-probs.tsv")
    df.columns = ["Posterior probability", "Standard Error"]
    df["Covariate"] = df.index

    # Here we look for the UKB field ID and unit measurement for every
    # covariate
    for cov in list(df["Covariate"]):
        # Translate covariate names to cleartext
        coded_cov = coding_covariates[cov]

        # All ICDs stem from the same column
        if cov.startswith("ICD"):
            covariate2fieldID[coded_cov] = 41202
            covariate2unit[coded_cov] = "binary"

        # So do all meds
        elif cov.startswith("dummy_MEDS"):
            covariate2fieldID[coded_cov] = 20003
            covariate2unit[coded_cov] = "binary"

        # Now other dummy variables
        elif cov.startswith("dummy"):

            # Unfortunately they were coded in a few different ways we are
            # trying out here
            versions = [coded_cov.split(":")[0],
                        coded_cov.split("  ")[0].replace("dummy ", ""),
                        " ".join(coded_cov.split(" ")[:-1])]

            # Switch to see if we found a correct version
            switch = True
            for version in versions:
                if version in covariate2fieldID:
                    covariate2fieldID[coded_cov] = covariate2fieldID[version]
                    switch = False
                    break
            # If all versions don't work just set to NA
            if switch:
                covariate2fieldID[coded_cov] = "NA"

            covariate2unit[coded_cov] = "binary"

    # Get a dictionary of all covariates and their posterior probabilities
    pop_dict = dict(zip(df["Covariate"].map(
        coding_covariates), df["Posterior probability"].round(3).astype(str)))
    df = df.fillna(0)

    # Get the beta values in
    df = df.merge(beta_df, on="Covariate", how="left")

    # Read the top 20 columns from the raw data for the contingency tables,
    # violin plots and variance
    data_df = get_input_data(df["Covariate"], input_table, pheno_table, exclude_eids)

    # Numbers of cases and controls
    phenotype_counts = data_df['Phenotype'].value_counts()
    number_cases = phenotype_counts.get(1,0)
    number_ctrls = phenotype_counts.get(0,0)
    if number_cases + number_ctrls != len(data_df['Phenotype']):
        raise Exception("Unexpected value in phenotype file")

    # We need the number of samples to calculate the p-value later. So here
    # we check if we subsampled and if not take all participants
    samp_size = number_cases + number_ctrls
    if number_subsample != 0:
        samp_size = number_cases + min(number_subsample, number_ctrls)

    print("Sample size:        {}".format(samp_size))
    print("Number of cases:    {}".format(number_cases))
    print("Number of controls: {}".format(number_ctrls))

    # Get the variance of the columns
    variance = data_df.var()

    # Make into DataFrame
    variance = pd.DataFrame.from_dict(
                {
                    "Covariate": list(variance.index),
                    "Variance": variance
                }
                )

    # Merge with the posterior probabilities
    df = df.merge(variance, on="Covariate", how="left")

    # Calculate the p-value from the posterior probability and the number of
    # samples
    df["p-value"] = [calc_p_value(float(x),
                                  float(samp_size),
                                  float(config['h']),
                                  float(config['nu']),
                                  float(config['mu']))
                     for x in df["Posterior probability"]]

    # Translate covariates to cleartext
    df["Covariate"] = df["Covariate"].map(coding_covariates)

    # Get units for the covariates in
    df["Unit"] = df["Covariate"].map(covariate2unit)

    # Get UKB field ID into the table
    df["Field ID"] = df["Covariate"].map(covariate2fieldID)

    # Get the proportion of phenotypical variance explained
    df["Proportion of Phenotypic Variance Explained"] = (
        df["Standard Deviation Beta"]**2 + df["Beta"]**2) * df["Variance"]

    # Round to 3 decimals
    df = df.round(3)

    # Translate factor levels into cleartext
    data_df.columns = [coding_covariates[x] for x in data_df.columns]

    # Re-index
    df.index = df["Covariate"]

    # Sort by posterior probabilities
    df = df.sort_values("Posterior probability", ascending=False)

    # Return the string that will be the output table, the posterior probability
    # dict, the beta dict, the dataframe with the raw data and the translation
    # of the factor levels
    return(pop_dict, df, beta_dict, data_df, coding_covariates, pheno_table)


def calc_p_value(PoP, n, h, nu, mu):
    """
    Function that calculates the *adjusted* p-value from posterior probability and parameters
    """

    posterior_odds = np.float64(PoP) / (1.0 - np.float64(PoP))
    chi2_stat = 2*np.log(posterior_odds / (nu * mu * np.sqrt(h / (n + h))))
    p_value = -1*np.log10(1-chi2.cdf(chi2_stat, 1))
    if p_value > (-1*np.log10(0.2)):
        return(round(p_value, 2))
    else:
        return("-")


def get_prediction():
    """
    Function that collects the predicted outputs of the logistic regression of
    every chain. It averages over all predictions to arrive at the result.
    """

    # Get all predictions
    prediction_files = glob.glob("doublethink-prediction.ID*.csv")
    prediction_files = [
        pd.read_csv(
            pred, index_col=0, header=0, names=["eid", "pheno_pred_chain"])
        for pred in prediction_files]

    # Get predictions together
    start_df = prediction_files[0]
    for df in prediction_files:
        start_df = start_df.merge(df, how="outer", on="eid")

    # Average prediction
    start_df.index = start_df["eid"]
    start_df.drop(columns="eid", inplace=True)
    start_df["pheno_predict"] = list(start_df.median(axis=1))
    start_df.reset_index(inplace=True)
    start_df = start_df[["eid", "pheno_predict"]]
    return(start_df)


def read_in_data_coding_ukb(showcase_file, encoding_file, coding_dict):
    """
    Function that reads in the Data_dictionary_showcase file from UKB and the
    Codings.tsv file to translate the integer labels of the factor levels
    to cleartext
    """

    # Return emtpy list if key not present
    coding2covariate = defaultdict(list)

    # Read in showcase files
    df = pd.read_csv(
        showcase_file,
        usecols=[
            "FieldID",
            "Units",
            "Field",
            "Coding"])
    df["Field"] = df["Field"].map(
        lambda x: re.sub(
            "[^0-9a-zA-Z]",
            " ",
            x.replace(
                "(",
                "").replace(
                ")",
                "")))

    # First codes as integers then strings
    df["Coding"] = df["Coding"].fillna(0).astype(int).astype(str)
    covariate2fieldID = dict(
        zip(df["Field"], df["FieldID"].fillna(0).astype(int)))
    covariate2unit = dict(zip(df["Field"], df["Units"].fillna("none")))

    # Fill the coding dictionary with which covariate is encoded with what coding
    # not only the levels of factors are giving numbers but also the codings
    # itselfs are encoded by numbers, so what translation to use for each column
    # So for example field pain types experienced uses coding 5 and in coding
    # 5 1 refers to no pain 2 refers to headaches and so on
    for coding, covariate in zip(df["Coding"], df["Field"]):
        coding2covariate[coding].append(covariate)

    # As the codings itself are given numbers, now note what the individual
    # codings are and save them at the integer that references their encoding
    # So coding 5 stands for 1 no pain, 2 headaches and so on. This is a dictionary
    # of dictionaries where the highest levels is the integer encoding of the codings
    # and beneath are the integer encodings of the factor levels
    with open(encoding_file) as inf:
        for line in inf:
            coding, integer, cleartext = line.strip().split("\t")
            # sometimes the codings are empty so skip

            try:
                float(integer)
            except ValueError:
                continue

            # negative level encodings are always didn't answer or wouldnt tell
            if float(integer) < 0:
                integer = "0"
                cleartext = "No answer"
            if coding not in coding_dict:
                coding_dict[coding] = {}
            if coding not in coding2covariate:
                continue

            # Get the covariates that are subject to the current encoding
            covariates = coding2covariate[coding]
            for covariate in covariates:
                # clean up the variable name
                covariate_clean = re.sub(
                    "[^0-9a-zA-Z]",
                    "_",
                    covariate.replace(
                        "(",
                        "").replace(
                        ")",
                        ""))

                covariate_clean = re.sub("_+", "_", covariate_clean)

                # Get them into the format we're used to
                covariate_coded = "dummy_{}_{}0".format(
                    covariate_clean, integer)

                covariate2fieldID[covariate_coded] = covariate2fieldID[covariate]
                covariate2unit[covariate_coded] = covariate2unit[covariate].capitalize(
                )

                # Save in dictionary
                coding_dict[covariate_coded.replace(
                    "__", "_")] = "{}: {}".format(covariate, cleartext)

    # Return the dictionary
    return(coding_dict, covariate2fieldID, covariate2unit)


def row_wrangler(row, single=False):
    """
    Function that takes a row of the posterior probability table and put it
    into the right format for the output table that gets written to disk
    """

    pop, pop_se, cov, beta, beta_se, var, pvalue, unit, fieldID, var_expl = row.tolist()
    pop = "{:.2f}({:.2f})".format(pop, pop_se)
    beta = "{:.2f}({:.2f})".format(beta, beta_se)

    # If it is a singular covariate
    if single == True:
        string = "{}\t\t\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            cov,
            pop,
            beta,
            var_expl,
            pvalue,
            unit,
            fieldID
        )

    # If it is a group of covariates
    else:
        string = "\t{}\t\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            cov,
            pop,
            beta,
            var_expl,
            pvalue,
            unit,
            fieldID
        )

    return(string)


def group_PoP(output_string, pop_df, pop_dict, beta_dict, coding_covariates, config):
    """
    Get the groups of posterior probabilities that are created by
    cluster_posterior_probabilities.py that gets called within the
    postprocess-mcmc.R script. Then add these to the table
    """

    groups = []
    group_pop_dict = []
    covariates_in_groups = []

    # Read out the grouped posterior probabilities
    with open("results.group-posterior-inclusion-probs.tsv") as inf:
        for line in inf.readlines()[1:]:

            # Get group and posterior probability
            group_unsplit, PoP = line.strip().split("\t")
            PoP = round(float(PoP), 3)

            # Split the group into individuals: delimiter defined postprocess-mcmc.R line 280
            group = group_unsplit.split("\037")

            # Trranslate factor levels to cleartext
            covariates_in_groups += group
            group = [coding_covariates[x] for x in group]
            group = [(x, float(pop_dict[x])) for x in group]

            group_pop_dict.append((group_unsplit, PoP))
            if len(group) > 30:
                continue

            # Sort the group by posterior probabilities
            group = sorted(group, key=lambda x: x[1], reverse=True)
            groups.append([float(PoP), [pop_df.loc[x[0]] for x in group]])
            # Skip clusters of covariates that are bigger than 20. These are
            # usually outgroups

    # We need this to figure out covariates that are not in groups
    covariates_in_groups = set(covariates_in_groups)

    # Interested in all variables and groups
    for cov in pop_df["Covariate"]:
        if cov not in covariates_in_groups:
            row = pop_df.loc[cov]
            groups.append([float(row["Posterior probability"]), [row]])

    # Sort the groups by the groupwise posterior probability
    group_pop_dict = sorted(group_pop_dict, key=lambda x: x[1], reverse=True)

    # Only output the biggest 10 groups * for plotting purposes *
    group_pop_dict = dict(group_pop_dict[:10])

    # Sort by posterior probability
    groups = sorted(groups, key=lambda x: x[0], reverse=True)

    output_string += """Group\tCovariate\tGroup posterior probability\tIndividual Posterior probability (s.e.)\tBeta(s.e.)\tVariance explained\t-log10 adjusted p-value\tUnits\tFieldID\n"""
    for i, group in enumerate(groups):
        if len(group) == 1:
            output_string += row_wrangler(group[0][1], True)
        else:
            output_string += "{}\t\t{:.2f}\t\t\t\t{}\n".format(
                i+1, group[0], calc_p_value(float(group[0]),
                                            float(samp_size),
                                            float(config['h']),
                                            float(config['nu']),
                                            float(config['mu'])))
            for g in group[1]:
                output_string += row_wrangler(g)

    # Return the string that will be output table and a dictionary containing
    # all the groups
    return(output_string, group_pop_dict)


def pop_cov(output_string, pop_dict, coding_covariates):
    """
    Function that reads in the covariance of the posterior probabilities from
    results.posterior-inclusion-cov.tsv
    """

    output_string += "\n___________POSTERIOR INCLUSION COVARIANCE____________\n\n"

    # Read in the table
    df = pd.read_table("results.posterior-inclusion-cov.tsv")

    # Covariates that are never included in the logistic regression, have
    # a PoP of NA and covariance of NA. We can throw those out
    covariates_not_used = df[df.isnull().all(axis=1)].index
    df.drop(
        columns=covariates_not_used,
        index=covariates_not_used,
        inplace=True)

    # Make a deep copy of the dataframe so we can use that later for the plot
    df_cov = df.copy()
    df_cov.columns = [coding_covariates[x] for x in df_cov.columns]
    df_cov.index = [coding_covariates[x] for x in df_cov.index]
    df = pd.DataFrame(np.tril(df.values), columns=df.columns, index=df.index)

    # Collect the covariates with the lowest (so highest negative) covariance
    # Those are the most closely related covariates
    top_correlations = []
    for col in df.columns:
        min_cov = df[col].min()
        min_cov_col = df[col].idxmin(skipna=True)
        top_correlations.append(
            (min_cov, coding_covariates[min_cov_col], coding_covariates[col]))

    # Sort the top negative covariates
    top_correlations = sorted(top_correlations, key=lambda x: x[0])[:25]

    # Write all the top negative covariance into the file
    output_string += "Covariate 1\tPoP 1\tCovariance\tCovariate 2\tPoP 2\n"
    for cov in top_correlations[:15]:
        cov_val, covariate_1, covariate_2 = cov
        cov_val = str(round(cov_val, 2))
        pop1 = pop_dict[covariate_1]
        pop2 = pop_dict[covariate_2]
        output_string += "\t".join([covariate_1,
                                    pop1, cov_val, covariate_2, pop2])
        output_string += "\n"

    # Return the table, the top correlations and the full dataframe of
    # covariance
    return(output_string, top_correlations, df_cov)


def read_in_config():
    """
    Read in the config file. This is the same config file that doublethink
    was started with. If there are no encoding files we can skip the encoding
    step
    """

    # Find the file that ends on .cfg
    input_cfg = glob.glob("*.cfg")[0]

    # If we can't find encoding files then ignore the encoding step
    encoding_file = False
    showcase_file = False
    # Default
    number_subsample = 0
    with open(input_cfg) as inf:
        for line in inf:
            # the input table
            if line.startswith("input_filename"):
                input_table = line.strip().split("=")[1].replace("\"", "")
            # the phenotype file
            elif line.startswith("cases_filename"):
                pheno_table = line.strip().split("=")[1].replace("\"", "")
            # the Data_dictionary_showcase.tsv file
            elif line.startswith("showcase_file"):
                showcase_file = line.strip().split("=")[1].replace("\"", "")
            # the encoding of the factor levels
            elif line.startswith("encoding_file"):
                encoding_file = line.strip().split("=")[1].replace("\"", "")
            elif line.startswith("controls_subsample"):
                number_subsample = int(line.strip().split("=")[1])
            # the exclude_filename file
            elif line.startswith("exclude_filename"):
                exclude_eids = line.strip().split("=")[1].replace("\"", "")

    # returns the 5 variables
    return(input_table, pheno_table, showcase_file, encoding_file, number_subsample, exclude_eids)


def get_input_data(top_posteriors, input_table, pheno_table, exclude_eids):
    """
    Read in the raw data from the top 20 covariates. We need this to plot
    violinplots and contingency table as well as calculating the variance.
    """

    # Read in the raw table
    df = pd.read_csv(
        input_table,
    )
    for col in df.columns:
        if col == "eid":
            continue

    df_pheno = pd.read_csv(pheno_table)
    df_pheno.columns = ["eid", "Phenotype"]
    # Binarize: 0 if less than median, 1 if greater than median
    pheno_original = df_pheno['Phenotype'].copy()
    pheno_median = df_pheno['Phenotype'].median()
    df_pheno['Phenotype'] = (df_pheno['Phenotype'] > pheno_median).astype(int)
    if not (pheno_original == df_pheno['Phenotype']).all():
        warnings.warn("get_input_data: df_pheno['Phenotype'] column was altered by binarization.")
    # Left merge
    df = df.merge(df_pheno, on="eid", how="left")
    
    # Remove the EIDs excluded
    eids = list(map(np.int64, open(exclude_eids).read().splitlines()))
    df_gd = df.loc[df['eid'].isin(eids) == False]

    return(df_gd)

# ___________________________________________________________PLOTTING FUNCTIONS


def plot_predictive_performance(prediction_table, pheno_table, pdf):
    """
    Function to plot the predictive performance of the logistic regression

    """

    df_pheno = pd.read_csv(pheno_table)
    df_pheno["eid"] = df_pheno["eid"].astype(int)
    # Binarize
    pheno_original = df_pheno["pheno"].copy()
    pheno_median = df_pheno["pheno"].median()
    df_pheno["pheno"] = (df_pheno["pheno"] > pheno_median).astype(int)
    if not (pheno_original == df_pheno["pheno"]).all():
        warnings.warn("plot_predictive_performance: df_pheno['pheno'] column was altered by binarization.")

    prediction_table["eid"] = prediction_table["eid"].astype(int)
    prediction_table = prediction_table.merge(
        df_pheno, on="eid", how="left")

    # Get area under the ROC
    false_positive_rate, true_positive_rate, thresholds = metrics.roc_curve(
        prediction_table["pheno"], prediction_table["pheno_predict"])
    precision, recall, thresholds = metrics.precision_recall_curve(
        prediction_table["pheno"], prediction_table["pheno_predict"])
    fscore = (2 * precision * recall) / (precision + recall)
    gmeans = np.sqrt(true_positive_rate * (1-false_positive_rate))
    youdenJ = true_positive_rate - false_positive_rate

    decision_threshold = thresholds[np.argmax(fscore)]
    print("Decision threshold: {}".format(decision_threshold))

    prediction_table["pheno_predict_rounded"] = 0
    prediction_table["pheno_predict_rounded"] = prediction_table['pheno_predict_rounded'].where(
        prediction_table['pheno_predict'] >= decision_threshold, 1)
    print(prediction_table["pheno_predict_rounded"].value_counts())
    # Plot prediction
    fig, axes = plt.subplots(3, 1)
    sns.histplot(
        x=prediction_table["pheno"],
        bins=np.linspace(0, 1, num=100),
        ax=axes[0]
    )
    axes[0].title.set_text("True Phenotypes")
    sns.histplot(
        x=prediction_table["pheno_predict"],
        bins=np.linspace(0, 1, num=100),
        ax=axes[1],
    )
    axes[1].title.set_text("Predicted Phenotypes")
    sns.histplot(
        x=prediction_table["pheno_predict_rounded"],
        bins=np.linspace(0, 1, num=100),
        ax=axes[2],
    )
    axes[2].title.set_text("Predicted Phenotypes Rounded")
    plt.tight_layout()
    # pdf.savefig(fig)
    plt.close()

    auc = metrics.roc_auc_score(
        prediction_table["pheno"],
        prediction_table["pheno_predict"])

    # Get confusion matrix
    TN, FP, FN, TP = metrics.confusion_matrix(
        prediction_table["pheno"],
        prediction_table["pheno_predict_rounded"]).ravel()

    # Sensitivity, hit rate, recall, or true positive rate
    TPR = TP/(TP+FN)
    # Specificity or true negative rate
    TNR = TN/(TN+FP)
    # Precision or positive predictive value
    PPV = TP/(TP+FP)
    # Negative predictive value
    NPV = TN/(TN+FN)
    # Fall out or false positive rate
    FPR = FP/(FP+TN)
    # False negative rate
    FNR = FN/(TP+FN)
    # False discovery rate
    FDR = FP/(TP+FP)
    # Overall accuracy
    ACC = (TP+TN)/(TP+FP+FN+TN)

    print("Here is a short accuracy report")
    print(metrics.classification_report(
        prediction_table["pheno"],
        prediction_table["pheno_predict_rounded"]))
    # Plot all scorers
    scorers = [ACC, TPR, TNR, PPV, NPV, FPR, FNR, FDR]
    scorer_names = [
        "Accuracy",
        "True positive rate",
        "True negative rate",
        "Positive predictive value",
        "Negative predictive value",
        "False positive rate",
        "False negative rate",
        "False discovery rate"]
    y_positions = range(len(scorers))
    cmap = matplotlib.cm.get_cmap('tab10')
    colour = np.linspace(0, 1, num=len(scorers))
    colour = [cmap(x) for x in colour]

    fig = plt.figure()
    ax = plt.gca()
    ax.barh(
        y_positions,
        scorers,
        align='center',
        color=colour,
        alpha=0.6)
    xticks = np.arange(0, 1.1, 0.1)
    ax.xaxis.set_ticks(xticks)
    plt.xticks(xticks, rotation=90)
    ax.set_yticks(y_positions, labels=scorer_names)
    plt.tight_layout()
    # pdf.savefig(fig)
    plt.close()

    # Plot ROC curve
    fig = plt.figure()
    sns.lineplot(
        x=false_positive_rate,
        y=true_positive_rate,
        )
    plt.plot([0, 1], [0, 1], ls="--", c="silver")
    # plt.axhline(true_positive_rate[np.argmax(fscore)], ls="dotted", c="dodgerblue")
    # plt.axvline(false_positive_rate[np.argmax(fscore)], ls="dotted",
    # c="dodgerblue")
    label = "AUC = {:.2f}\n".format(auc)
    text_box = AnchoredText(label, frameon=True, loc=4, pad=0.5)
    plt.setp(text_box.patch, facecolor='white', alpha=0.5)
    plt.gca().add_artist(text_box)
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    pdf.savefig(fig)

    return(pdf)


def plot_posteriors(pdf, pop_dict, pop_df, beta_dict, data_df):
    """
    Function that takes the table of top posteriors created by PoP_table and
    plots a barchart of posteriors, a plot of betas, a barchart of explained variance
    and a violin plot of the data distribution
    """
    # Y positions for the plots
    pop_df["y"] = range(pop_df.shape[0])

    # Colourmap so every covariate has its own colour
    cmap = matplotlib.cm.get_cmap('jet')
    colour = np.linspace(0, 1, num=pop_df.shape[0])
    colour = [cmap(x) for x in colour]

    # Make figure with 4 columns
    fig, axes = plt.subplots(
        nrows=1,
        ncols=4,
        sharey=True,
        figsize=(15, 15),
        gridspec_kw={'width_ratios': [4, 2.5, 2.5, 5]}
        )

    # Horizontal barplot of the posterior probabilities use standard errors for
    # errorbars on the first column of the plot
    axes[0].barh(
        pop_df["y"],
        pop_df["Posterior probability"].astype(float),
        xerr=pop_df["Standard Error"],
        align='center',
        alpha=0.5,
        color=colour)

    # Format labels so that underscores are spaces again
    labels = [x.replace("_", " ").replace(": ", "\n")
              for x in pop_df["Covariate"]]

    axes[0].set_yticks(pop_df["y"], labels=labels)
    axes[0].title.set_text("Top Posterior Probabilities")

    # Plot the betas with the associated standard deviation on the second
    # column
    pop_df = pop_df.fillna(0)
    axes[1].errorbar(x=pop_df["Beta"],
                     y=pop_df["y"],
                     xerr=pop_df["Standard Deviation Beta"],
                     ecolor=colour,
                     linestyle="",
                     elinewidth=3,
                     alpha=0.5)
    axes[1].scatter(x=pop_df["Beta"],
                    y=pop_df["y"],
                    color=colour)

    # Plot zero line so the directionality can be more easily established
    axes[1].axvline(0, color="grey", ls="--")

    # Put labels at bottom
    axes[1].tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left=False,        # ticks along the bottom edge are off
        right=False,       # ticks along the top edge are off
        labelright=False)  # labels along the bottom edge are off
    axes[1].title.set_text("Directionality of Effect")

    # Plot proportion variance explained in the third column
    axes[2].barh(
        pop_df["y"],
        pop_df["Proportion of Phenotypic Variance Explained"].astype(float),
        # xerr=pop_df["Standard Deviation Proportion of Phenotypic Variance
        # Explained"].astype(float),
        align='center',
        alpha=0.5,
        color=colour)
    axes[2].title.set_text("Proportion of \n variance explained")

    # Make a seperate dataframe for the violinplots as seaborn requires a weird
    # format
    x_violin = []
    y_violin = []
    pheno_violin = []
    for i, col in enumerate(pop_df["Covariate"]):
        if col == "eid":
            continue
        # Min max scale the data as otherwise we can't use the same x-axis for
        # the violinplots
        minmax_scaled = data_df[col].sub(
            data_df[col].min()).div(
            (data_df[col].max() - data_df[col].min()))

        x_violin += list(minmax_scaled)
        y_violin += [str(i+1)] * data_df.shape[0]
        pheno_violin += list(data_df["Phenotype"])

    # Combine all columns into dataframe
    df_violin = pd.DataFrame(
        list(zip(x_violin, y_violin, pheno_violin)),
        columns=["Value", "Covariate", "Phenotype"])

    # Order them by the correct y order
    order = list(map(str, range(1, len(labels)+1)))

    # Make a violinplot which is split by phenotype
    sns.violinplot(
        data=df_violin,
        x="Value",
        y="Covariate",
        hue="Phenotype",
        split=True,
        bw=0.2,
        orient="h",
        gridsize=1000,
        cut=0,
        scale_hue=True,
        scale="area",
        ax=axes[3],
        linecolor="white",
        hue_order=[1, 0],
        inner=None,
        palette=[".5", ".8"],
        order=order
        )

    # Change the colour of the histogram to resemble the other columns
    colour = [val for val in colour for _ in (0, 1)]
    for i, violin in enumerate(axes[3].findobj(PolyCollection)):
        if i % 2 != 0:
            rgb = 0.5 + 0.4 * np.array(colour[i])  # make whiter
        else:
            rgb = colour[i]
        violin.set_facecolor(rgb)

    axes[3].set_yticks(pop_df["y"], labels=labels)
    axes[3].set_ylabel("")

    # Turn all the axis labels off
    axes[3].tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left=False,        # ticks along the bottom edge are off
        right=False,       # ticks along the top edge are off
        labelright=False)  # labels along the bottom edge are off
    axes[3].title.set_text("Violin Plots Data Distribution")
    axes[3].legend(loc='lower right')
    plt.tight_layout()
    plt.savefig("results.top-posteriors-barplot.png")
    pdf.savefig(fig)

    return(pdf)


def plot_contingency(pdf, data_df, pop_df, pop_dict, coding_covariates):
    """
    Function that takes the binary columns within the top posterior probabilities
    and plots the contingency tables as heatmaps
    """
    # Check which columns are binary
    binary_cols = [
        col for col in data_df if np.isin(
            data_df[col].unique(), [
                0, 1]).all()]

    # Not phenotype though
    binary_cols.remove("Phenotype")
    binary_cols = set(pop_df["Covariate"]).intersection(set(binary_cols))

    # Get the posterior probabilities for sorting
    binary_cols = [(x, float(pop_dict[x])) for x in binary_cols]

    # Sort by posterior probability
    binary_cols = sorted(binary_cols, key=lambda x: x[1], reverse=True)
    contingencies = []
    contingencies_annot = []

    # Make contingency tables normalised within columns for colour purposes
    for col, pop in binary_cols:
        contingency = pd.crosstab(
            data_df[col].astype(int),
            data_df["Phenotype"],
            normalize="columns",
            )
        contingency = contingency.applymap(lambda x: 1-x)
        contingency_annot = pd.crosstab(
            data_df[col].astype(int),
            data_df["Phenotype"],
            )

        # Get the right tables in
        contingency.index = [
            "{} - 0".format(coding_covariates[col]),
            "PoP = {} - 1".format(pop)]
        contingencies.append(contingency)
        contingencies_annot.append(contingency_annot)

    # Make all contingencies into one long dataframe
    contingencies = pd.concat(contingencies, axis=0)
    contingencies_annot = pd.concat(contingencies_annot, axis=0).astype(int)

    # Make line between different covariates
    lines = list(range(0, len(contingencies)*2, 2))
    fig = plt.figure(figsize=(10, 15))

    # Make the heatmap
    sns.heatmap(
        contingencies,
        annot=contingencies_annot,
        # square=True,
        cmap="coolwarm",
        cbar=False,
        fmt="g",
        linewidth=.5)

    # Plot the divider lines
    for line in lines:
        plt.axhline(line, xmin=0, xmax=2, c="white", lw=3)
    ax = plt.gca()
    ax.set_xlabel('Phenotype')
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    plt.tight_layout()
    plt.savefig("results.top-posteriors-contingencies-heatmap.png")

    pdf.savefig(fig)

    return(pdf)


def plot_group(pdf, group_pop_dict, pop_dict, coding_covariates):
    """
    Function that takes the top scoring groups of covariates and plots them
    as grouped barcharts
    """

    # Categorical colours
    cmap = matplotlib.cm.get_cmap('tab10')
    colour = np.linspace(0, 1, num=10)
    colour = [cmap(x) for x in colour]

    # Sort groups by posterior probability
    groups = sorted(group_pop_dict.items(), key=lambda x: x[1], reverse=True)
    x = []
    x_group = []
    titles = []
    colours = []
    counter = 0
    colour_counter = 0
    y_mins = []
    y_maxs = []

    # Iterate through groups to collect all the values necessary for the plot
    for group, group_pop in groups:
        # When we have 10 groups stop
        if colour_counter == 10:
            break

        # Individual members of the group: delimiter defined postprocess-mcmc.R line 280
        group_split = group.split("\037")

        # Translate factor levels to cleartext
        group_split = [coding_covariates[x] for x in group_split]
        # Get individual posteriors
        group_split = [(x, pop_dict[x]) for x in group_split]
        # Sort by posterior probability within group
        group_split = sorted(group_split, key=lambda x: x[1], reverse=True)

        # If group bigger than 15 skip, because these are usually outliers
        if len(group_split) > 15:
            continue

        # Adjust the boundaries of the outer plot
        y_mins.append(counter - 0.5)

        # Get posterior probability of individual group members
        for cov, pop in group_split:
            x.append(float(pop))
            titles.append(cov)
            colours.append(colour[colour_counter])
            # This is for the position in the plot
            counter += 1

        # Append x value for the outer plot
        x_group.append(float(group_pop))
        y_maxs.append(counter-0.5)

        # Add an empty row in the plot to seperate the groups
        x.append(0)
        titles.append("")
        colours.append("white")
        colour_counter += 1
        counter += 1

    y = range(len(titles))
    fig, ax = plt.subplots(figsize=(10, 6))
    counter = 0

    # Make the outer group plot
    for x1, y1, y2 in zip(x_group, y_mins, y_maxs):
        ax.fill_between(
            x=[0, x1],
            y1=[y1, y1],
            y2=[y2, y2],
            color=colour[counter],
            alpha=0.25
        )
        counter += 1

    # Make the inner barplots for individual covariates
    ax.barh(
        y,
        x,
        align='center',
        color=colours,
        alpha=0.6)
    titles = [x.replace("_", " ") for x in titles]

    xticks = np.arange(0, 1.1, 0.1)
    ax.xaxis.set_ticks(xticks)
    plt.xticks(xticks, rotation=90)
    ax.set_yticks(y, labels=titles)

    # Remove change labelsize
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        right=False,       # ticks along the bottom edge are off
        left=False,        # ticks along the top edge are off
        labelleft=True,
        labelsize=7,
        )  # labels along the bottom edge are off
    ax.invert_yaxis()
    plt.title("Grouped Posterior Probabilities")
    plt.tight_layout()
    plt.savefig("results.group-posteriors-barplot.png")
    pdf.savefig(fig)

    return(pdf)


def plot_cov(pdf, cov_df, top_correlations, pop_dict, beta_dict):
    """
    Function that takes the covariance of posterior inclusion and plots a
    clustermap. This is a heatmap, where columns and rows are clustered by
    correlation
    """

    # Get the covariances from the most correlated covariates
    _, covs1, covs2 = zip(*top_correlations[:15])
    covs = list(covs1) + list(covs2)
    covs = set(covs)

    # Get only the top interactions from the covariance matrix
    cov_df = cov_df[covs].T
    cov_df = cov_df[covs]

    # Colourmap for the posterior probability, we want to plot this on the top
    cmap = matplotlib.cm.get_cmap('Greys')

    # Colourmap for betas we want to plot this on the side
    cmap2 = matplotlib.cm.get_cmap('RdYlGn')

    # get the posterior probabilities
    pops = [float(pop_dict[x]) for x in cov_df.columns]

    # Get the betas
    betas = [float(beta_dict[x][0]) for x in cov_df.columns]

    # Get colours
    colours = [cmap(x) for x in pops]
    colours2 = [cmap2(x) for x in betas]

    # Minmax scale the covariance for correct clustering
    df = MinMaxScaler().fit_transform(cov_df)

    # Translate factor levels
    covs = [coding_covariates[x].replace("_", " ") for x in covs]
    df = pd.DataFrame(df, columns=covs, index=covs)

    # Make sure the distance matrix is square and symmetric
    df = (df+df.T)/2
    np.fill_diagonal(df.values, 0)

    # Calculate linkage
    lnkg = linkage(squareform(df), method='average')

    # Plot clustermap
    c = sns.clustermap(
        cov_df,
        cmap="bwr",
        row_linkage=lnkg,
        col_linkage=lnkg,
        row_colors=colours,
        col_colors=colours2,
        square=True
        )

    # Disable colourbar
    c.cax.set_visible(False)
    ax = c.ax_heatmap

    ax.set_xticklabels(covs)
    ax.set_yticklabels(covs)
    ax.tick_params(
        axis='both',
        which='both',
        labelsize=7,
        right=False,
        bottom=False)
    plt.ylabel('Posterior Probability')
    plt.tight_layout()
    plt.title("Clustermap of Top Covariances")
    fig = plt.gcf()
    plt.savefig("results.posterior-inclusion-cov-clustermap.png")
    pdf.savefig(fig)

    return(pdf)


def read_config(name):
    """
    Read in the config file. This is the same config file that doublethink
    was started with. Read all the variables
    """

    # Find the file that ends on .cfg
    input_cfg = name + '.cfg'

    config = {}
    exec(open(input_cfg).read(), config)

    # Compute hyper-parameters, including maximum number of variables (nu)
    with open(config['columns_filename'], "rbU") as f:
        nu = sum(1 for _ in f)
    config['nu'] = nu-1                                  # Maximum number of variables
    config['mu'] = config['mu_prior']                    # Prior odds of variable inclusion
    config['h'] = config['h_prior']                      # Scale parameter for coefficients

    return(config)


if __name__ == "__main__":
    """
    Main function
    """
    name = os.getcwd().split("/")[-1]
    
    # Read in the config file
    config = read_config(name)

    # This will become the output table
    output_string = ""

    # This is the pdf we drop all our plots in
    pdf = matplotlib.backends.backend_pdf.PdfPages(
        "results.{}.plots.pdf".format(name))

    # Read in posterior probabilities
    pop_dict, pop_df, beta_dict, data_df, coding_covariates, pheno_table = PoP_table(config)
    # Read in predictions from the logistic regression
    prediction_table = get_prediction()

    # Plot predictive performance
    pdf = plot_predictive_performance(prediction_table, pheno_table, pdf)

    # Plot posterior probabilities
    pdf = plot_posteriors(pdf, pop_dict, pop_df.head(10), beta_dict, data_df)
    # Plot contingency table
    #pdf = plot_contingency(
    #    pdf, data_df, pop_df.head(10),
    #    pop_dict, coding_covariates)

    # Get groupwise posterior probabilities
    output_string, group_pop_dict = group_PoP(
        output_string, pop_df, pop_dict, beta_dict, coding_covariates, config)

    # Plot the groupwise posterior probabilities
    pdf = plot_group(pdf, group_pop_dict, pop_dict, coding_covariates)

    # Get covariances
    output_string, top_correlations, df_cov = pop_cov(
        output_string, pop_dict, coding_covariates)

    # Plot covariances
    pdf = plot_cov(pdf, df_cov, top_correlations, pop_dict, beta_dict)

    # Write the table to file
    with open("results.{}.summary.tsv".format(name), "w+") as outf:
        outf.write(output_string.replace("_", " "))
    pdf.close()
