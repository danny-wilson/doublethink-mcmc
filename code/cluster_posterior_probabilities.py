"""
This script takes the table of covariance of the posterior probabilities output
by doublethink and clusters them into groups of covariates, so their posteriors
can be combined. OPTICS is used as it can use precomputed distances, handle
uneven clusters and discard outliers. Outliers in our case are all the
covariates which are not used in the model, which is the majority of
the covariates.
"""

import argparse
import pandas as pd
import numpy as np

from sklearn.cluster import OPTICS
from sklearn.preprocessing import MinMaxScaler


def process_cov(cov_inf):
    """
    Load the covariance matrix and transform into distance matrix
    """
    # Read in matrix
    cov = pd.read_csv(cov_inf, sep="\t")
    covariate_names = cov.columns

    # If a covariate is not used in the modell all the covariances are NA
    # So drop those as they're not interesting
    covariates_not_used = cov[cov.isnull().all(axis=1)].index
    cov.drop(
        columns=covariates_not_used,
        index=covariates_not_used,
        inplace=True)

    # We are min max scaling later and the diagonal should be 0, but i'm using
    # the minimum distance here so it doesn't influence the scaling and becomes
    # 0 after the scaling
    colnames = cov.columns
    np.fill_diagonal(cov.values, np.nanmin(cov.to_numpy()))

    # NAs get created when the covariate is not used in the model, that can be
    # Maximum distance
    cov.fillna(np.nanmax(cov.to_numpy()), inplace=True)

    # Min max scale the whole thing
    #scaler = MinMaxScaler()
    #cov = scaler.fit_transform(cov)
    cov = (cov-np.nanmin(cov.to_numpy()))/(np.nanmax(cov.to_numpy())-np.nanmin(cov.to_numpy()))

    # Get the column names back in
    #cov = pd.DataFrame(cov, index=colnames, columns=colnames)

    # Make the distance matrix symmetric
    cov = (cov+cov.T)/2
    np.fill_diagonal(cov.values, 0)

    return(cov, covariate_names)


def cluster_cov(cov):
    """
    Function that uses OPTICS clustering to cluster the covariates by distance
    """

    # Initialise clusters
    clus = OPTICS(
            metric="precomputed",
            min_samples=2,
            cluster_method="xi",
            n_jobs=1)

    # cluster away
    clus.fit(cov)
    labels = clus.labels_

    return(labels)


if __name__ == "__main__":
    # This is the main function
    # Get command line arguments
    parser = argparse.ArgumentParser(
        description="Correlate posterior probability by negative covelation")
    parser.add_argument(
        "-cov", type=str, default="results.posterior-inclusion-cov.tsv",
        help="Negative covelation matrix. Default name from doublethink is results.posterior-inclusion-cor.tsv ")
    args = parser.parse_args()

    # Get the covariance matrix from the doublethink output and reformulate distance
    # matrix
    cov, covariate_names = process_cov(args.cov)

    # Use the distance to cluster the covariates
    clusters = cluster_cov(cov)

    # Dictionary registering the cluster for every covariate
    covariate2cluster = zip(cov.columns, clusters)

    # Write the clusters to file
    covariate_names = np.array(covariate_names)
    out_array = np.full(covariate_names.shape, -2)
    with open("results.clusters.out", "w+") as outf:
        for covariate, clus in covariate2cluster:
            out_array[covariate_names == covariate] = clus
            np.savetxt("results.clusters.out", out_array, delimiter=",", fmt= "%i")

