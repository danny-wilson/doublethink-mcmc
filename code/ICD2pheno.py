"""
This is a script that you can feed a list of binary columns and it returns
a phenotype table that combines all positive instances across all columns
listed. Typically this would be to combine multiple ICD codes into a phenotype
table representing one disease. However, this can be used with any binary column

Make sure to include the excluded columns file created at the end in the
exclude_columns variable of the config
"""

import pandas as pd
import sys
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Make a list of binary columns into phentype")
    parser.add_argument(
        "-i",
        type=str,
        help="input datatable")
    parser.add_argument(
        "-icd",
        type=str,
        help="file with list of binary columns to create phenotype from")
    global args
    args = parser.parse_args()

    # Get the list of binary columns
    with open(sys.args.icd, "r") as inf:
        ICD_list = inf.readlines()
    # If this was created from the column encoding there might be a comma in here
    ICD_list = [x.split(",")[0] for x in ICD_list]

    # Read in data
    df = pd.read_hdf(args.i)

    # There might be noise columns in here for earlier machine learning purposes
    noise_cols = [x for x in df.columns if "noise" in x]
    exclude_list = ICD_list + noise_cols

    # get all the phenotype columns
    df = df[["eid"] + ICD_list]

    # Sum all the positive instances
    df["pheno"] = df[ICD_list].sum(axis=1)
    df = df[["eid", "pheno"]]
    phenos = df["pheno"]
    # cut the phenotype sum at 1
    phenos[phenos > 1] = 1
    df["pheno"] = phenos
    print("Here are the counts for the phenotypes:")
    print(df["pheno"].value_counts())

    # Give outfile a name
    out_title = sys.argv[1].split("_")[0]
    out_csv = out_title + "_pheno.csv"
    out_exclude = out_title + "_exclude_cols.txt",
    df.to_csv(out_csv, index=False)
    print("""Dumping phenotype csv in: {}\n Dumping columns to exlude in : {}
    \n Please include the latter in the config file under the exclude_columns
    variable""".format(out_csv, out_exclude))
    with open(out_exclude, "w+") as outf:
        outf.write("\n".join(exclude_list))
