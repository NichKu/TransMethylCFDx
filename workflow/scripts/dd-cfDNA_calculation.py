"""
dd-cfDNA Calculation Script

This script calculates the donor-derived cell-free DNA (dd-cfDNA) fraction from pileup files.
It supports processing multiple samples, applying quality control and generating summary statistics.
The script can also use recipient genomic information for samples with %dd-cfDNA above 20%.

Features:
- Parses genomic mapping and pileup files.
- Calculates dd-cfDNA fractions with or without recipient genomic information.
- Applies quality control, including outlier detection and filtering.
- Supports Gaussian Mixture Model (GMM) fitting for noise adjustment and donor genotype identification.
- Generates plots for visualizing results.
- Outputs results in CSV format.

Dependencies:
- Python 3.10 or higher
- pandas
- matplotlib
- scikit-learn

Usage:
Run the script with the required arguments. For example:

    python dd-cfDNA_calculation.py -dir <input_directory> -o <output_file> -c <min_coverage>

    -h or --help for more information on available options.


Author: Nicholas Kueng
Date: 2024-10-04

License: GNU General Public License v3.0
"""

import argparse
import glob
import os
import warnings
import pandas as pd
import matplotlib.pyplot as plt
from statistics import mean
from sklearn.mixture import GaussianMixture
warnings.filterwarnings("ignore")

def parse_genomic_mapping_file(directory: str, file_ending: str):
    """
    Parses the genomic mapping file in the specified directory.
    :param directory: Directory containing the genomic mapping file.
    :param file_ending: The ending of the genomic mapping file.
    :return: DataFrame containing the genomic mapping data or an error message.
    """

    try:
        # Find the file
        file_list = glob.glob(os.path.join(directory, file_ending))

        # Check if there is exactly one file with the given ending
        if len(file_list) == 1:
            file_path = file_list[0]
            # Read the CSV file using Pandas
            genomic_mapping = pd.read_csv(file_path)
            return genomic_mapping

        elif len(file_list) > 1:
            return f"Error: More than one file found with the ending '{file_ending}'."
        else:
            return f"Error: No file found with the ending '{file_ending}'."

    except FileNotFoundError:
        return "Error: The directory does not exist."
    except pd.errors.EmptyDataError:
        return "Error: The file is empty."
    except pd.errors.ParserError:
        return "Error: The file could not be parsed."
    except Exception as e:
        return f"An unexpected error occurred: {e}"


def find_file_by_prefix(directory: str, file_prefix: str):
    """
    Finds a file in the specified directory that starts with the given prefix.
    :param directory: Directory to search for the file.
    :param file_prefix: Prefix of the file name to search for.
    :return: Path to the file if found, otherwise an error message.
    """
    try:
        # Ensure the directory exists
        if not os.path.isdir(directory):
            return "Error: The specified directory does not exist."

        # Create a pattern to match files starting with the given prefix
        pattern = os.path.join(directory, f'{file_prefix}*')

        # Find all files matching the pattern
        file_list = glob.glob(pattern)
        if len(file_list) == 1:
            return file_list[0]  # Return the path of the single matching file
        elif len(file_list) > 1:
            return "Error: More than one file found with the specified prefix."
        else:
            return "Error: No files found with the specified prefix."

    except Exception as e:
        return f"Error: An unexpected error occurred: {e}"


def parse_pileup_file(path_to_pileup_file: str):

    df_pileup = pd.read_csv(path_to_pileup_file, sep="\t")

    return df_pileup

def mean_coverage(df: pd.DataFrame) -> float:
    """
    calculates the mean coverage of the pileup file
    :param df: a dataframe with SNP coverage data
    :return: the mean coverage
    """
    return mean(df["cov"])

def uniformity_percentage(df: pd.DataFrame) -> float:
    """
    calculates the uniformity as the percentage of SNPs that have a coverage below 20% of the overall mean coverage.
    :param df: a dataframe with SNP coverage data
    :return: the uniformity percentage
    """
    overall_mean_coverage = df["cov"].mean()
    threshold = 0.2 * overall_mean_coverage
    above_threshold_count = df[df["cov"] > threshold].shape[0]
    total_snps = df.shape[0]
    uniformity = (above_threshold_count / total_snps) * 100
    return uniformity

def quality_metrics(df: pd.DataFrame) -> dict:
    """
    calculates the quality metrics of a pileup file
    :param df: a dataframe
    :return: a dictionary with the quality metrics
    """
    mean_cov = mean_coverage(df)
    uniformity = uniformity_percentage(df)
    metrics = {"mean_coverage": mean_cov, "uniformity": uniformity}
    return metrics

def outlier_detection(df: pd.DataFrame, zscore: int) -> tuple:
    """
    detects outliers in the dataframe based on the zscore threshold.
    :param df: dataframe
    :return: dataframe with the outliers removed
    """
    # Calculate the mean and standard deviation of the dd-cfDNA% values
    mean_val = df["dd-cfDNA%"].mean()
    std_dev = df["dd-cfDNA%"].std()

    # Remove rows where the dd-cfDNA% value is more than "zscore" standard deviations away from the mean
    df_wo_outliers = df[df["dd-cfDNA%"] - mean_val < zscore * std_dev]

    num_outliers = len(df) - len(df_wo_outliers)

    return df_wo_outliers, num_outliers

def determine_recipient_homozygous_base(row: pd.Series) -> str:
    """
    determines the homozygous base in the recipient
    :param row:
    :return:
    """
    if row["A"] > row["T"]:
        base = "A"
    elif row["A"] < row["T"]:
        base = "T"
    else:
        base = "Error"
    return base


def calculate_ddcfdna(row: pd.Series, genomic=False) -> pd.Series:
    """
    calculates the fraction of dd-cfDNA with the lower frequency allele as the donor allele
    :param row: a row of a dataframe
    :param genomic: boolean, if True, the recipient genotype is used to determine the dd-cfDNA fraction
    :return: a row of a dataframe with the dd-cfDNA fraction
    """

    if genomic:
        homozygous = row["HomozygousBaseRecipient"].upper()
        if homozygous in ["A", "T"]:
            donor_allele = "T" if homozygous == "A" else "A"
            row["dd-cfDNA%"] = (row[donor_allele] / row["Reads"]) * 100
        else:
            row["dd-cfDNA%"] = "Error: The recipient is not homozygous at this position for A or T."
    else:
        donor_allele = "A" if row["T"] > row["A"] else "T"
        row["dd-cfDNA%"] = (row[donor_allele] / row["Reads"]) * 100
    return row


def filter_rows(df):
    """
    removes the top row if it is more than 3 times the average of the next three rows
    :param df: dataframe
    :return: filtered dataframe
    """

    nrows_remove = 0
    while True:
        # Sort the dataframe in descending order of dd-cfDNA%
        df = df.sort_values(by=["dd-cfDNA%"], ascending=False).reset_index(drop=True)

        # Check if there are at least 3 rows to compare
        if len(df) < 3:
            break  # Not enough data to apply the rule

        # Calculate the average of the next two highest values
        average_of_next_two = df.loc[1:5, "dd-cfDNA%"].mean()

        # Check if the highest value is more than twice the average of the next two
        if df.loc[0, "dd-cfDNA%"] > 1.4 * average_of_next_two:
            # Remove the top row
            df = df.iloc[1:].reset_index(drop=True)
            nrows_remove += 1
        else:
            break  # Condition is not met, exit the loop
    return df, nrows_remove


def gmm_fitting(df: pd.DataFrame) -> tuple:
    """
    performs the GMM fitting to determine the corrected dd-cfDNA fraction
    :param df: dataframe with the %dd-cfDNA fraction
    :return: results of the GMM fitting
    """
    ddcfdna_frac = df["dd-cfDNA%"].values.reshape(-1, 1)

    try:
        gmm = GaussianMixture(n_components=3, max_iter=200)
        unadjusted = mean(df["dd-cfDNA%"]) * 2
        gmm.fit(ddcfdna_frac)
        labels = gmm.predict(ddcfdna_frac)
        df_sum = pd.DataFrame(gmm.means_, columns=["means"])
        df_sum["weights"] = gmm.weights_
        df_sum = df_sum.sort_values(by=["means"], ignore_index=True)

        average_frac = ((df_sum.loc[1]["means"] * 2 * df_sum.loc[1]["weights"]
                         + df_sum.loc[2]["means"] * df_sum.loc[2]["weights"])
                        / (df_sum.loc[1]["weights"] + df_sum.loc[2]["weights"]))

        noise = df_sum.loc[0]["means"]
        noise_adj = average_frac - noise
        no_snp = len(df["dd-cfDNA%"])

        return unadjusted, average_frac, noise_adj, no_snp, labels, ddcfdna_frac, noise, ""

    except:
        if len(df) == 0:
            return "NA", "NA", "NA", "NA", "NA", "NA", "NA", "too low coverage"
        else:
            return "NA", "NA", "NA", "NA", "NA", "NA", "NA", "other"


def plot_results(df: pd.DataFrame, labels, ddcfdna_frac, plot_output_dir: str, filename: str):
    """
    plot the results of the %dd-cfDNA calculation of the respective sample
    :param df: dataframe with the %dd-cfDNA fraction and the donor allele count
    :param labels:
    :param ddcfdna_frac:
    :return: nothing
    """

    plt.figure(1)
    ax = plt.subplot(121)
    plt.plot(sorted(df["dd-cfDNA%"]), linestyle='-', marker='.', markersize=3)
    plt.ylabel('Fraction', fontsize=6)
    ax.tick_params(labelsize=5)
    plt.title('Donor Allele Fraction', fontsize=8)

    fractions = pd.DataFrame({'dd-cfDNA%': list(ddcfdna_frac.reshape(1, -1)[0]), 'cluster': list(labels)}).sort_values(
        by=["dd-cfDNA%"], ascending=True)
    fractions.reset_index(drop=True, inplace=True)

    ax = plt.subplot(122)
    plt.plot(fractions[fractions["cluster"] == 0]["dd-cfDNA%"], linestyle='-', marker='.', markersize=3, color="r")
    plt.plot(fractions[fractions["cluster"] == 1]["dd-cfDNA%"], linestyle='-', marker='.', markersize=3, color="g")
    plt.plot(fractions[fractions["cluster"] == 2]["dd-cfDNA%"], linestyle='-', marker='.', markersize=3, color="b")

    ax.tick_params(labelsize=5)
    plt.title('dd-cfDNA fraction cluster', fontsize=8)
    plt.suptitle(filename, fontsize=10, fontweight="bold")
    plt.savefig(plot_output_dir + filename + "_dd-cfDNAfraction.png", dpi=300)
    plt.clf()


def determine_ddcfdna(df_cfdna, min_cov_ddcfdna, homozygous_threshold, zscore):
    """
    determines the dd-cfDNA fraction of a sample
    :param df_cfdna: dataframe with the pileup results
    :param min_cov: minimum coverage of a SNP required to consider the fraction for the final calculation
    :param homozygous_threshold: threshold to filter heterozygous SNPs out
    :param remove_top_percent: top x% of SNPs to be removed
    :return: a lot of different information...
    """

    # calculate quality metrics
    metrics = quality_metrics(df_cfdna)

    # calculate dd-cfDNA fraction for each SNP
    df_cfdna = df_cfdna.apply(calculate_ddcfdna, axis=1)

    # filter by coverage
    df_cfdna_filt = df_cfdna[df_cfdna["Reads"] >= min_cov_ddcfdna]

    # filter heterozygous SNPs out
    df_cfdna_filt = df_cfdna_filt[df_cfdna_filt["dd-cfDNA%"] < homozygous_threshold]

    # remove top x% of SNPs.
    df_cfdna_filt, nrows_remove = outlier_detection(df_cfdna_filt, zscore)

    # perform GMM fitting to determine the corrected dd-cfDNA fraction
    unadjusted, average_frac, noise_adj, no_snp, labels, ddcfdna_frac, noise, error = gmm_fitting(df_cfdna_filt)

    return unadjusted, average_frac, noise_adj, noise, no_snp, error, labels, \
        ddcfdna_frac, df_cfdna_filt, nrows_remove, metrics


def determine_ddcfdna_with_recipient_info(df_cfdna, df_genomic_recipient, min_cov_ddcfdna, min_cov_recip,
                                          homozygous_threshold_recip, zscore):
    """
    determines the dd-cfDNA fraction of a sample using the genomic information of the recipient
    :param df_cfdna: dataframe with the pileup results
    :param df_genomic_recipient: dataframe with the pileup results of the recipient genomic DNA or equivalent
    :param min_cov_ddcfdna: minimum coverage of a SNP required to consider the fraction for the final calculation
    :param min_cov_recip: minimum coverage of a SNP for homozygous recipien SNP identification
    :param homozygous_threshold: threshold to filter heterozygous SNPs out
    :param remove_top_percent: top x% of SNPs to be removed
    :return: a lot of different information...
    """

    # calculate quality metrics
    metrics = quality_metrics(df_cfdna)

    # calculate recipient genotype for each SNP
    df_genomic_recipient = df_genomic_recipient.apply(calculate_ddcfdna, axis=1)

    # remove heterozygous SNPs in recipient
    df_genomic_recipient_homo = df_genomic_recipient[df_genomic_recipient["dd-cfDNA%"] < homozygous_threshold_recip]

    # filter by coverage
    df_genomic_recipient_homo_filt = df_genomic_recipient_homo[df_genomic_recipient_homo["Reads"] >= min_cov_recip]

    # determine homozygous base in recipient
    df_genomic_recipient_homo_filt["HomozygousBaseRecipient"] = \
        df_genomic_recipient_homo_filt.apply(determine_recipient_homozygous_base, axis=1)

    # keep relevant columns only
    df_genomic_recipient_homo_filt = \
        df_genomic_recipient_homo_filt[["Chr", "Position", "ref", "HomozygousBaseRecipient"]]

    # filter by coverage
    df_cfdna_filt = df_cfdna[df_cfdna["Reads"] >= min_cov_ddcfdna]

    # remove SNPs from donor df where recipient is heterozygous
    df_cfdna_filt = df_cfdna_filt.merge(df_genomic_recipient_homo_filt, how='inner', on=["Chr", "Position"])

    # calculate dd-cfDNA fraction for each SNP
    df_cfdna_filt = df_cfdna_filt.apply(calculate_ddcfdna, args=(True,), axis=1)

    # outlier detection
    df_cfdna_filt, nrows_remove = outlier_detection(df_cfdna_filt, zscore)

    # remove top SNP if 3 times higher than average of previous three SNPs.
    #df_cfdna_filt, nrows_remove = filter_rows(df_cfdna_filt)

    # print(df_cfdna_filt[["Chr", "Position", "dd-cfDNA%"]])
    # perform GMM fitting to determine the corrected dd-cfDNA fraction
    unadjusted, average_frac, noise_adj, no_snp, labels, ddcfdna_frac, noise, error = gmm_fitting(df_cfdna_filt)

    return unadjusted, average_frac, noise_adj, noise, no_snp, error, labels, \
        ddcfdna_frac, df_cfdna_filt, nrows_remove, metrics


def check_directory(raw_path: str):
    """
    checks if the input directory is a valid directory
    :param raw_path: the raw path
    :return: the path with the correct format or an error message
    """
    if not os.path.isdir(raw_path):
        raise argparse.ArgumentTypeError('"{}" is not an existing directory'.format(raw_path))
    elif raw_path[-1] != "/":
        raw_path += "/"
        return raw_path
    elif raw_path[-1] == "/":
        return raw_path
    else:
        raise argparse.ArgumentTypeError('"{}" is not a valid directory'.format(raw_path))


def check_output_file(filename: str):
    """
    checks if the output file is a valid filename and not a directory.
    :param filename: a filename
    :return: the filename or error message
    """
    if os.path.isdir(filename):
        raise argparse.ArgumentTypeError('"{}" is an existing directory. Enter a filename without extention.'
                                         .format(filename))
    elif filename.count(".") > 0:
        raise argparse.ArgumentTypeError('"{}" is not a valid filename. Enter a filename without extention or dots.'
                                         .format(filename))
    return filename


def create_output_dir(args):
    """
    creates the output directory
    :param args: the parsed arguments
    :return: output directory or error message
    """
    if args.output_dir is None:
        output_dir = args.input_dir + args.output_file + ".csv"
        print("The output file will be saved to:", args.input_dir + args.output_file + ".csv")
        return output_dir

    elif os.path.isdir(args.output_dir):
        output_dir = args.output_dir + args.output_file + ".csv"
        print("The output file will be saved to:", args.output_dir + args.output_file + ".csv")
        return output_dir

    else:
        raise argparse.ArgumentTypeError('"{}" is not an existing directory'.format(args.output_dir))


def check_plot_output_dir(args):
    """
    checks if the plot output directory is specified if the plot flag is added
    :param args: the parsed arguments
    :return: error message if the plot output directory is not specified
    """
    if args.plot_results and not args.plot_output_dir:
        raise argparse.ArgumentTypeError("You added the -plot flag, but did not specify the directory where the plots")


def init_arg_parser():
    parser = argparse.ArgumentParser()
    return parser


def take_SNPs_subset(df_cfdna_file, df_subset_file):
    """
    takes the subset of SNPs from the cfDNA file that are in the subset file
    :param df_cfDNA_file: dataframe with the pileup results of the cfDNA
    :param df_subset_file: dataframe with the subset of SNPs
    :return: the subset of SNPs
    """
    # version 3 panel
    df_subset = df_cfdna_file.merge(df_subset_file, how='inner', on=["Chr", "Position"])

    return df_subset

def exclude_SNPs_subset(df_cfdna_file, df_subset_file):
    """
    takes the subset of SNPs from the cfDNA file that are in the subset file
    :param df_cfDNA_file: dataframe with the pileup results of the cfDNA
    :param df_subset_file: dataframe with the subset of SNPs
    :return: the subset of SNPs
    """

    # remove susbitutous SNPs
    matching_rows = df_cfdna_file.merge(df_subset_file, on=["Chr", "Position"], how="inner")

    # Exclude matching rows
    df_subset = df_cfdna_file[~df_cfdna_file.set_index(["Chr", "Position"]).index.isin(
        matching_rows.set_index(["Chr", "Position"]).index)]

    return df_subset


def add_arguments(parser):
    parser.add_argument("-dir", "--input_dir", type=check_directory, required=True,
                        help="The directory containing the pileup files.")
    parser.add_argument("-o", "--output_file", type=check_output_file, required=True,
                        help="The name of the output file.")
    parser.add_argument("-e_in", "--pileup_file_ending", required=False, default="_cons.pileup.tsv",
                        type=str, help="The name ending of the pileup file. Default='%(default)s'")
    parser.add_argument("-o_dir", "--output_dir", type=check_directory, required=False,
                        help="The path where the output file should be saved. If not specified, the file will be saved"
                             "to the input directory containing the pileup files.")
    parser.add_argument("-c", "--min_cov", required=False, type=int, default=300,
                        help="The minimum coverage required for a SNP to be considered for the final calculation.")
    parser.add_argument("-c_recip", "--min_cov_genotyping", required=False, type=int, default=50,
                        help="The minimum coverage required for a SNP to be considered for the recipient genotyping.")
    parser.add_argument("-z", "--outlier_zscore", required=False, type=float, default=3,
                        help="The z-score threshold to remove outliers.")
    parser.add_argument("-gmap", "--genomic_mapping_file_ending", required=False, default="_genomic_mapping.csv",
                        type=str, help="The name of the genomic mapping file. Default='%(default)s'")
    parser.add_argument("-SNPthres", "--homozygous_threshold", required=False, default=30,
                        type=float, help="The frequency threshold to filter heterozygous SNPs out. Default=%(default)s")
    parser.add_argument("-SNPthres_recip", "--homozygous_threshold_recip", required=False, default=10,
                        type=float, help="The frequency threshold to filter heterozygous SNPs out in the recipient genotyping."
                                         "Default=%(default)s")
    parser.add_argument("-plot", "--plot_results", required=False, default=False, action='store_true',
                        help="Add this flag if the results should be plotted. If this flag is added, the directory"
                             "where the plots should be save has to be saved using the --plot_output_dir flag.")
    parser.add_argument("-plotout", "--plot_output_dir", type=check_directory, required=False,
                        help="The directory where the plots should be saved.")
    parser.add_argument("-info", "--save_sample_results", required=False, default=False, action='store_true',
                        help="Add this flag if the results for each sample should be saved in a separate csv file.")
    parser.add_argument("-top", "--remove_top_percent", required=False, default=0.05, type=float,
                        help="The top x of SNPs will be removed from the calculation. "
                             "This is to avoid outliers. Default=%(default)s")
    return parser


def main():
    """
    main function executing the whole workflow for the calculation of the dd-cfDNA fraction
    saves the csv file with the results to the specified output directory and takes the arguments from the command line
    :return: nothing
    """

    # Read the subset of SNPs to be removed that have shown to be problematic in the past
    df_subset_remove = pd.read_csv(os.path.dirname(__file__) + "/exclude_SNPs_chr_position.csv")

    # Initialize the parser
    parser = init_arg_parser()

    # Add arguments
    add_arguments(parser)

    # Parse arguments
    args = parser.parse_args()

    # Check if the plot output directory is specified if the plot flag is added
    check_plot_output_dir(args)

    # Create the output directory
    output_path_file = create_output_dir(args)

    df_results = pd.DataFrame()
    genomic_mapping_ending = "*" + args.genomic_mapping_file_ending

    genomic_mapping_df = parse_genomic_mapping_file(args.input_dir, genomic_mapping_ending)

    # Check if the genomic mapping DataFrame is a string indicating an error
    if type(genomic_mapping_df) == str and genomic_mapping_df.startswith("Error: No file found with the ending"):
        genomic_mapping_avail = False
    else:
        genomic_mapping_avail = True

    sample_counter = 0  # Counter for the number of samples processed
    last_percent_reported = 0  # Track the last reported percentage
    total_samples = len([filename for filename in os.listdir(args.input_dir) if
                         filename.endswith(args.pileup_file_ending)])  # Total number of samples
    print("")
    print(f"Total number of samples to process: {total_samples}")
    print("Starting the dd-cfDNA fraction calculation...")
    # Iterate through each pileup file in the input directory
    for filename in os.listdir(args.input_dir):
        if filename.endswith(args.pileup_file_ending):  # and filename.split("_cons", 1)[0] == "OLT29_2p_S2": #in genomic_mapping_df["cfDNA_file_name"].values:
            sample_name = filename.split(".", 1)[0]

            # Check if the sample has genomic information available
            if genomic_mapping_avail and sample_name in genomic_mapping_df["cfDNA_file_name"].values:
                genomic_sample_name = genomic_mapping_df[
                    genomic_mapping_df["cfDNA_file_name"] == sample_name]["genomicDNA_file_name"].iloc[0]

                genomic_file_name = genomic_sample_name
                genomicfile_path_or_error = find_file_by_prefix(args.input_dir, genomic_file_name)

                if genomicfile_path_or_error.startswith("Error"):  # Log error message in the final output table
                    df_results_sample = pd.DataFrame({"SampleID": [sample_name],
                                                      "UnadjustedFraction": "",
                                                      "AverageFraction": "",
                                                      "NoiseAdjustedFraction": "",
                                                      "SNPsUsedCount": "",
                                                      "SNPsRemovedCount": "",
                                                      "MeanCoverageAllSNPs": "",
                                                      "Uniformity": "",
                                                      "Noise": "",
                                                      "with_genomic_info": "",
                                                      "Error": genomicfile_path_or_error})
                    df_results = pd.concat([df_results, df_results_sample], axis=0)
                    continue  # Continue with the next file

                file_path = os.path.join(args.input_dir, filename)
                df_cfdna = parse_pileup_file(file_path)
                df_genomic = parse_pileup_file(genomicfile_path_or_error)

                # Remove SNPs that are in the subset of SNPs to be removed
                df_cfdna = exclude_SNPs_subset(df_cfdna, df_subset_remove)

                (unadjusted, average_frac, noise_adj, noise, no_snp, error, labels, ddcfdna_frac,
                 df_fraction, nrows_remove, metrics) = \
                    determine_ddcfdna_with_recipient_info(df_cfdna=df_cfdna,
                                                          df_genomic_recipient=df_genomic,
                                                          min_cov_ddcfdna=args.min_cov,
                                                          min_cov_recip=args.min_cov_genotyping,
                                                          homozygous_threshold_recip=args.homozygous_threshold_recip,
                                                          zscore=args.outlier_zscore)

                if error == "" and args.plot_results:
                    plot_results(df_fraction, labels, ddcfdna_frac, args.plot_output_dir, sample_name)

                df_results_sample = pd.DataFrame({"SampleID": [sample_name],
                                                  "UnadjustedFraction": [unadjusted],
                                                  "AverageFraction": [average_frac],
                                                  "NoiseAdjustedFraction": [noise_adj],
                                                  "SNPsUsedCount": [no_snp],
                                                  "SNPsRemovedCount": [nrows_remove],
                                                  "MeanCoverageAllSNPs": [metrics["mean_coverage"]],
                                                  "Uniformity": [metrics["uniformity"]],
                                                  "Noise": [noise],
                                                  "with_genomic_info": "yes",
                                                  "Error": [error]})
                df_results = pd.concat([df_results, df_results_sample], axis=0)
                if args.save_sample_results:
                    df_fraction.to_csv(args.output_dir + sample_name + "_dd-cfDNAfraction.csv", index=False)

            else:
                # If genomic information is not available, process the sample without it
                file_path = os.path.join(args.input_dir, filename)
                df_cfdna = parse_pileup_file(file_path)

                df_cfdna = exclude_SNPs_subset(df_cfdna, df_subset_remove)

                (unadjusted, average_frac, noise_adj, noise, no_snp, error, labels, ddcfdna_frac,
                 df_fraction, nrows_remove, metrics) = \
                    determine_ddcfdna(df_cfdna=df_cfdna,
                                      min_cov_ddcfdna=args.min_cov,
                                      homozygous_threshold=args.homozygous_threshold,
                                      zscore=args.outlier_zscore)

                if error == "" and args.plot_results:
                    plot_results(df_fraction, labels, ddcfdna_frac, args.plot_output_dir, sample_name)

                df_results_sample = pd.DataFrame({"SampleID": [sample_name],
                                                  "UnadjustedFraction": [unadjusted],
                                                  "AverageFraction": [average_frac],
                                                  "NoiseAdjustedFraction": [noise_adj],
                                                  "SNPsUsedCount": [no_snp],
                                                  "SNPsRemovedCount": [nrows_remove],
                                                  "MeanCoverageAllSNPs": [metrics["mean_coverage"]],
                                                  "Uniformity": [metrics["uniformity"]],
                                                  "Noise": [noise],
                                                  "with_genomic_info": "no",
                                                  "Error": [error]})
                df_results = pd.concat([df_results, df_results_sample], axis=0)
                if args.save_sample_results:
                    df_fraction.to_csv(args.output_dir + sample_name + "_dd-cfDNAfraction.csv", index=False)

            sample_counter += 1  # Increment the sample counter

            # Calculate the percentage of processed samples
            percent_processed = (sample_counter / total_samples) * 100

            # Print progress only if it is 5% more than the last reported percentage
            if percent_processed - last_percent_reported >= 5:
                print(f"Processed {percent_processed:.0f}% of samples")
                last_percent_reported = percent_processed

    df_results.to_csv(output_path_file, index=False)
    print(f"The script finished successfully. The results are saved to {output_path_file}")


if __name__ == "__main__":
    main()