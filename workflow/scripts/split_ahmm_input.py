"""
Extract each sample from an Ancestry HMM format panel file into its own file,
retaining the first 7 rows, which contain site position information and source
population panel allele counts.

Columns (unlabelled) are:
    chrom pos ref0A ref0a ref1A ref1a distance(M) sample0A sample0a ...
Where A is the count of reference alleles and a is the alternate allele count.
"""

import argparse
import csv


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", type=str, required=True,
        help="Input Ancestry HMM-format panel file")
    parser.add_argument(
        "-s", "--samples", type=str, required=True,
        help="Sample file used to generate `input`")
    parser.add_argument("--sample", required=True,)
    parser.add_argument("-o", "--output", required=True)
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    # Read sample files
    samples = []
    with open(args.samples, "r") as fin:
        for line in csv.reader(fin, delimiter="\t"):
            if line[1] == "admixed":
                samples.append(line[0])

    sample_idx = samples.index(args.sample)
    get_idx = list(range(7)) + [2 * sample_idx + 7, 2 * sample_idx + 8]

    with open(args.input, "r") as fin:
        for line in csv.reader(fin, delimiter="\t"):
            line_out = "\t".join([line[i] for i in get_idx])
            print(line_out)

    with open(args.output, "w") as fout:
        fout.write(f"{args.sample}\t2\n")

