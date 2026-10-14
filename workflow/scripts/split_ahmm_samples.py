"""
Extract each sample from an Ancestry HMM format panel file into its own file,
retaining the first 7 rows, which contain site position information and source
population panel allele counts.

Columns expected in a panel file are not labeled. Columns are:
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
    parser.add_argument(
        "-o", "--output", type=str, required=True,
        help="Output prefix. Output filenames are {output}.{sample}.tsv")
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()

    # Scan input file
    contents = []
    with open(args.input, "r") as fin:
        for line in csv.reader(fin, delimiter="\t"):
            contents.append(line)

    # Read sample files
    samples = []
    with open(args.samples, "r") as fin:
        for line in csv.reader(fin, delimiter="\t"):
            if line[1] == "admixed":
                samples.append(line[0])

    for i, sample in enumerate(samples):
        output_file = f"{args.output}.{sample}.tsv"
        idx = list(range(7)) + [2 * i + 7, 2 * i + 8]
        with open(output_file, "w") as fout:
            for line in contents:
                line_out = "\t".join([line[j] for j in idx]) + "\n"
                fout.write(line_out)
        sample_file = f"{args.output}.{sample}.sample"
        with open(sample_file, "w") as fout:
            fout.write(f"{sample}\t2")




