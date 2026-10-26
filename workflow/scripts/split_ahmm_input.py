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
import pathlib


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", type=str, required=True,
        help="Input Ancestry HMM-format panel file")
    parser.add_argument(
        "-s", "--samples", type=str, required=True,
        help="Sample file used to generate `input`")
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

    indices = dict()
    outputs = dict()
    for sample_idx, sample in enumerate(samples):
        indices[sample_idx] = [2 * sample_idx + 7, 2 * sample_idx + 8]
        outputs[sample_idx] = f"{args.output}{sample}.tsv"
        sample_file = f"{args.output}{sample}.samples"
        with open(sample_file, "w") as fout:
            fout.write(f"{sample}\t2\n")

    for sample_idx in outputs:
        file = outputs[sample_idx]
        if pathlib.Path(file).exists():
            raise ValueError(f"output file {file} already exists!")

    with open(args.input, "r") as fin:
        for line in csv.reader(fin, delimiter="\t"):
            constant = "\t".join(line[:7])
            for sample_idx in indices:
                get_idx = indices[sample_idx]
                line_out = (
                    constant + "\t" + "\t".join([line[i] for i in get_idx]) + "\n")
                with open(outputs[sample_idx], "a") as fout:
                    fout.write(line_out)

