"""
Merge several VCF-derived, tab-separated files with these columns (and no
header line):
    <chrom> <position> <ref> <alt> <ac>
Into a file with this format:
    <chrom> <position> <ref> <alt0> <ac0> <alt1> <ac1> ...
"""

import argparse
import collections
import csv
import gzip


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", nargs="*", required=True)
    parser.add_argument("--labels", nargs="*", required=True)
    parser.add_argument("-o", "--output", required=True)
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    assert len(args.input) == len(args.labels)
    cols = ["CHROM", "POS", "REF"]
    for label in args.labels:
       cols += [f"{x}_{label}" for x in ["ALT", "AC"]]
    header = "\t".join(cols) + "\n"
    # Read data
    make_func = lambda: {l: [".", "."] for l in args.labels}
    raw_data = collections.defaultdict(make_func)
    for label, in_file in zip(args.labels, args.input):
        with open(in_file, "r") as fin:
            for line in csv.reader(fin, delimiter="\t"):
                chrom, pos, ref, alt, ac = line
                raw_data[(chrom, pos, ref)][label] = [alt, ac]
    # Write data
    with open(args.output, "w") as fout:
        fout.write(header)
        for key in raw_data:
            chrom, pos, ref = key
            info = [chrom, pos, ref]
            for label in args.labels:
                info += raw_data[key][label]
            line = "\t".join(info) + "\n"
            fout.write(line)

