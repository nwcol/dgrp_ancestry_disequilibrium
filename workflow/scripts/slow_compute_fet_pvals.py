"""
Read a table with ancestry assignments and compute Fisher exact test (FET)
p-values for this comparison:
    00  01
    10  11
Where the first and second characters are ancestry at two loci and 0, 1 are
ancestry assignments.
"""

import argparse
import csv
import numpy as np
import scipy


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True,
        help="")
    parser.add_argument("-arm1", "--arm1", required=True,)
    parser.add_argument("-arm2", "--arm2", required=True,)
    parser.add_argument("-b", "--bin-output", required=True,
        help="Output file (.tsv) for binned FET p-values")
    parser.add_argument("-s", "--sig-output", required=True,
        help="Output file (.tsv) for sites with threshold-passing p-values")
    parser.add_argument("--threshold", type=float, default=1e-3)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    # Set up p-value bins
    n_bins = 1000
    p_bins = np.zeros(n_bins + 1, dtype=np.int64)
    # Load and cache input data
    data1 = dict()
    data2 = dict()
    with open(args.input, "r") as fin:
        for line in csv.reader(fin, delimiter="\t"):
            # Skip header line
            if line[0] == "chrom":
                continue
            arm, pos = line[:2]
            if arm == args.arm1:
                to_data = data1
            elif arm == args.arm2:
                to_data = data2
            else:
                continue
            anc = np.array([int(a) for a in line[3:]])
            to_data[(arm, pos)] = anc
    # Scan locus pairs
    with open(args.sig_output, "w") as fout:
        for site1, anc1 in data1.items():
            for site2, anc2 in data2.items():
                n00 = np.sum((anc1 == 0) & (anc2 == 0))
                n01 = np.sum((anc1 == 0) & (anc2 == 1))
                n10 = np.sum((anc1 == 1) & (anc2 == 0))
                n11 = np.sum((anc1 == 1) & (anc2 == 1))
                tab = np.array([[n00, n01], [n10, n11]])
                p_fet = scipy.stats.fisher_exact(tab, alternative="greater").pvalue
                p_idx = int(p_fet * n_bins)
                p_bins[p_idx] += 1
                # Print site information if p is significant enough
                # Info: arm1 pos1 arm2 pos2 n00 n01 n10 n11 p_fet
                if p_fet < args.threshold:
                    info = [site1[0], site1[1], site2[0], site2[1],
                            n00, n01, n10, n11, p_fet]
                    line_out = "\t".join([str(x) for x in info]) + "\n"
                    fout.write(line_out)
    # Output p-value bins and statistics
    with open(args.bin_output, "w") as fout:
        for ii, count in enumerate(p_bins):
            bin_start = ii / n_bins
            line_out = f"{bin_start}\t{count}\n"
            fout.write(line_out)

