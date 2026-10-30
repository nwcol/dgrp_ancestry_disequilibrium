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
    parser.add_argument("-i1", "--input1", required=True,)
    parser.add_argument("-i2", "--input2", required=True,)
    parser.add_argument("-b", "--bin-output", required=True,
        help="Output file (.tsv) for binned FET p-values")
    parser.add_argument("-s", "--sig-output", required=True,
        help="Output file (.tsv) for sites with threshold-passing p-values")
    parser.add_argument("--threshold", type=float, default=1e-3)
    return parser.parse_args()


def load_input(filename):
    data = dict()
    with open(filename, "r") as fin:
        for line in csv.reader(fin, delimiter="\t"):
            # Skip header line
            if line[0] == "chrom":
                continue
            arm, pos = line[:2]
            anc = np.array([int(a) if a != "0.5" else -1 for a in line[3:]],
                           dtype=np.int16)
            data[(arm, pos)] = anc
    return data


if __name__ == "__main__":
    args = parse_args()
    # Set up p-value bins
    n_bins = 1000
    p_bins = np.linspace(0, 1, n_bins + 1)
    binned_p = np.zeros(n_bins, dtype=np.int64)
    # Load and cache input data
    data1 = load_input(args.input1)
    data2 = load_input(args.input2)
    # Scan locus pairs
    sites2 = list(data2.keys())
    arm2_arr = np.array([data2[s] for s in sites2], dtype=np.int16)
    tab_arr = np.zeros((len(arm2_arr), 2, 2), dtype=np.int64)
    with open(args.sig_output, "w") as fout:
        for site1, anc1 in data1.items():
            tab_arr[:, 0, 0] = np.sum((anc1 == 0) & (arm2_arr == 0), axis=1)
            tab_arr[:, 0, 1] = np.sum((anc1 == 0) & (arm2_arr == 1), axis=1)
            tab_arr[:, 1, 0] = np.sum((anc1 == 1) & (arm2_arr == 0), axis=1)
            tab_arr[:, 1, 1] = np.sum((anc1 == 1) & (arm2_arr == 1), axis=1)
            p_fets = np.array(
                [scipy.stats.fisher_exact(tab, alternative="greater").pvalue
                for tab in tab_arr]
            )
            binned_p += np.histogram(p_fets, bins=p_bins)[0]
            for ii in np.where(p_fets <= args.threshold)[0]:
                info = [site1[0], site1[1], sites2[ii][0], sites2[ii][1],
                        p_fets[ii], tab_arr[ii, 0, 0], tab_arr[ii, 0, 1],
                        tab_arr[ii, 1, 0], tab_arr[ii, 1, 1]]
                line_out = "\t".join([str(x) for x in info]) + "\n"
                fout.write(line_out)
    # Output p-value bins and statistics
    with open(args.bin_output, "w") as fout:
        for ii, count in enumerate(binned_p):
            bin_start = p_bins[ii]
            line_out = f"{bin_start}\t{count}\n"
            fout.write(line_out)
