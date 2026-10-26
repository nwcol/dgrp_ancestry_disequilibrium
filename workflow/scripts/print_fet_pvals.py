"""
Read a table with ancestry assignments and compute Fisher exact test (FET)
p-values for this comparison:
    00  01
    10  11
Where the first and second characters are ancestry at two loci and 0, 1 are
ancestry assignments.
"""

import argparse
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "input", required=True,
        help="")
    parser.add_argument("-arm1", "--arm1", required=True,)
    parser.add_argument("-arm2", "--arm2", required=True,)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    # Set up p-value bins

    # Load and cache input data


    # Scan locus pairs
    for ii in range(len(arm1_data)):
        for jj in range(len(arm2_data)):
            n00 = 0
            n01 = 0
            n10 = 0
            n11 = 0
            for kk in range(n_samples):
                
                if data[ii][kk] == -1 or data[jj][kk] == -1:
                    continue
               
                if data[ii][kk] == 0:
                    if data[jj][kk] == 0:
                        n00 += 1
                    else:
                        n01 += 1
                else:
                    if data[jj][kk] == 0:
                        n10 += 1
                    else:
                        n11 += 1

            tab = np.array([[n00, n01], [n10, n11]])
    # Output format: locus1 locus2 n00 n01 n10 n11 pFET
    

    # Output p-value bins and statistics


