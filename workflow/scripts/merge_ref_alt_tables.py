"""
"""

import argparse
import csv


def get_args();
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", nargs="*", required=True)
    parser.add_argument("--labels", nargs="*", required=True)
    parser.add_argument("-o", "--output", required=True)
    return parser.parse_args()


if __name__ == "__main__":
    

    cols = ["chrom", "pos"]
    rep_cols = ["ref", "alt", "Nalt"]
    for label in args.labels:
       cols += [f"{x}_{label}" for x in rep_cols]

 



