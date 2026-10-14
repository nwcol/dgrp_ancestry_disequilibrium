"""
Modify a tab-separated recombination map by:
- Renaming 'Arm' column to 'Chrom'
- Relabeling chromosomes using a tab-separated mapping file (--labels)
- Renaming 'c (cM/Mb)' column to 'cM/Mb'.
"""

import pandas
import argparse


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-l", "--labels", required=True)
    parser.add_argument("-o", "--output", requred=True)
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    col_names = ["Chrom", "Start", "Stop", "cM/Mb"]
    df = pandas.read_csv(args.input, sep="\\s+", names=col_names, skiprows=1)
    chr_map_df = pandas.read_csv(
        args.labels, sep="\\s+", names=["label", "code"])
    chr_map = {int(chr_map_df.loc[i, "code"]): chr_map_df.loc[i, "label"]
               for i in range(len(chr_map_df))}
    df.Chrom = [chr_map[x] for x in df["Chrom"]]
    df.to_csv(args.output, sep="\t", index=0)


