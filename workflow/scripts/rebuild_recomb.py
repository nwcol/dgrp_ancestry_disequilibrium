"""
Modify a tab-separated recombination map by:
- Extending the map (with rate 0.01 cM/Mb, or r=1e-10) to chromosome start/end
- Renaming 'Arm' column to 'Chrom'
- Relabeling chromosomes using a tab-separated mapping file (--labels)
- Renaming 'c (cM/Mb)' column to 'cM/Mb'.
"""

import pandas
import argparse


# Recombination rate for missing regions
C_MISS = 0.01


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-l", "--labels", required=True)
    parser.add_argument("-g", "--genome", required=True,
                        help="Genome file with tab-separated chrom label, length")
    parser.add_argument("-o", "--output", required=True)
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
    # Extend the map to chromosome start/end positions
    genome = pandas.read_csv(
        args.genome, sep="\\s+", names=["Chrom", "Length"], index_col="Chrom")
    extended_sub_dfs = []
    for chrom in set(df["Chrom"]):
        sub_df = df[df["Chrom"] == chrom]
        map_start = sub_df.iloc[0]["Start"]
        map_end = sub_df.iloc[-1]["Stop"]
        chrom_start = 1
        chrom_end = genome.loc[chrom]["Length"]
        data = []
        if map_start > chrom_start:
            first_row = pandas.DataFrame(
                [[chrom, chrom_start, map_start - 1, C_MISS]],
                columns=col_names
            )
            data.append(first_row)
        data.append(sub_df)
        if chrom_end > map_end:
            last_row = pandas.DataFrame(
                [[chrom, map_end + 1, chrom_end, C_MISS]], columns=col_names)
            data.append(last_row)
        extended_sub_dfs.append(pandas.concat(data))
    out_df = pandas.concat(extended_sub_dfs)
    out_df.to_csv(args.output, sep="\t", index=0)


