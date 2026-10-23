"""
Read a VCF file and convert it to Ancestry HMM input format, then print to
stdout.

Adapted from vcf2ahmm.py on the Ancestry HMM github:
https://github.com/russcd/Ancestry_HMM/tree/master

Changes from original
---------------------
- Added support for gzipped VCF files
- Added support for recombination maps
- Added long forms for argument names
- Allow reading phased genotypes

Recombination map format
------------------------
Chrom\tStart\tStop\tc(cM/Mb)
2L\t1\t10000\t0.01
2L\t10001\t20000\t1.00
...

Positions (Start/Stop) are treated as 1-indexed and inclusive. Entries in
"Chrom" must match entries in the VCF "CHROM" column.
"""

import argparse
import collections
import csv
import gzip
import numpy as np
import pandas


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v", "--vcf",
        required=True,
        type=str,
        help="Input vcf file",
    )
    parser.add_argument(
        "-s", "--samples",
        required=True,
        type=str,
        help="Sample-to-population mapping file"
    )
    parser.add_argument(
        "-r", "--rec",
        required=True,
        type=str,
        help="Uniform recombination rate OR path to tab-separated recombination map file"
    )
    parser.add_argument(
        "-g",
        type=int,
        default=0,
        help="Boolean (0 = use sample reads| 1 = use genotypes)"
    )
    parser.add_argument(
        "-m",
        type=int,
        default=1000,
        help="Minimum distance between successive snps (int,bp)"
    )
    parser.add_argument(
        "-o",  "--out",
        type=str,
        default="ahmm.ploidy",
        help="Ploidy file for ahmm input"
    )
    parser.add_argument(
        "--min_samples",
        type=int,
        default=1,
        help="Minimum number of reference population samples (int)"
    )
    parser.add_argument(
        "--min_diff",
        type=float,
        default=0.0,
        help="minimum allele frequency difference between populations (float)"
    )
    return parser.parse_args()


def load_rec_map(filename, chrom):
    """
    Load a recombination map interpolation function.
    """
    df = pandas.read_csv(filename, sep="\\s+")
    # Subset to `chrom` and extract interval endpoints
    sub_df = df[df["Chrom"] == chrom]
    map_xs = [1]
    map_ys = [0]
    for _, row in sub_df.iterrows():
        start = row["Start"]
        stop = row["Stop"]
        # Convert to M/bp
        rate = row["cM/Mb"] * 1e-8
        map_span = (stop - start + 1) * rate
        map_y = map_ys[-1] + map_span
        map_xs.append(stop)
        map_ys.append(map_y)
    # Set up a function to interpolate site map coordinates
    rec_map = lambda x: np.interp(x, map_xs, map_ys)
    return rec_map


if __name__ == "__main__":
    args = get_args()
    # Get population-to-sample mappings
    pop_sample_map = collections.defaultdict(list)
    ordered_pops = []
    with open(args.samples) as file:
        for line in file:
            samp, pop = line.split()
            pop_sample_map[pop].append(samp)
            if pop not in ordered_pops and pop != "admixed":
                ordered_pops.append(pop)

    # Prepare recombination map(s)
    try:
        float(args.rec)
        use_uniform_r = True
    except:
        use_uniform_r = False

    if use_uniform_r:
        r = float(args.rec)

    # Keep track of chromosome, position
    last_position = 0
    last_map_coord = 0
    current_chrom = "."

    # Allow gzipped VCF
    if args.vcf.endswith(".gz"):
        open_func = lambda f: gzip.open(f, "rt")
    else:
        open_func = lambda f: open(f, "r")

    with open_func(args.vcf) as fin:
        for line in csv.reader(fin, delimiter="\t"):
            # Header information
            if line[0].startswith("##"):
                continue

            if line[0] == "#CHROM":
                header = line
                # Build a population-to-sample index mapping
                # These index the line, not the samples in the line. Hence
                # the index of the first sample column is 9.
                pop_idx_map = {
                    pop: [header.index(s) for s in pop_sample_map[pop]]
                    for pop in pop_sample_map}
                continue

            # Skip sites that are not biallelic SNPs
            if len(line[3]) != 1 or len(line[4]) != 1:
                continue

            if current_chrom != line[0]:
                current_chrom = line[0]
                last_position = 0
                last_map_coord = 0
                if not use_uniform_r:
                    rec_map = load_rec_map(args.rec, current_chrom)

            # Skip if the last position is too close
            elif int(line[1]) - last_position < args.m:
                continue

            # Get population allele counts
            alts = collections.defaultdict(int)
            tots = collections.defaultdict(int)
            for pop in ordered_pops:
                for idx in pop_idx_map[pop]:
                    # Assume genotype is the first three sample chars
                    genotype = line[idx][:3]
                    # Replace phased GT symbol
                    if "|" in genotype:
                        genotype = genotype.replace("|", "/")
                    if genotype == "0/0":
                        tots[pop] += 2
                    elif genotype == "0/1" or genotype == "1/0":
                        alts[pop] += 1
                        tots[pop] += 2
                    elif genotype == "1/1":
                        alts[pop] += 2
                        tots[pop] += 2
                    # Use haploid genotype if needed
                    elif genotype.startswith("0"):
                        tots[pop] += 1
                    elif genotype.startswith("1"):
                        alts[pop] += 1
                        tots[pop] += 1

            # Make sure there are enough sampled alleles
            skip = False
            freqs = []
            for pop, tot in tots.items():
                freqs.append(alts[pop]/tot)
                if tot < args.min_samples:
                    skip = True
            if skip:
                continue

            # Skip if allele frequency difference is not great enough
            for i in range(0, len(freqs) - 1):
                for j in range(i + 1, len(freqs)):
                    if abs(freqs[i] - freqs[j]) < args.min_diff:
                        skip = True
            if skip:
                continue

            # Start printing site information
            print(line[0], line[1], end="\t", sep="\t")
            for pop in ordered_pops:
                tot = tots[pop]
                alt = alts[pop]
                ref = tot - alt
                print(ref, alt, end="\t", sep="\t")

            # Compute map distance to the last site
            if use_uniform_r:
                map_coord = int(line[1]) * r
            else:
                position = int(line[1])
                map_coord = rec_map(position)
            map_dist = map_coord - last_map_coord
            print(map_dist, end="")

            # Print counts for each admixed sample
            for idx in pop_idx_map["admixed"]:
                # Use genotype information
                if args.g == 1:
                    # For diploids, GT should be the first three characters X/X
                    genotype = line[idx][:3]
                    # Replace phased GT symbol
                    if "|" in genotype:
                        genotype = genotype.replace("|", "/")
                    if genotype == "0/0":
                        print("\t", 2, "\t", 0, end="", sep="")
                    elif genotype == "0/1":
                        print("\t", 1, "\t", 1, end="", sep="")
                    elif genotype == "1/1":
                        print("\t", 0, "\t", 2, end="", sep="")
                   # Haploids
                    elif genotype.startswith("0"):
                        print("\t", 1, "\t", 0, end="", sep="")
                    elif genotype.startswith("1"):
                        print("\t", 0, "\t", 1, end="", sep="")
                    else:
                        print("\t", 0, "\t", 0, end="", sep="")

                # Use read count information
                else:
                    fields = line[idx].split(":")
                    # Assume sample format is GT:AD:...
                    if len(fields) > 1:
                        ad = fields[1].split(",")
                        print("\t", ad[0], "\t", ad[1], end="", sep="")
                    else:
                        print("\t", "0", "\t", "0", end="", sep="")

            # Newline
            print()

            # Update position information
            last_position = int(line[1])
            last_map_coord = map_coord

