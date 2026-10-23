"""
Adapted from vcf2ahmm.py on the AncestryHMM github:
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
import csv
import gzip
import numpy as np
import pandas


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--vcf", type=str,
                        help="input vcf file", required=True)
    parser.add_argument("-s", "--samples", type=str,
                        help="sample to population file", required=True)
    parser.add_argument("-r", "--rec", type=str,
                        help="uniform recombination rate OR path to tab-separated recombination map file", required=True)
    parser.add_argument("-g", type=int,
                        help="boolean (0 = use sample reads| 1 = use genotypes)", default=0)
    parser.add_argument("-m", type=int,
                        help="minimum distance between successive snps (int,bp)", default=1000)
    parser.add_argument("-o",  "--out", type=str,
                        help="ploidy file for ahmm input", default="ahmm.ploidy")
    parser.add_argument("--min_total", type=int,
                        help="minimum number of reference population samples (int)", default=10)
    parser.add_argument("--min_diff",type=float,
                        help="minimum allele frequency difference between a pair of populations (float)", default=0.1)
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    # Get sample-to-population mappings
    sample_map = {}
    panel_pops = set()
    ordered_samples = []
    with open(args.samples) as file:
        for line in file:
            (key, val) = line.split()
            sample_map[key] = val
            panel_pops.add(val)
            if val == "admixed":
                ordered_samples.append(key)
    # Sort panel population indices
    panel_pops = sorted([p for p in panel_pops if p != "admixed"])

    ### IDs
    sample_fields = []

    ### keep track of where we are on the chromosome
    last_position = -1000000
    last_map_coord = 0
    current_chrom = "NA"

    ### have we printed the ploidy file yet
    ploidy_print = 0

    # Load recombination map file
    if args.rec.isnumeric():
        r = float(args.rec)
        rec_maps = None
    else:
        # If a recombination map was provided; make a dictionary of interpolation
        # functions (one for each chromosome)
        r = None
        rec_df = pandas.read_csv(args.rec, sep="\\s+")
        rec_maps = {}
        for chrom in set(rec_df["Chrom"]):
            sub_df = rec_df[rec_df["Chrom"] == chrom]
            starts = np.array(sub_df["Start"])
            stops = np.array(sub_df["Stop"])
            # Convert rates from cM/Mb to M/bp
            rates = np.array(sub_df["cM/Mb"]) * 1e-8
            # Calculate sites spanned by each interval
            spans_bp = stops - starts + 1
            spans_M = spans_bp * rates
            map_M = np.cumsum(spans_M)
            map_M = np.concatenate([[0], map_M])
            map_bp = np.concatenate([[starts[0]], stops])
            # Interpolate site map coordinates
            rec_maps[chrom] = lambda x: np.interp(x, map_bp, map_M)

    ## read the vcf file
    if args.vcf.endswith(".gz"):
        open_func = lambda f: gzip.open(f, "rt")
    else:
        open_func = lambda f: open(f, "r")

    with open_func(args.vcf) as fin:
        ### split on tab
        for line in csv.reader(fin, delimiter="\t"):
            ### parse the header line
            if (line[0].startswith('#CHROM')):

                ### get our info
                header = line
                admixed_idx = [header.index(sample) for sample in
                               ordered_samples if sample_map[sample] == "admixed"]
                continue

            ### skip lines that start with ##
            if (line[0].startswith('##')):
                continue

            ### remove sites that are not biallelic snps
            if (len(line[3]) != 1 or len(line[4]) != 1):
                continue

            ### now read the line and decide what to do with it
            if (current_chrom != line[0]):
                current_chrom = line[0]
                last_position = -1000000
                last_map_coord = 0

            ### if position too close skip
            elif (int(line[1]) - last_position < args.m):
                continue

            ### otherwise get population allele frequencies
            alts = {}
            totals = {}
            for i in range(9, len(line)):
                if (sample_idx[i] in sample_map):
                    if sample_map[sample_idx[i]] not in totals:
                        totals[sample_map[sample_idx[i]]] = 0
                        alts[sample_map[sample_idx[i]]] = 0
                    # Assume the genotype is the first three sample chars
                    genotype = line[i][:3]
                    if "|" in genotype:
                        genotype.replace("|", "/")
                    if genotype == "0/0":
                        totals[sample_map[sample_idx[i]]] += 2
                    elif genotype == "0/1":
                        alts[sample_map[sample_idx[i]]] += 1
                        totals[sample_map[sample_idx[i]]] += 2
                    elif genotype == "1/1":
                        alts[sample_map[sample_idx[i]]] += 2
                        totals[sample_map[sample_idx[i]]] += 2

                    ### read haploid genotypes as needed
                    elif genotype.startswith('1'):
                        alts[sample_map[sample_idx[s]]] += 1
                        totals[sample_map[sample_idx[s]]] += 1
                    elif genotype.startswith('0'):
                        totals[sample_map[sample_idx[s]]] += 1

            ### make sure we have enough samples at this position
            too_few = 0
            freqs = []
            for population, total in totals.items():
                if (population != "admixed"):
                    freqs.append(float(alts[population])/float(total))
                    if (total < args.min_total):
                        too_few += 1
            if (too_few > 0):
                continue

            ### check allele frequencies
            diff = 0
            for p1 in range(0, len(freqs)-1):
                for p2 in range(p1+1,len(freqs)):
                    if (abs(freqs[p1] - freqs[p2]) > args.min_diff):
                        diff += 1
            if (diff == 0):
                continue

            ### if we got this far, we're printing the site
            print(line[0], line[1], end="\t", sep="\t")
            for population in panel_pops:
                total = totals[population]
                if (population != "admixed"):
                    ref = total - alts[population]
                    print (ref, alts[population], end="\t", sep="\t")

            # Calculate map distance to the previous site
            if rec_maps is None:
                rec = float(float(line[1]) - last_position) * args.rec
            else:
                position = float(line[1])
                map_coord = rec_maps[current_chrom](position)
                rec = map_coord - last_map_coord
            print (rec, end="")

            # Now print the read counts for each admixed sample ...
            for idx in admixed_idx:
                # From genotypes
                if (args.g == 1):
                    genotype = line[idx][:3]
                    if "|" in genotype:
                        genotype = genotype.replace("|", "/")
                    if genotype == "0/0":
                        print("\t", 2, "\t", 0, end="", sep="")
                    elif genotype == "0/1":
                        print("\t", 1, "\t", 1, end="", sep="")
                    elif genotype == "1/1":
                        print("\t", 0, "\t", 2, end="", sep="")
                    elif genotype.startswith('0'):
                        print("\t", 1, "\t", 0, end="", sep="")
                    elif genotype.startswith('1'):
                        print("\t", 0, "\t", 1, end="", sep="")
                    else:
                        print("\t", 0, "\t", 0, end="", sep="")
                    if (ploidy_print == 0):
                        print(sample_idx[s], 2, sep="\t", file=open(args.out, "a"))

                # From read counts
                else:
                    fields = line[idx].split(sep=":")
                    if (len(fields) > 2):
                        ad = fields[1].split(sep=",")
                        print ("\t", ad[0], "\t", ad[1], end="", sep="")
                    else:
                        print ("\t", 0, "\t", 0,  end="", sep="")
                    if (ploidy_print == 0):
                        print(sample_idx[s], 2, sep="\t", file=open(args.out, "a"))

            ## newline
            print()

            ### update last position after printing
            last_position = int(line[1])
            last_map_coord = map_coord
            ploidy_print = 1
