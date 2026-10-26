"""
Take sample-specific Ancestry HMM output files and aggregate them together
into a table holding maximum a posteriori ancestry `calls` for each sample.
All input files must have the same length and the same chrom/position entry
at each site.

Input file format (one/sample)
------------------------------
chrom   position    2,0     1,1             0,2
2L      5845        1       1.51388e-07     2.25402e-14
...

Output file format
------------------
chrom   position    sample0     sample1    ...
2L      5845        0           -1         ...
...

Here, -1 denotes a site where `heterozygous` ancestry 1,1 has the highest
posterior probability. We ignore these calls and focus on heterozygous ones.
"""

import argparse
import csv
import numpy as np


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, nargs="*",
                        help="Single-sample Ancestry HMM output files")
    parser.add_argument("-s", "--samples", required=True,
                        help="Sample file used to generate inputs")
    parser.add_argument("-o", "--output", required=True,
                        help="Output filename")
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()

    # Read sample names
    samples = []
    with open(args.samples, "r") as fin:
        for line in csv.reader(fin, delimiter="\t"):
            if line[1] == "admixed":
                samples.append(line[0])
    samples = samples[:len(args.input)]

    # Codes for two ancestries. I use `codes[np.argmax(probs)]`, where `probs`
    # are site posterior probabilities, to retrieve ancestry codes.
    codes = [0, 0.5, 1]
    # Load input files
    calls = dict()
    for input_file in args.input:
        with open(input_file, "r") as fin:
            for line in csv.reader(fin, delimiter="\t"):
                # Skip the header line
                if line[0] == "chrom":
                    continue
                chrom, pos = line[:2]
                if (chrom, pos) not in calls:
                    calls[(chrom, pos)] = list()
                probs = [float(x) for x in line[2:5]]
                code = codes[np.argmax(probs)]
                calls[(chrom, pos)].append(code)

    # Check sample counts
    n_samples = len(samples)
    for chrom_pos in calls:
        assert n_samples == len(calls[chrom_pos])

    # Emit an output file
    with open(args.output, "w") as fout:
        header = "chrom\tposition\t" + "\t".join(samples) + "\n"
        fout.write(header)
        for (chrom, pos), _calls in calls.items():
            line = "\t".join([chrom, pos] + [str(x) for x in _calls]) + "\n"
            fout.write(line)

