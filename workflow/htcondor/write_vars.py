# chromosome arg segments
Xs = [f"X_{i}" for i in range(11)]
L2s = [f"2L_{i}" for i in range(11)]
R2s = [f"2R_{i}" for i in range(12)]
L3s = [f"3L_{i}" for i in range(13)]
R3s = [f"3R_{i}" for i in range(15)]

with open("variables.txt", "w") as fout:
    for i, arms1 in enumerate([Xs, L2s, R2s, L3s, R3s]):
        for arms2 in [Xs, L2s, R2s, L3s, R3s][i + 1:]:
            for arm1 in arms1:
                for arm2 in arms2:
                    varz = [
                        f"DGRP_{arm1}",
                        f"DGRP_{arm2}",
                        f"p_bins_{arm1}_{arm2}.tsv",
                        f"sites_{arm1}_{arm2}.tsv"
                    ]
                    fout.write(",".join([str(x) for x in varz]) + "\n")
