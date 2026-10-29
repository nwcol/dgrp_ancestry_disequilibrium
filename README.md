### Mosaic ancestry in the Drosophila Genetic Reference Panel (DGRP-3)

This repository contains a Snakemake pipeline and Python scripts to infer
ancestry tracts in the admixed *Drosophila melanogaster* Genetic Reference
Panel version 3 (DGRP-3) (Lyman et al, 
[2026 https://doi.org/10.64898/2026.05.29.728783](https://doi.org/10.64898/2026.05.29.728783)).
By inferring ancestry across the genome and computing metrics of association
between ancestry at pairs of loci, we can learn about the extent of `ancestry
epistasis' and putative incompatabilities between the parent populations of
DGRP-3.


## Requirements

To obtain the required Python packages, create a virtual environment and
install them using pip:
```
python -m venv my/venv
pip install -r requirements.txt
```

You will also need Ancestry HMM, a command-line tool for inferring genomic
ancestry tracts.
Install Ancestry HMM by following the directions found here:
[https://github.com/russcd/Ancestry\_HMM/tree/master](https://github.com/russcd/Ancestry_HMM/tree/master).


## Running the pipeline

To run the pipeline, first download the DGRP-3 VCF file from 
[https://quantgenet.msu.edu/dgrp3/downloads.html](https://quantgenet.msu.edu/dgrp3/downloads.html)
and place it in the [data/vcfs](data/vcfs) directory.
Also download the Zambia (ZI), Ethiopia (ET), and France (FR) ancestral
reference panel from XXX [not currently available online] and place it in the
same directory.
With these files in hand, you can run the [Snakefile](workflow/Snakefile).
This will simplify/merge the DRGP-3 and panel VCF files, process genotypes into
the input format required by Ancestry HMM, run ancestry HMM, then analyze
inferred ancestry associations.


