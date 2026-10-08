#!/bin/bash
set -e
zcat ../data/DGRP-3_final_biallelic-new.vcf.gz | head -n 20000 > short_DGRP.vcf
bgzip -f short_DGRP.vcf
bcftools index short_DGRP.vcf.gz
# extract sample names
bcftools query -l short_DGRP.vcf.gz > DGRP.samples

zcat ../data/ZI_ET_FR_combined_all_num_biallelic.vcf.gz | head -n 20000 > short_ZI_ET_FR_comb.vcf
bgzip -f short_ZI_ET_FR_comb.vcf
bcftools index short_ZI_ET_FR_comb.vcf.gz
# get sample names for FR and ZI
bcftools query -l short_ZI_ET_FR_comb.vcf.gz > ZI_ET_FR.samples
grep FR ZI_ET_FR.samples > FR.samples
grep ZI ZI_ET_FR.samples > ZI.samples
cat ZI.samples FR.samples > ZI_FR.samples

# Trim DGRP panel file
bcftools annotate -x INFO,^FORMAT/GT,FORMAT/AD,FORMAT/GQ short_DGRP.vcf.gz \
    -o trimmed_DGRP.vcf.gz

# subset and trim panel VCF file
bcftools view -S ZI_FR.samples short_ZI_ET_FR_comb.vcf.gz | \
    bcftools annotate -x INFO,^FORMAT/GT,FORMAT/AD,FORMAT/GQ \
    -o trimmed_ZI_FR_comb.vcf.gz

# Need to fix contig IDs for this to work!
# Merge trimmed VCF files
bcftools index trimmed_DGRP.vcf.gz
bcftools index trimmed_ZI_FR_comb.vcf.gz
bcftools merge trimmed_DGRP.vcf.gz trimmed_ZI_FR_comb.vcf.gz | \
    bcftools view -m2 -M2 -v snps -o merged_ZI_FR_DGRP.vcf.gz
