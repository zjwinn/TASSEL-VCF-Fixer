# TASSEL-VCF-Fixer
The Trait Analysis by aSSociation, Evolution and Linkage (TASSEL) 5.0 genotyping-by-sequencing (GBS) standalone pipeline version 2.0 produces a variant calling format (VCF) file which is always coded a major/minor allele. This program, written in Bash, is made to fix this and recode the alleles in reference/alternative. The program has three dependencies:

-vcftools
-bcftools
-bedtools

These combined tools are used to reference an indexed reference sequence file (.fa) and identify which single nucleotide polymorphisms (SNPs) in a VCF are coded properly in reference/alternate format and which must be changed to refelect a true reference/alternate issue. For help on running the program, refer to the usage function by running the shell script (.sh) with the help option (-h,-help). 
