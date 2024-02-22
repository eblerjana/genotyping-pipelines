# SV count comparisons

Given multi-sample VCFs with SV calls and genotypes, this pipeline compares the number of SV alleles carried by each sample. Callsets are compared based on the intersection of the samples they contain.

## What this pipeline can do

The pipeline compares different multi-sample SV callsets based on the number of SV alleles carried by each sample. Due to different ways of representing SV events, especially SVs in more complex or repetitive sequence contexts, comparing different SV callsets is not trivial. This pipeline therefore compares callsets based on the number of SV alleles carried by each sample, circumventing the need of directly comparing SV records across callsets.

## How to set up

Provide paths to the input data in the ``config/config.yaml``:

``` bat

# name of the output folder to write all results to
results: "result"

# callsets to consider. vcf specifies the path to the VCF file, collapse specifies whether or not to collapse overlapping
# variant alleles prior to evaluation

callsets:
 <callsetname1>:
   vcf: "path/to/callset1.vcf.gz"
   collapse: True
   reference: "path/to/reference.fa"
 <callsetname2>:
   vcf: "path/to/callset2.vcf.gz"
   collapse: False
   reference: "path/to/reference.fa"


# file with information on sample and read data in PED format. Required columns (in this order):
# SampleID	Superpopulation
# SampleID: name of the child sample
# Superpopulation: which superpopulation (AFR, AMR, EAS, EUR, SAS)
populations: "path/to/populations.tsv"

```
For each callset, provide a vcf file containing SV calls. When ``collapse`` is set to `` True ``, truvari collapse is run on the callset prior to analysis in order to merge similar alleles. If set to `` False ``, the calls are used as they are. In addition to the callsets, a tsv file needs to be provided (``populations``) that specifies the population each sample belongs to (format described above). Make sure to list all samples here that occur in the vcf files.

Note that the comparison is done based on the intersection of samples present in all vcf files to make sure results are comparable. This pipeline is not intended for cases with disjoint sample sets.

## How to run the pipeline

``` bat

snakemake -j <nr_cores>

```
