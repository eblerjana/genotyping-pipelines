# PanGenie genotyping pipeline

This pipeline can be used to perform different genotyping experiments with PanGenie. 


## What the pipeline can do

**leave-one-out experiment**: PanGenie's genotyping performance is evaluated by removing one of the panel samples from the input VCF, and then genotyping it using a panel containing the remaining n-1 samples only. PanGenie predicts a genotype for each input variant. The resulting VCF is converted into the biallelic VCF representation, based on the annotations that it contains in the INFO field. This results in a genotyped version of the variants in the biallelic VCF. The genotypes are next compared to the ground truth genotypes of the left out sample using the weighted genotype concordance as a metric.

**population-typing**: PanGenie is run on a large set of samples (e.g. the 1000 Genomes samples). The multiallelic VCF is used as input for PanGenie, producing genotype predictions for all variants it contains. Genotyped VCFs are again converted to the biallelic VCF representation, based on the annotations in the multiallelic VCF, producing genotypes for all alleles contained in the biallelic VCF. Results are then evaluated based on trio information and different statistics, and a positive set (=strict) as well as a filtered set (=lenient) are computed containing subsets of variants genotyped with high quality.


## How to set up

In order to run the pipeline, paths to required input data must be provided in the config file as explained below. Required is a PanGenie-ready input VCF (annotated multi-allelic and bi-allelic versions), as provided below in section "Existing datasets", or as produced by this pipeline: https://bitbucket.org/jana_ebler/vcf-merging/src/master/. The following sections in the config file need to be filled:


```yaml

callsetname:
  # PanGenie-ready, multi-allelic VCF file
  multi: "/path/to/graph.vcf"
  # PanGenie-ready, bi-allelic VCF file
  bi: "/path/to/callset.vcf"
  # reference genome in FASTA format
  reference: "/path/to/reference.fa"
  # variants contained in the callset. Options are: snp|indels|large-deletion|large-insertion|large-complex
  variants:
    - snp 
    - indels
    - large-deletion
    - large-insertion
    - large-complex
  # repeat annotations in BED format (see resources/ folder for GRCh38 and CHM13-based annotations that can be used here)
  repeat_regions: "/path&to/annotations.bed"
  # if leave-one-out experiment shall be run, specify samples to run it on. Otherwise, leave empty.
  leave_one_out_samples: []
   
# file with information on sample and read data in PED format. Required columns (in this order):
# FamilyID	SampleID	FatherID	MotherID	Population	Superpopulation	Sample_Illumina
# FamilyID: specifies the name of a trio
# SampleID: name of the child sample
# FatherID: name of the father (0 if not available)
# MotherID: name of the mother (0 if not available)
# Population: which population the sample is from
# Superpopulation: which superpopulation (AFR, AMR, EAS, EUR, SAS)
# Sample_Illumina: path to a FASTA/FASTQ file with Illumina reads of the child sample
reads: "resources/genotyping-pilot-reads.tsv"


# PanGenie command. Different versions can be run by listing several commandlines.
pangenie:
 pangenie.v100.subsample14: "singularity exec --bind /:/hilbert container-main.sif PanGenie"
 pangenie.v201.subsample5: "singularity exec --bind /:/hilbert eblerjana_eblerjana_pangenie-v2.1.0.sif PanGenie -a 5"


# Downsampling coverages for leave-one-out experiment. If reads shall not be downsampled, leave empty.
downsampling: []

```

For example:

```yaml

callsets:
 cactus-hg38:
   multi: "/gpfs/project/projects/medbioinf/users/ebler/minigraph-cactus-paper/genotyping-experiments-hg38/results/data/vcf/cactus-100000/cactus_filtered_ids.vcf.gz"
   bi: "/gpfs/project/projects/medbioinf/users/ebler/minigraph-cactus-paper/genotyping-experiments-hg38/results/data/vcf/cactus-100000/cactus_filtered_ids_biallelic.vcf.gz"
   reference: "/gpfs/project/projects/medbioinf/users/ebler/minigraph-cactus-paper/genotyping-experiments-hg38/results/data/fasta/hg38.fa"
   variants:
     - snp 
     - indels
     - large-deletion
     - large-insertion
     - large-complex
   repeat_regions: "resources/ucsc-simple-repeats.bed"
   leave_one_out_samples:
     - HG00438
     - HG02717
     - HG00733
     - NA20129
     - HG03453
   
reads: "genotyping-reads.tsv"


# PanGenie command. Example on how to specify singularity containers to be used.
pangenie:
 pangenie.v100: "singularity exec --bind /:/hilbert pangenie-v100.sif PanGenie"
 pangenie.v201: "singularity exec --bind /:/hilbert pangenie-v210.sif PanGenie"


# Downsampling coverages for leave-one-out experiment. This example downsamples reads to 10x and 20x and runs leave-one-out experiments on these coverages in addition to the full coverage data.
downsampling:
 - 10
 - 20
```

## Existing datasets


We have already produced input reference panels for several datasets from high-quality, haplotype-resolved assemblies that can be used as input to PanGenie and this pipeline. They are available from here:

| Dataset | PanGenie input VCF        |  Callset VCF         | 
|-------------| :-------------: |:-------------:| 
| HGSVC-GRCh38 (freeze3, 64 haplotypes) | [graph-VCF](https://zenodo.org/record/7763717/files/pav-panel-freeze3.vcf.gz?download=1) | [callset-VCF](https://zenodo.org/record/7763717/files/pav-calls-freeze3.vcf.gz?download=1) | 
| HGSVC-GRCh38 (freeze4, 64 haplotypes) |  [graph-VCF](https://zenodo.org/record/7763717/files/pav-panel-freeze4.vcf.gz?download=1)     | [callset-VCF](https://zenodo.org/record/7763717/files/pav-calls-freeze4.vcf.gz?download=1) | 
| HPRC-GRCh38 (88 haplotypes) | [graph-VCF](https://zenodo.org/record/6797328/files/cactus_filtered_ids.vcf.gz?download=1)     |  [callset-VCF](https://zenodo.org/record/6797328/files/cactus_filtered_ids_biallelic.vcf.gz?download=1)    | 
| HPRC-CHM13 (88 haplotypes) | [graph-VCF](https://zenodo.org/record/7839719/files/chm13_cactus_filtered_ids.vcf.gz?download=1) | [callset-VCF](https://zenodo.org/record/7839719/files/chm13_cactus_filtered_ids_biallelic.vcf.gz?download=1)   | 



## How to run the pipeline

Paths to input files needed must be specified in the config file: `` config/config.yaml `` as explained in the previous section.
The whole pipeline can then be run using the following command:

``  snakemake --use-conda -j <number of cores>  `` 

Alternatively, different parts of the pipeline can be run individually:

* ``  snakemake leave_one_out --use-conda -j <number of cores>  ``  runs the PanGenie leave-one-out experiment.
* ``  snakemake population_typing --use-conda -j <number of cores>  `` runs PanGenie on all 3,202 genomes of the 1000 Genomes Project and analyzes the results. This also includes defining a positive (=strict) and filtered (=lenient) set of genotypes.
