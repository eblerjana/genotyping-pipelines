# Create PanGenie-ready VCF from Minigraph-Cactus VCF

This pipeline was used in the HPRC paper (https://www.nature.com/articles/s41586-023-05896-x) in order to generate PanGenie-ready VCFs from the Minigraph-Cactus VCFs for human (produced based on ``vg decompose`` and were **not** decomposed with vcfwave).

## What the pipeline can do

It performs the following steps:

* remove bubbles for which more than 20% of the haplotypes carry a missing allele (".")
* decompose bubbles to determine all nested variant alleles and annotate VCF accordingly

The pipeline is designed for human data.

## How to set up

Specify paths to the required input files in the `` config.yaml `` file:

``` bat

# name of the folder to write the results to
results: "results"

# callsets generated by Minigraph-Cactus pipeline (already filtered with vcfbub)
callsets:
 <callsetname>:
  # Minigraph-Cactus VCF
  vcf: "/path/to/vcf.vcf.gz"
  # Minigraph-Cactus graph GFA
  gfa: "/path/to/gfa.gfa"
  # TSV file specifying sex (1=male, 2=female) of each sample in input VCF. Format: <sample-name> <sex>. 
  sample_info: "/path/to/sample-info.tsv"
  # prefix used for chromosome names in FASTA file. 
  reference_prefix: "chr"
  # Choose between: CHM13 | GRCh38
  reference_version: "CHM13"    

```


## Required input data

### VCF
VCF file produced by the [Minigraph-Cactus pipeline](https://github.com/ComparativeGenomicsToolkit/cactus). The VCF must be preprocessed by `` vcfbub `` but **must not** be decomposed with vcfwave. Such VCFs can be produced with Minigraph-Cactus using the ``--vcf`` option. **Note:** VCFs decomposed with vcfwave (option ``--vcfwave`` in Minigraph-Cactus) are **not** compatible with this pipeline!


### GFA
corresponding GFA (gzip compressed) from which the vcf was generated


## Outputs

* a **multi-allelic** graph VCF called `` {results}/vcf/{callsetname}/{callsetname}_filtered_ids.vcf`` representing bubbles in the graph

* a **bi-allelic** callset VCF called `` {results}/vcf/{callsetname}/{callsetname}_filtered_ids_biallelic.vcf`` representing variant alleles contained in the graph


### How to use the output VCFs with PanGenie

The **multi-allelic** VCF can be used as input to PanGenie (https://github.com/eblerjana/pangenie/tree/master), e.g.:

``` bat
PanGenie-index -v {results}/vcf/{callsetname}/{callsetname}_filtered_ids.vcf -r <reference-genome> -t 24 -o index
PanGenie -f index -i <input-reads> -o pangenie -j 24 -t 24
```

VCFs produced by this pipeline contain special annotations in the INFO field ("ID"). In the multi-allelic VCF, each record defines a bubble in the pangenome graph. Each allele of a bubble is annotated by a sequence of IDs in the ID field, separated by a colon, which define variants nested inside of these bubbles. The bi-allelic VCF contains one record for each such individual ID. Both VCFs provide different representation of the same genetic variation present in the graph. After genotyping the bubbles with the command above, one can convert the bubble genotypes into genotypes for all nested variants using the script [convert-to-biallelic.py](https://bitbucket.org/jana_ebler/hprc-experiments/src/master/genotyping-experiments/workflow/scripts/convert-to-biallelic.py) and the following command:

``` bat

cat pangenie_genotyping.vcf | python3 convert-to-biallelic.py {results}/vcf/{callsetname}/{callsetname}_filtered_ids_biallelic.vcf > pangenie_genotyping_biallelic.vcf

```

This is often useful for evaluation of the genotypes, since bubbles are often really large and might contain hundreds of variants nested inside. In the bi-allelic VCF, all these nested variants will be listed, which makes comparison to existing callsets possible.

**Note:** do not use the bi-allelic VCF as input to PanGenie since PanGenie needs the bubble information encoded in the multi-allelic VCF. Instead, use the stradegy above to derive genotypes for all variants in the bi-allelic VCF.



## How to run the pipeline

Prepare the `` config.yaml `` file by adding paths to the input VCF/GFA, as well as information on these datasets (like the number of samples, reference genome, etc.). Then run the pipeline:

``` bat
snakemake -j <nr_cores>

```

## NOTE

This pipeline was specifically designed for the MC-graphs computed for HPRC. It requires the reference genome to be represented as a linear path through the graph, i.e. the reference genome is not allowed to loop back to itself. **Thus, this pipeline in its present form does not work for graphs produced by PGGB**.
