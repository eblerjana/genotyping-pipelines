# Create PanGenie-ready VCF from Minigraph-Cactus VCF

This pipeline was used in the HPRC paper (https://www.nature.com/articles/s41586-023-05896-x) in order to generate PanGenie-ready VCFs from the Minigraph-Cactus VCFs. It performs the following steps:

* filter out all non-top level bubbles (LV>1)
* remove bubbles for which more than 20% of the haplotypes carry a missing allele (".")
* decompose bubbles to determine all nested variant alleles and annotate VCF accordingly

## Required inputs

**vcf**: VCF file produced from the MC graph using `` vg decompose `` (as done in the HPRC paper)   
**gfa**: corresponding GFA from which the vcf was generated

## Outputs

* a **multi-allelic** graph VCF called `` results/vcf/{callsetname}/{callsetname}_filtered_ids.vcf`` that can be used as input to PanGenie

``` bat

PanGenie -i <input-reads> -v results/vcf/{callsetname}/{callsetname}_filtered_ids.vcf -r <reference-genome> -o pangenie -j 24 -t 24

```

* a **bi-allelic** callset VCF called `` results/vcf/{callsetname}/{callsetname}_filtered_ids_biallelic.vcf`` that can be used to transform the genotypes computed by PanGenie for all bubbles to genotypes for all nested variant alleles using the command below and this script [convert-to-biallelic.py](https://bitbucket.org/jana_ebler/hprc-experiments/src/master/genotyping-experiments/workflow/scripts/convert-to-biallelic.py):

``` bat

cat pangenie_genotyping.vcf | python3 convert-to-biallelic.py {callsetname}_filtered_ids_biallelic.vcf > pangenie_genotyping_biallelic.vcf

```

## How to run

Prepare the `` config.yaml `` file by adding paths to the input VCF/GFA, as well as information on these datasets (like the number of samples, reference genome, etc.). Then run the pipeline:

``` bat
snakemake -j <nr_cores>

```

## NOTE

This pipeline was specifically designed for the MC-graphs computed for HPRC. It requires the reference genome to be represented as a linear path through the graph, i.e. the reference genome is not allowed to loop back to itself. **Thus, this pipeline in its present form does not work for graphs produced by PGGB**.
