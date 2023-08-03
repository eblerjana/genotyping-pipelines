# Create PanGenie-ready VCF from Minigraph-Cactus VCF

This pipeline was used in the HPRC paper (https://www.nature.com/articles/s41586-023-05896-x) in order to generate PanGenie-ready VCFs from the Minigraph-Cactus VCFs. It performs the following steps:

* filter out all non-top level bubbles (LV>1)
* remove bubbles for which more than 20% of the haplotypes carry a missing allele (".")
* decompose bubbles to determine all nested variant alleles and annotate VCF accordingly

## Required inputs

**vcf**: VCF file produced from the MC graph using `` vg decompose `` (as done in the HPRC paper)   
**gfa**: corresponding GFA from which the vcf was generated

## Outputs

* a **multi-allelic** graph VCF called `` results/vcf/{callsetname}/{callsetname}_filtered_ids.vcf`` representing bubbles in the graph

* a **bi-allelic** callset VCF called `` results/vcf/{callsetname}/{callsetname}_filtered_ids_biallelic.vcf`` representing variant alleles contained in the graph


### How to use the output VCFs with PanGenie

The **multi-allelic** VCF can be used as input to PanGenie (https://github.com/eblerjana/pangenie/tree/master), e.g.:

``` bat

PanGenie -i <input-reads> -v results/vcf/{callsetname}/{callsetname}_filtered_ids.vcf -r <reference-genome> -o pangenie -j 24 -t 24

```

VCFs produced by this pipeline contain special annotations in the INFO field ("ID"). In the multi-allelic VCF, each record defines a bubble in the pangenome graph. Each allele of a bubble is annotated by a sequence of IDs in the ID field, separated by a colon, which define variants nested inside of these bubbles. The bi-allelic VCF contains one record for each such individual ID. Both VCFs provide different representation of the same genetic variation present in the graph. After genotyping the bubbles with the command above, one can convert the bubble genotypes into genotypes for all nested variants using the script [convert-to-biallelic.py](https://bitbucket.org/jana_ebler/hprc-experiments/src/master/genotyping-experiments/workflow/scripts/convert-to-biallelic.py) and the following command:

``` bat

cat pangenie_genotyping.vcf | python3 convert-to-biallelic.py {callsetname}_filtered_ids_biallelic.vcf > pangenie_genotyping_biallelic.vcf

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