configfile: "config/config.yaml"


sample_index_unrelated = config['sample_index_unrelated']
sample_index_related = config['sample_index_related']
read_dir_unrelated = config['read_dir_unrelated']
read_dir_related = config['read_dir_related']
pangenie = config['pangenie']

# read paths to sequencing data from files
lsamples = {}

# read 2504 unrelated samples
for line in open(sample_index_unrelated, 'r'):
	if line.startswith('study'):
		continue
	fields = line.split('\t')
	run_id = fields[3]
	sample_name = fields[2]
	lsamples[sample_name] = read_dir_unrelated + sample_name + '_' + run_id + '.fasta.gz'

# read 689 unrelated samples
for line in open(sample_index_related, 'r'):
	if line.startswith('study'):
		continue
	fields = line.split('\t')
	run_id = fields[3]
	sample_name = fields[2]
	lsamples[sample_name] = read_dir_related + sample_name + '_' + run_id + '.fasta.gz'


variants = ['snp', 'indels', 'large-deletion', 'large-insertion', 'large-complex']
lsources = [s for s in config['graph_vcf'].keys()]


################################################################
######   prepare input panel and ground truth genotypes  #######
################################################################

# remove positions that are ".|." in the left out sample. These cannot be used for evaluation, as the true genotype is unknown
# for now, also remove CHM13, because PanGenie cannot handle haploid samples
# remove sample from the panel
rule remove_missing:
	input:
		lambda wildcards: config['graph_vcf'][wildcards.source] if wildcards.callset == 'multi' else config['biallelic_vcf'][wildcards.source] 
	output:
		"results/leave-one-out/{source}/{sample}/preprocessed-vcfs/{source}_{callset}_no-missing.vcf"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=20000
	priority: 1
	wildcard_constraints:
		callset = "bi|multi"
	shell:
		"zcat {input} | python3 workflow/scripts/remove-missing.py {wildcards.sample} > {output}"


rule prepare_panel:
	input:
		vcf="results/leave-one-out/{source}/{sample}/preprocessed-vcfs/{source}_multi_no-missing.vcf.gz"
	output:
		"results/leave-one-out/{source}/{sample}/input-panel/panel-{sample}.vcf"
	conda:
		"../envs/genotyping.yml"
	priority: 1
	log:
		"results/leave-one-out/{source}/{sample}/input-panel/panel-{sample}.log"
	resources:
		mem_total_mb=20000
	shell:
		"bcftools view --samples ^{wildcards.sample} {input.vcf} | bcftools view --min-ac 1 2> {log} 1> {output}"



# extract ground truth genotypes for sample
rule prepare_truth:
	input:
		vcf= "results/leave-one-out/{source}/{sample}/preprocessed-vcfs/{source}_bi_no-missing.vcf.gz"
	output:
		"results/leave-one-out/{source}/{sample}/truth/truth-{sample}.vcf"
	conda:
		"../envs/genotyping.yml"
	priority: 1
	resources:
		mem_total_mb=20000
	log:
		"results/leave-one-out/{source}/{sample}/truth/truth-{sample}.log"
	shell:
		"bcftools view --samples {wildcards.sample} {input.vcf} 2> {log} 1> {output}"



########################################################
##################    run PanGenie    ##################
########################################################


# run pangenie
rule pangenie:
	input:
		reads= lambda wildcards: lsamples[wildcards.sample],
		fasta= lambda wildcards: config['reference'][wildcards.source],
		vcf="results/leave-one-out/{source}/{sample}/input-panel/panel-{sample}.vcf",
	output:
		reads = temp("results/leave-one-out/{source}/{sample}/pangenie/reads.fa"),
		path_segments = temp("results/leave-one-out/{source}/{sample}/pangenie/pangenie-{sample}_path_segments.fasta"),
		genotyping = "results/leave-one-out/{source}/{sample}/pangenie/pangenie-{sample}_genotyping.vcf"
	log:
		"results/leave-one-out/{source}/{sample}/pangenie/pangenie-{sample}.log"
	threads: 24
	resources:
		mem_total_mb=100000,
		runtime_hrs=3,
		runtime_min=1
	priority: 1
	params:
		out_prefix="results/leave-one-out/{source}/{sample}/pangenie/pangenie-{sample}"
	shell:
		"""
		gunzip -c {input.reads} > {output.reads}
		module load Singularity
		(/usr/bin/time -v {pangenie} -i {output.reads} -v {input.vcf} -r {input.fasta} -o {params.out_prefix} -s {wildcards.sample} -j {threads} -t {threads} -g ) &> {log}
		"""


########################################################
##################    Evaluation      ##################
########################################################


rule alleles_per_bubble:
	input:
		lambda wildcards: config['graph_vcf'][wildcards.source]
	output:
		plot="results/leave-one-out/{source}/alleles-per-bubble.pdf",
		bed="results/leave-one-out/{source}/complex-bubbles.bed"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=20000,
		runtime_hrs=1
	shell:
		"zcat {input} | python3 workflow/scripts/variant-statistics.py {output.plot} 1 > {output.bed}"


# prepare beds for biallelic and complex graph regions
rule prepare_beds:
	input:
		bed="results/leave-one-out/{source}/complex-bubbles.bed",
		fai= lambda wildcards: config['reference'][wildcards.source] + ".fai"
	output:
		bed="results/leave-one-out/{source}/biallelic-bubbles.bed",
		tmp=temp("results/leave-one-out/{source}/biallelic-bubbles.fai")
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		sort -k1,1d -k 2,2n -k 3,3n {input.fai} > {output.tmp}
		bedtools complement -i {input.bed} -g {output.tmp} > {output.bed}
		"""


# convert genotyped VCF to biallelic representation
rule convert_genotypes_to_biallelic:
	input:
		vcf="results/leave-one-out/{source}/{sample}/pangenie/pangenie-{sample}_genotyping.vcf.gz",
		biallelic= lambda wildcards: config['biallelic_vcf'][wildcards.source]
	output:
		"results/leave-one-out/{source}/{sample}/pangenie/pangenie-{sample}_genotyping-biallelic.vcf"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=30000
	priority: 1
	shell:
		"zcat {input.vcf} | python3 workflow/scripts/convert-to-biallelic.py {input.biallelic} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n \"}}' > {output}"


# determine untypable ids
rule untypable_ids:
	input:
		biallelic= lambda wildcards: config['biallelic_vcf'][wildcards.source]
	output:
		lists="results/leave-one-out/{source}/untypable-ids/{sample}-untypable.tsv",
		summary="results/leave-one-out/{source}/untypable-ids-{sample}.tsv"
	params:
		out="results/leave-one-out/{source}/untypable-ids"
	conda:
		"../envs/genotyping.yml"
	shell:
		"zcat {input.biallelic} | python3 workflow/scripts/untypable-ids-single.py {params.out} {wildcards.sample} > {output.summary}"


# determine untypable IDs
rule remove_untypable:
	input:
		vcf="results/leave-one-out/{source}/{sample}{other}.vcf",
		ids="results/leave-one-out/{source}/untypable-ids/{sample}-untypable.tsv"
	output:
		vcf="results/leave-one-out/{source}/{sample}{other}-typable-{vartype}.vcf.gz",
		tbi="results/leave-one-out/{source}/{sample}{other}-typable-{vartype}.vcf.gz.tbi"
	wildcard_constraints:
		sample = "|".join([s for s in lsamples.keys()]),
		vartype="|".join(variants)
	resources:
		mem_total_mb=20000,
		runtime_hrs=1
	priority: 1
	shell:
		"""
		cat {input.vcf} | python3 workflow/scripts/skip-untypable.py {input.ids} | python3 workflow/scripts/extract-varianttype.py {wildcards.vartype} | bgzip -c > {output.vcf}
		tabix -p vcf {output.vcf}
		"""

rule rtg_format:
	input:
		lambda wildcards: config['reference'][wildcards.source]
	output:
		directory("results/leave-one-out/{source}/SDF")
	resources:
		mem_total_mb=20000
	priority: 1
	shell:
		'rtg format -o {output} {input}'


def region_to_bed(wildcards):
	if wildcards.regions == "biallelic":
		return "results/leave-one-out/{source}/biallelic-bubbles.bed".format(source=wildcards.source)
	if wildcards.regions == "multiallelic":
		return "results/leave-one-out/{source}/complex-bubbles.bed".format(source=wildcards.source)
	assert(False)


# precision-recall
rule vcfeval:
	input:
		callset="results/leave-one-out/{source}/{sample}/pangenie/pangenie-{sample}_genotyping-biallelic-typable-{vartype}.vcf.gz",
		callset_tbi="results/leave-one-out/{source}/{sample}/pangenie/pangenie-{sample}_genotyping-biallelic-typable-{vartype}.vcf.gz.tbi",
		baseline="results/leave-one-out/{source}/{sample}/truth/truth-{sample}-typable-{vartype}.vcf.gz",
		baseline_tbi="results/leave-one-out/{source}/{sample}/truth/truth-{sample}-typable-{vartype}.vcf.gz.tbi",
		regions= region_to_bed,
		sdf="results/leave-one-out/{source}/SDF"
	output:
		summary="results/leave-one-out/{source}/{sample}/evaluation/precision-recall-typable/pangenie/{regions}_{vartype}/summary.txt"
	conda:
		"../envs/genotyping.yml"
	priority: 1
	wildcard_constraints:
		sample = "|".join(lsamples),
		regions = "biallelic|multiallelic",
		vartype = "|".join(variants)
	params:
		tmp = "results/leave-one-out/{source}/{sample}/evaluation/precision-recall-typable/pangenie/{regions}_{vartype}_tmp",
		outname = "results/leave-one-out/{source}/{sample}/evaluation/precision-recall-typable/pangenie/{regions}_{vartype}",
		which = "--all-records"
	resources:
		mem_total_mb=20000,
		runtime_hrs=1,
		runtime_min=40
	shell:
		"""
		rtg vcfeval -b {input.baseline} -c {input.callset} -t {input.sdf} -o {params.tmp} --ref-overlap --evaluation-regions {input.regions} {params.which} --Xmax-length 30000 > {output.summary}.tmp
		mv {params.tmp}/* {params.outname}/
		mv {output.summary}.tmp {output.summary}
		rm -r {params.tmp}
		"""


# determine the variants that went into re-typing per category
rule collected_typed_variants:
	input:
		callset="results/leave-one-out/{source}/{sample}/preprocessed-vcfs/{source}_bi_no-missing.vcf.gz",
		regions= region_to_bed,
		ids="results/leave-one-out/{source}/untypable-ids/{sample}-untypable.tsv"
	output:
		"results/leave-one-out/{source}/{sample}/genotyped-ids/{regions}_{vartype}.tsv"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		sample = "|".join(lsamples),
		regions = "biallelic|multiallelic",
		vartype = "|".join(variants)
	resources:
		mem_total_mb=50000
	priority: 1
	shell:
		"zcat {input.callset} | python3 workflow/scripts/skip-untypable.py {input.ids} | python3 workflow/scripts/extract-varianttype.py {wildcards.vartype} | bedtools intersect -header -a - -b {input.regions} -u -f 0.5 | python3 workflow/scripts/get_ids.py > {output}"


# compute concordances
rule genotype_concordances:
	input:
		callset="results/leave-one-out/{source}/{sample}/pangenie/pangenie-{sample}_genotyping-biallelic-typable-{vartype}.vcf.gz",
		callset_tbi="results/leave-one-out/{source}/{sample}/pangenie/pangenie-{sample}_genotyping-biallelic-typable-{vartype}.vcf.gz.tbi",
		baseline="results/leave-one-out/{source}/{sample}/truth/truth-{sample}-typable-{vartype}.vcf.gz",
		baseline_tbi="results/leave-one-out/{source}/{sample}/truth/truth-{sample}-typable-{vartype}.vcf.gz.tbi",
		regions= region_to_bed,
		typed_ids = "results/leave-one-out/{source}/{sample}/genotyped-ids/{regions}_{vartype}.tsv"
	output:
		tmp_vcf1=temp("results/leave-one-out/{source}/{sample}/evaluation/concordance/pangenie_{regions}_{vartype}_base.vcf"),
		tmp_vcf2=temp("results/leave-one-out/{source}/{sample}/evaluation/concordance/pangenie_{regions}_{vartype}_call.vcf"),
		summary="results/leave-one-out/{source}/{sample}/evaluation/concordance/pangenie/{regions}_{vartype}/summary.txt"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		sample = "|".join(lsamples),
		regions = "biallelic|multiallelic",
		vartype = "|".join(variants)
	log:
		"results/leave-one-out/{source}/{sample}/evaluation/concordance/pangenie/{regions}_{vartype}/summary.log"
	resources:
		mem_total_mb=40000,
		runtime_hrs=0,
		runtime_min=40
	priority: 1
	shell:
		"""
		bedtools intersect -header -a {input.baseline} -b {input.regions} -u -f 0.5 | bgzip > {output.tmp_vcf1}
		bedtools intersect -header -a {input.callset} -b {input.regions} -u -f 0.5 | bgzip > {output.tmp_vcf2}
		python3 workflow/scripts/genotype-evaluation.py {output.tmp_vcf1} {output.tmp_vcf2} {input.typed_ids} --qual 0 2> {log} 1> {output.summary}
		"""


########################################################
##################     Plotting      ###################
########################################################

# collect results across all samples
rule collect_results:
	input:
		lambda wildcards: expand("results/leave-one-out/{{source}}/{sample}/evaluation/{{metric}}/pangenie/{{regions}}_{{vartype}}/summary.txt", sample = config['leave_one_out_samples'][wildcards.source])
	output:
		"results/leave-one-out/{source}/plots/{metric}/pangenie/{metric}_{regions}_{vartype}.tsv"
	params:
		samples = lambda wildcards: ','.join(config['leave_one_out_samples'][wildcards.source]),
		outfile = 'results/leave-one-out/{source}/plots/{metric}/pangenie/{metric}'
	priority: 1
	shell:
		"python3 workflow/scripts/collect-results.py {wildcards.metric} {wildcards.source} {params.samples} {wildcards.regions} -variants {wildcards.vartype} -folder results/leave-one-out/ -outfile {params.outfile}"



# plot results using graph regions (biallelic/multiallelic) as stratifications
rule plotting_graph:
	input:
		lambda wildcards: expand("results/leave-one-out/{{source}}/plots/{{metric}}/pangenie/{{metric}}_{regions}_{vartype}.tsv", regions=["biallelic", "multiallelic"], vartype=config['variants'][wildcards.source])
	output:
		"results/leave-one-out/{source}/plots/{metric}/pangenie/{metric}_graph.pdf"
	wildcard_constraints:
		metric="concordance|precision-recall-typable"
	priority: 1
	conda:
		"../envs/genotyping.yml"
	params:
		variants= lambda wildcards: ' '.join([v for v in config['variants'][wildcards.source]])
	shell:
		"python3 workflow/scripts/plot-results.py -files {input} -outname {output} -sources biallelic multiallelic -metric {wildcards.metric} -variants {params.variants}"


rule compress_vcf:
	input:
		"results/leave-one-out/{filename}.vcf"
	output:
		vcf="results/leave-one-out/{filename}.vcf.gz",
		tbi="results/leave-one-out/{filename}.vcf.gz.tbi"
	priority: 1
	shell:
		"""
		bgzip -c {input} > {output.vcf}
		tabix -p vcf {output.vcf}
		"""
