configfile: "config/config.yaml"

run_pilot = config['run_pilot']
trios = config['pilot_trios'] if run_pilot else config['trios'] 
sample_index_unrelated = config['sample_index_unrelated']
sample_index_related = config['sample_index_related']
repeats_bed = "results/data/bed/ucsc-simple-repeats.merged.bed"

samples_1000G = []

for s in open(sample_index_unrelated, 'r'):
	sample_name = s.split()[2]
	samples_1000G.append(sample_name)

for s in open(sample_index_related, 'r'):
	sample_name = s.split()[2]
	samples_1000G.append(sample_name)

# put the samples with additional data in list as well
for sample in config['additional_read_data'].keys():
	samples_1000G.append(sample)


rule collect_panel_samples:
	output:
		"results/population-typing/{source}/sample-index-panel.tsv"
	run:
		with open(output[0], "w") as outfile:
			for s in config['panel_samples'][wildcards.source]:
				outfile.write(s + '\n')
			

rule collect_samples:
	input:
		unrelated=sample_index_unrelated,
		related=sample_index_related,
		panel="results/population-typing/{source}/sample-index-panel.tsv"
	output:
		all="results/population-typing/{source}/sample-index.tsv",
		related="results/population-typing/{source}/sample-index-related.tsv",
		unrelated="results/population-typing/{source}/sample-index-unrelated.tsv"
	shell:
		"""
		grep -v "study_" {input.unrelated} | cut -f3 | cat - {input.panel} | sort | uniq > {output.unrelated}
		grep -v "study_" {input.related} | cut -f3 > {output.related}
		cat {output.unrelated} {output.related} > {output.all}
		"""


####################################################################
# extract necessary subsets of samples (e.g. all unrelated samples
####################################################################

rule extract_samples:
	input:
		vcf="results/population-typing/{source}/merged-vcfs/whole-genome/all-samples_bi_all.vcf.gz",
		samples="results/population-typing/{source}/sample-index-unrelated.tsv"
	output:
		"results/population-typing/{source}/merged-vcfs/whole-genome/unrelated-samples_bi_all.vcf.gz"
	resources:
		mem_total_mb = 50000,
		runtime_hrs = 10
	conda:
		"../envs/genotyping.yml"
	shell:
		"bcftools view --samples-file {input.samples} --force-samples {input.vcf} | bgzip -c > {output}"



################################################################
# HWE testing and analysis
################################################################

# compute filtered callsets
rule filtered_callsets:
	input:
		vcf="results/population-typing/{source}/merged-vcfs/whole-genome/{population}_bi_all.vcf.gz",
		filters="results/population-typing/{source}/evaluation/statistics/plot_bi_all_filters.tsv"
	output:
		"results/population-typing/{source}/merged-vcfs/filtered/{population}_bi_all_{filter}.vcf.gz"
	resources:
		mem_total_mb=20000,
		runtime_hrs=10,
		runtime_min=59
	wildcard_constraints:
		filter="unfiltered|lenient|strict",
		population="all-samples|unrelated-samples"
	shell:
		"zcat {input.vcf} | python3 workflow/scripts/select_ids.py {input.filters} {wildcards.filter} | bgzip -c > {output}"


rule test_hwe:
	input:
		vcf="results/population-typing/{source}/merged-vcfs/filtered/unrelated-samples_bi_all_{filter}.vcf.gz",
		bed= lambda wildcards: config['repeat_regions'][wildcards.source] if wildcards.region == "repeat" else "results/population-typing/{source}/bed/non-repetitve-regions.bed"
	output:
		hwe="results/population-typing/{source}/evaluation/hwe/unrelated-samples_{filter, unfiltered|lenient|strict}_{varianttype}-{region}.hwe",
		vcf="results/population-typing/{source}/evaluation/hwe/unrelated-samples_{filter, unfiltered|lenient|strict}_{varianttype}-{region}.vcf"
	params:
		prefix = "results/population-typing/{source}/evaluation/hwe/unrelated-samples_{filter, unfiltered|lenient|strict}_{varianttype}-{region}",
		bedtools = lambda wildcards: "-u" if wildcards.region == "repeat" else "-f 1.0 -u" 
	conda:
		'../envs/vcftools.yml'
	wildcard_constraints:
		filter="unfiltered|lenient|strict",
		region="repeat|nonrep",
		varianttype="|".join(['snp', 'indels', 'small-deletion', 'small-insertion', 'small-complex', 'midsize-deletion', 'midsize-insertion', 'large-complex', 'large-deletion', 'large-insertion', 'large-complex'])
	resources:
		mem_total_mb=30000,
		runtime_hrs=2,
		runtime_min=59
	shell:
		"""
		zcat {input.vcf} | python3 workflow/scripts/extract-varianttype.py {wildcards.varianttype} | bedtools intersect -header -a - -b {input.bed} {params.bedtools} > {output.vcf}	
		vcftools --vcf {output.vcf} --hardy --max-missing 0.9 --out {params.prefix}
		"""

# get all regions outside of repeat regions
rule get_non_repeat_regions:
	input:
		fai= lambda wildcards: config['reference'][wildcards.source] + '.fai',
		repeats= lambda wildcards: config['repeat_regions'][wildcards.source]
	output:
		subset_bed=temp("results/population-typing/{source}/bed/bed-tmp.bed"),
		fai=temp("results/population-typing/{source}/bed/fai-tmp.bed"),
		bed=temp("results/population-typing/{source}/bed/non-repetitve-regions.bed")
	conda:
		"../envs/genotyping.yml"
	shell:
		"bash workflow/scripts/non-repetitive-regions.sh {input.repeats} {output.subset_bed} {input.fai} {output.fai} {output.bed}"


def hwe_statistics_files(wildcards):
	files = []
	varis = [v for v in config['variants'][wildcards.source]]
	for var in varis:
		for reg in ["repeat", "nonrep"]:
			files.append("results/population-typing/{source}/evaluation/hwe/unrelated-samples_{filter}_{varianttype}-{region}.hwe".format(source=wildcards.source, filter=wildcards.filter, varianttype=var, region=reg))
	return files


def hwe_statistics_labels(wildcards):
	labels=[]
	varis = [v for v in config['variants'][wildcards.source]]
	for var in varis:
		for reg in ["repeat", "nonrep"]:
			labels.append(var + '-' + reg)
	return labels
	

rule compute_hwe_statistics:
	input:
		hwe_statistics_files
	output:
		tsv="results/population-typing/{source}/evaluation/statistics/{filter}/unrelated-samples_{filter}_hwe.tsv"
	log:
		"results/population-typing/{source}/evaluation/statistics/{filter}/unrelated-samples_{filter}.log"
	params:
		outname="results/population-typing/{source}/evaluation/statistics/{filter}/unrelated-samples_{filter}",
		labels=hwe_statistics_labels
	conda:
		'../envs/genotyping.yml'
	resources:
		mem_total_mb=20000,
		runtime_hrs=2
	shell:
		"python3 workflow/scripts/hwe.py {input} --labels {params.labels} --outname {params.outname} &> {log}"


###################################################
# compute mendelian consistency
###################################################

# count variants mendelian consistent in 0,1,2,...,nr_trios trios
rule check_consistent_trios:
	input:
		vcf="results/population-typing/{source}/merged-vcfs/whole-genome/all-samples_bi_all.vcf.gz",
		ped=trios,
		samples= lambda wildcards: config['pilot_samples'] if run_pilot else "results/population-typing/{source}/sample-index.tsv"
	output:
		variant_stats="results/population-typing/{source}/evaluation/statistics/all/mendelian-statistics_bi_all.tsv",
		trio_stats="results/population-typing/{source}/evaluation/statistics/all/trio-statistics_bi_all.tsv"
	log:
		"results/population-typing/{source}/evaluation/statistics/all/mendelian-statistics_bi_all.log"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=300000,
		runtime_hrs=96,
	#	runtime_hrs=23,
		runtime_min=59
	shell:
		"python3 workflow/scripts/mendelian-consistency.py statistics -vcf {input.vcf} -ped {input.ped} -samples {input.samples} -table {output.variant_stats} -column-prefix pangenie > {output.trio_stats}"


###################################################
# compute allele frequency/genotype statistics
###################################################

rule compute_statistics:
	input:
		vcf="results/population-typing/{source}/merged-vcfs/whole-genome/unrelated-samples_bi_all.vcf.gz",
		vcf_all="results/population-typing/{source}/merged-vcfs/whole-genome/all-samples_bi_all.vcf.gz",
		panel= lambda wildcards: config['biallelic_vcf'][wildcards.source]
	output:
		"results/population-typing/{source}/evaluation/statistics/all/genotyping-statistics_bi_all.tsv"
	conda:
		'../envs/genotyping.yml'
	resources:
		mem_total_mb=200000,
		runtime_hrs=96,
		runtime_min=59
	shell:
		"python3 workflow/scripts/collect-vcf-stats.py {input.panel} {input.vcf} {input.vcf_all} > {output}"



#########################################
# self-genotyping evaluation
#########################################

	
# genotyping concordance for each ID (over all samples)
rule genotype_concordance_variants:
	input:
		computed="results/population-typing/{source}/merged-vcfs/whole-genome/all-samples_bi_all.vcf.gz",
		true= lambda wildcards: config['biallelic_vcf'][wildcards.source]
	output:
		"results/population-typing/{source}/evaluation/statistics/all/self_bi_all_variant-stats.tsv"
	params:
		file_prefix="results/population-typing/{source}/evaluation/statistics/all/self_bi_all",
		column_prefix="pangenie_self-genotyping",
		# restrict to panel samples that are part of 1000G samples
		samples= lambda wildcards: ','.join( list(set(config['panel_samples'][wildcards.source]) & set(samples_1000G) ) )
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=500000,
		runtime_hrs=40,
		runtime_min=59
	log:
		"results/population-typing/{source}/evaluation/statistics/all/self_bi_all_variant-stats.log"
	shell:
		"python3 workflow/scripts/genotype-concordance-variant.py {input.true} {input.computed} {params.file_prefix} {params.samples} {params.column_prefix} &> {log}"


#################################################
# make a table containing var ID and bubble ID
#################################################

rule variant_id_to_bubble:
	input:
		lambda wildcards: config['biallelic_vcf'][wildcards.source]
	output:
		"results/population-typing/{source}/evaluation/statistics/all/bubble-statistics_bi_all.tsv"
	shell:
		"zcat {input} | python3 workflow/scripts/id_to_bubble.py > {output}"


#################################################
# annotate variants by BED file
#################################################

rule annotate_variants:
	input:
		vcf= lambda wildcards: config['biallelic_vcf'][wildcards.source],
		bed= lambda wildcards: config['repeat_regions'][wildcards.source]
	output:
		"results/population-typing/{source}/evaluation/statistics/all/annotations_bi_all.tsv"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=70000,
		runtime_hrs=1,
		runtime_min=59
	shell:
		"bedtools annotate -i {input.vcf} -files {input.bed} | python3 workflow/scripts/annotate_repeats.py -names repeats -format tsv > {output}"


#################################################
# plot the results
#################################################


rule merge_table:
	input:
		"results/population-typing/{source}/evaluation/statistics/all/genotyping-statistics_bi_all.tsv",
		"results/population-typing/{source}/evaluation/statistics/all/mendelian-statistics_bi_all.tsv",
		"results/population-typing/{source}/evaluation/statistics/all/self_bi_all_variant-stats.tsv",
		"results/population-typing/{source}/evaluation/statistics/all/bubble-statistics_bi_all.tsv",
		"results/population-typing/{source}/evaluation/statistics/all/annotations_bi_all.tsv"
	output:
		"results/population-typing/{source}/evaluation/statistics/all/summary_bi_all.tsv"
	conda:
		"../envs/plotting.yml"
	resources:
		mem_total_mb=50000,
		runtime_hrs=5,
		runtime_min=1
	shell:
		"python3 workflow/scripts/merge-tables.py {input} {output}"



rule plot_statistics:
	input:
		"results/population-typing/{source}/evaluation/statistics/all/summary_bi_all.tsv"
	output:
#		expand("results/population-typing/{{source}}/evaluation/statistics/plot_bi_all_{vartype}_{filter}_{region}.pdf", vartype=["snps", "indels", "large_insertions", "large_deletions", "large_complex"], filter=['unfiltered', 'strict'], region=["all-regions", "repeat-regions", "nonrepeat-regions"]),
#		expand("results/population-typing/{{source}}/evaluation/statistics/plot_bi_all_{vartype}_{filter}_{region}.pdf", vartype=["large_insertions", "large_deletions", "large_complex"], filter=['lenient_-0.5'], region=["all-regions", "repeat-regions", "nonrepeat-regions"]),
		"results/population-typing/{source}/evaluation/statistics/plot_bi_all_filters.tsv"
	params:
		outprefix="results/population-typing/{source}/evaluation/statistics/plot_bi_all"
	log:
		"results/population-typing/{source}/evaluation/statistics/plot_bi_all.log"
	conda:
		'../envs/plotting.yml'
	resources:
		mem_total_mb=200000,
		runtime_hrs=8,
		runtime_min=59
	shell:
		"python3 workflow/scripts/analysis.py {input} {params.outprefix} &> {log}"

