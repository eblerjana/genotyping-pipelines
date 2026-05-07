callset_to_prefix = {
	"panel": "panel",
	"pangenie_all-samples_unfiltered": "pangenie-all",
	"pangenie_unrelated-samples_unfiltered": "pangenie-unrelated"
}

rule filter_collect_samples:
	"""
	Collect different sample sets.
	"""
	input:
		SAMPLE_SHEET
	output:
		unrelated = "{results}/filtering/unrelated-samples.tsv",
		all = "{results}/filtering/all-samples.tsv",
		panel = "{results}/filtering/overlap-panel-samples.tsv"
	run:
		with open(output.all, 'w') as all_samples, open(output.unrelated, 'w') as unrelated_samples, open(output.panel, 'w') as panel_samples, open(input[0], 'r') as infile:
			for line in infile:
				if line.startswith('#'):
					continue
				fields = line.strip().split()
				sample = fields[1]
				if fields[7] == "nan":
					continue
				all_samples.write(sample + '\n')
				if (fields[2] == '0') and (fields[3] == '0'):
					# unrelated sample
					unrelated_samples.write(sample + '\n')
				if sample in PANEL_SAMPLES:
					# overlap with panel samples
					panel_samples.write(sample + '\n')



rule filter_extract_unrelated_samples:
	"""
	Extract only unrelated (parent) samples.
	"""
	input:
		vcf = "{results}/genotyping/pangenie_all-samples_unfiltered.vcf.gz",
		samples = "{results}/filtering/unrelated-samples.tsv"
	output:
		"{results}/genotyping/pangenie_unrelated-samples_unfiltered.vcf.gz"
	benchmark:
		"{results}/genotyping/pangenie_unrelated-samples_unfiltered.benchmark.txt"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_mb = 80000,
		walltime = "12:00:00"
	shell:
		"""
		bcftools view --samples-file {input.samples} --force-samples {input.vcf} -Oz -o {output}
		tabix -p vcf {output}
		"""


rule filter_extract_chromosome:
	"""
	Extract a chromosome from the VCF
	"""
	input:
		lambda wildcards: PANEL_BI if wildcards.callset == "panel" else "{results}/genotyping/{callset}.vcf.gz"
	output:
		temp("{results}/filtering/chromosome-wise/{callset}_{chrom}.vcf.gz")
	benchmark:
		"{results}/filtering/chromosome-wise/{callset}_{chrom}.benchmark.txt"
	wildcard_constraints:
		callset = "panel|pangenie_all-samples_unfiltered|pangenie_unrelated-samples_unfiltered"
	resources:
		mem_mb = 50000,
		walltime = "02:00:00"
	conda:
		"../envs/genotyping.yml"
	params:
		tags = lambda wildcards: "AN,AC,AF,AC_Hom,AC_Het" if wildcards.callset == "panel" else "AN,AC,AF,AC_Hom,AC_Het,'HGQ:1=count(GQ>=200)'" 
	shell:
		"""
		bcftools view -r {wildcards.chrom} {input} | bcftools +fill-tags -Oz -o {output} -- -t {params.tags}
		"""


rule filter_compute_mendelian_consistency:
	"""
	Compute mendelian consistency.
	"""
	input:
		vcf = "{results}/filtering/chromosome-wise/pangenie_all-samples_unfiltered_{chrom}.vcf.gz",
		trios = SAMPLE_SHEET,
		samples = "{results}/filtering/all-samples.tsv"
	output:
		"{results}/filtering/pangenie_all-samples_unfiltered_mendelian-consistency_{chrom}.tsv"	
	conda:
		"../envs/genotyping.yml"
	log:
		"{results}/filtering/pangenie_all-samples_unfiltered_mendelian-consistency_{chrom}.log"
	benchmark:
		"{results}/filtering/pangenie_all-samples_unfiltered_mendelian-consistency_{chrom}.benchmark.tsv"
	resources:
		mem_mb = 30000,
		walltime = "08:00:00"
	shell:
		"""
		python3 workflow/scripts/evaluate-mendelian-consistency.py statistics -vcf {input.vcf} -ped {input.trios} -samples {input.samples} -table {output} -column-prefix pangenie &> {log}
		"""


rule filter_compute_genotype_statistics:
	"""
	Compute genotype statistics.
	"""
	input:
		"{results}/filtering/chromosome-wise/{callset}_{chrom}.vcf.gz"
	output:
		"{results}/filtering/{callset}_genotype-statistics_{chrom}.tsv"
	wildcard_constraints:
		callset = "panel|pangenie_all-samples_unfiltered|pangenie_unrelated-samples_unfiltered"
	benchmark:
		"{results}/filtering/{callset}_genotype-statistics_{chrom}.benchmark.txt"
	resources:
		mem_mb = 30000,
		walltime = "08:00:00"
	params:
		flag = lambda wildcards: "" if wildcards.callset == "panel" else "--genotyping-stats",
		column_prefix = lambda wildcards: callset_to_prefix[wildcards.callset]
	shell:
		"""
		zcat {input} | python3 workflow/scripts/collect-vcf-stats.py -outname {output} -column-prefix {params.column_prefix} {params.flag}
		"""

rule filter_self_genotyping_statistics:
	"""
	Evaluate genotyping of samples in panel.
	"""
	input:
		genotypes = "{results}/filtering/chromosome-wise/pangenie_all-samples_unfiltered_{chrom}.vcf.gz",
		panel = "{results}/filtering/chromosome-wise/panel_{chrom}.vcf.gz",
		samples = "{results}/filtering/overlap-panel-samples.tsv"
	output:
		"{results}/filtering/pangenie_all-samples_unfiltered_self-genotyping_{chrom}.tsv"
	log:
		"{results}/filtering/pangenie_all-samples_unfiltered_self-genotyping_{chrom}.log"
	benchmark:
		"{results}/filtering/pangenie_all-samples_unfiltered_self-genotyping_{chrom}.log"
	resources:
		mem_mb = 60000,
		walltime = "10:00:00",
	shell:
		"""
		python3 workflow/scripts/evaluate-self-genotyping.py {input.panel} {input.genotypes} {output} {input.samples} pangenie_self-genotyping &> {log}
		"""

rule filter_concat_tables:
	"""
	Concat chromosome-wise tables.
	"""
	input:
		expand("{{results}}/filtering/{{table}}_{chrom}.tsv", chrom = CHROMOSOMES)
	output:
		"{results}/filtering/{table}.tsv"
	wildcard_constraints:
		table = "pangenie_all-samples_unfiltered_mendelian-consistency|panel_genotype-statistics|pangenie_all-samples_unfiltered_genotype-statistics|pangenie_unrelated-samples_unfiltered_genotype-statistics|pangenie_all-samples_unfiltered_self-genotyping"
	shell:
		"""
		head -n 1 {input[0]} > {output}; tail -n +2 -q {input} >> {output}
		"""


rule filter_variant_id_to_bubble:
	"""
	Map variant IDs to bubbles.
	"""
	input:
		PANEL_MULTI
	output:
		"{results}/filtering/panel_bubble-statistics.tsv"
	benchmark:
		"{results}/filtering/panel_bubble-statistics.benchmark.txt"
	shell:
		"zcat {input} | python3 workflow/scripts/id_to_bubble.py > {output}"




rule filter_merge_tables:
	input:
		"{results}/filtering/pangenie_all-samples_unfiltered_mendelian-consistency.tsv",
		"{results}/filtering/panel_genotype-statistics.tsv",
		"{results}/filtering/pangenie_all-samples_unfiltered_genotype-statistics.tsv",
		"{results}/filtering/pangenie_unrelated-samples_unfiltered_genotype-statistics.tsv",
		"{results}/filtering/pangenie_all-samples_unfiltered_self-genotyping.tsv",
		"{results}/filtering/panel_bubble-statistics.tsv"
	output:
		"{results}/filtering/all_genotyping_statistics.tsv"
	benchmark:
		"{results}/filtering/all_genotyping_statistics.benchmark.txt"
	conda:
		"../envs/plotting.yml"
	resources:
		mem_mb = 150000,
		walltime = "05:00:00"
	shell:
		"""
		python3 workflow/scripts/merge-tables.py {input} {output}
		"""



rule filter_perform_regression:
	input:
		"{results}/filtering/all_genotyping_statistics.tsv"
	output:
		filters = "{results}/filtering/pangenie_filters.tsv",
		regression = "{results}/filtering/pangenie_regression.tsv"
	params:
		outprefix = "{results}/filtering/pangenie",
		threshold= 5 if len(ILLUMINA.keys()) < 1000 else 20
	log:
		"{results}/filtering/pangenie_regression.log"
	benchmark:
		"{results}/filtering/pangenie_regression.benchmark.txt"
	conda:
		'../envs/plotting.yml'
	resources:
		mem_mb=200000,
		walltime = "12:00:00"
	shell:
		"python3 workflow/scripts/regression.py -t {input} -o {params.outprefix} -n {params.threshold} --regression-only &> {log}"


rule filter_plot_statistics:
	input:
		"{results}/filtering/pangenie_regression.tsv"
	output:
		"{results}/filtering/pangenie_plot.log"
	benchmark:
		"{results}/filtering/pangenie_plot.benchmark.txt"
	params:
		outprefix="{results}/filtering/pangenie",
		threshold= 5 if len(ILLUMINA.keys()) < 1000 else 20
	conda:
		'../envs/plotting.yml'
	resources:
		mem_mb = 200000,
		walltime = "05:00:00"
	shell:
		"python3 workflow/scripts/regression.py -t {input} -o {params.outprefix} -n {params.threshold} --plot-only &> {output}"



rule filter_final_callsets:
	"""
	Create VCF with filtered variants.
	"""
	input:
		vcf = "{results}/genotyping/pangenie_all-samples_unfiltered.vcf.gz",
		filters = "{results}/filtering/pangenie_filters.tsv",
		fai = REFERENCE + ".fai"
	output:
		tmp = temp("{results}/filtering/tmp-pangenie_all-samples_filtered.vcf.gz"),
		vcf = "{results}/filtering/pangenie_all-samples_filtered.vcf.gz"
	benchmark:
		"{results}/filtering/pangenie_all-samples_filtered.benchmark.txt"
	resources:
		mem_mb = 20000,
		walltime = "12:00:00"
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		zcat {input.vcf} | python3 workflow/scripts/select_ids.py {input.filters} lenient | bgzip -c > {output.tmp} 
		bcftools reheader --fai {input.fai} {output.tmp} > {output.vcf}
		tabix -p vcf {output.vcf}
		"""
