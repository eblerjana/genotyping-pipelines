

rule mc_correct_sex_chromosomes:
	"""
	Correct genotypes on human sex chromosomes (X + Y).
	Assumes chromosomes are named as chrX, chrY or X,Y.
	"""
	input:
		vcf = MC_VCF,
		sample_info = SAMPLE_SHEET
	output:
		temp("{results}/mc-vcf/mc_corrected-sex-chromosomes.vcf.gz")
	log:
		"{results}/mc-vcf/mc_corrected-sex-chromosomes.log"
	benchmark:
		"{results}/mc-vcf/mc_corrected-sex-chromosomes.benchmark.txt"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_mb=20000,
		walltime="01:59:00"
	shell:
		"""
		python3 workflow/scripts/correct-sex-chromosomes.py {input.vcf} {input.sample_info} 2> {log} | bcftools +fill-tags -Oz -o {output} -- -t AN,AC,AF
		"""


rule mc_filter_vcf:
	"""
	Remove sites that:
	- are covered by missing alleles in more than MIN_FRAC haplotypes
	- contain Ns in their sequences
	"""
	input:
		"{results}/mc-vcf/mc_corrected-sex-chromosomes.vcf.gz"
	output:
		temp("{results}/mc-vcf/mc_filtered.vcf")
	conda:
		"../envs/genotyping.yml"
	log:
		"{results}/mc-vcf/mc_filtered.log"
	benchmark:
		"{results}/mc-vcf/mc_filtered.benchmark.txt"
	resources:
		mem_mb=20000,
		walltime="01:59:00"
	params:
		exclude = ','.join(SAMPLES_TO_EXCLUDE),
		min_frac = MIN_FRAC
	shell:
		"bcftools view --samples ^{params.exclude} --force-samples  {input} | bcftools view --min-ac 1 | python3 workflow/scripts/filter-vcf.py {params.min_frac} 2> {log} 1> {output}"


rule mc_trim_alt_alleles:
	"""
	Remove alternative alleles that are not covered by any haplotype.
	"""
	input:
		"{results}/mc-vcf/mc_filtered.vcf"
	output:
		temp("{results}/mc-vcf/mc_filtered_trim.vcf")
	benchmark:
		"{results}/mc-vcf/mc_filtered_trim.benchmark.txt"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_mb=20000,
		walltime="00:30:00"
	shell:
		"bcftools view --trim-alt-alleles {input} > {output}"


rule mc_annotate_vcf:
	"""
	Decompose bubbles, annotate VCF with resulting variant IDs
	and create equivalent biallelic VCF encoding decomposed alleles.
	"""
	input:
		vcf = "{results}/mc-vcf/mc_filtered_trim.vcf",
		gfa = MC_GFA
	output:
		multi = "{results}/mc-vcf/mc_filtered_ids.vcf.gz",
		multi_tmp = temp("{results}/mc-vcf/mc_filtered_ids-tmp.vcf"),
		biallelic = "{results}/mc-vcf/mc_filtered_ids_biallelic.vcf.gz",
		bi_tmp = temp("{results}/mc-vcf/mc_filtered_ids-tmp_biallelic.vcf")
	log:
		"{results}/mc-vcf/mc_decomposition.log"
	benchmark:
		 "{results}/mc-vcf/mc_decomposition.benchmark.txt"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_mb=200000,
		walltime="08:59:00"
	params:
		outname = "{results}/mc-vcf/mc_filtered_ids-tmp"
	shell:
		"""
		python3 workflow/scripts/annotate_vcf.py -vcf {input.vcf} -gfa {input.gfa} -o {params.outname} &> {log}
		cat {output.multi_tmp} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n\"}}' | bgzip > {output.multi}
		cat {output.bi_tmp} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n\"}}' | bgzip > {output.biallelic}
		tabix -p vcf {output.multi}
		tabix -p vcf {output.biallelic}
		"""

rule mc_norm_biallelic:
	"""
	Normalize decomposed VCF.
	"""
	input:
		biallelic = "{results}/mc-vcf/mc_filtered_ids_biallelic.vcf.gz",
		reference = REFERENCE,
	output:
		"{results}/mc-vcf/mc_filtered_ids_biallelic_norm.vcf.gz"
	log:
		"{results}/mc-vcf/mc_normalization.log"
	benchmark:
		"{results}/mc-vcf/mc_normalization.benchmark.txt"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_mb = 80000,
		walltime = "06:00:00"
	shell:
		"""
		bcftools norm -f {input.reference} {input.biallelic} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n\"}}'| bgzip > {output}
		tabix -p vcf {output}
		"""
