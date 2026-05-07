

rule shapeit_extract_males:
	"""
	Collect names of all male samples.
	"""
	input:
		SAMPLE_SHEET
	output:
		"{results}/phasing/haploid-samples.txt"
	shell:
		"awk '$5==1' {input} > {output}"



rule shapeit_extract_chromosome:
	"""
	Extract specific chromosome and set low quality genotypes
	to missing.
	"""
	input:
		"{results}/filtering/{callset}.vcf.gz"
	output:
		"{results}/phasing/vcf/{callset}_{chrom}.vcf.gz"
	conda:
		"../envs/shapeit.yaml"
	benchmark:
		"{results}/phasing/vcf/{callset}_{chrom}.benchmark.txt"
	resources:
		mem_mb = 70000,
		walltime = "05:00:00"
	threads: 10
	shell:
		"""
		bcftools view -r {wildcards.chrom} {input} --threads {threads} | bgzip > {output}
		tabix -p vcf {output}
		"""

rule shapeit_prepare_trios:
	"""
	Prepare a file with trio information.
	"""
	input:
		SAMPLE_SHEET
	output:
		"{results}/phasing/trios.ped"
	shell:
		"""
		awk '($3!=\"0\") && ($4!=\"0\")' {input} | cut -f2,3,4 > {output}
		"""


rule shapeit_phase_common:
	"""
	Phase a chromosome using shapeit.
	"""
	input:
		vcf = "{results}/phasing/vcf/{callset}_{chrom}.vcf.gz",
		fam = "{results}/phasing/trios.ped",
		map = lambda wildcards: MAPS[wildcards.chrom],
		haploids = "{results}/phasing/haploid-samples.txt"
	output:
		"{results}/phasing/{callset}_shapeit_{chrom}.bcf"
	log:
		"{results}/phasing/{callset}_shapeit_{chrom}.log"
	benchmark:
		"{results}/phasing/{callset}_shapeit_{chrom}-benchmark.txt"
	conda:
		"../envs/shapeit.yaml"
	threads: 32
        resources:
		mem_mb = 500000, # 100000
		walltime = "30:00:00"
	params:
		haploids = lambda wildcards: "--haploids " + "{results}/phasing/haploid-samples.txt".format(results = wildcards.results)  if ("X" in wildcards.chrom) or ("Y" in wildcards.chrom) else ""
	shell:
		"""
		SHAPEIT5_phase_common --input {input.vcf} --pedigree {input.fam} --region {wildcards.chrom} {params.haploids} --map {input.map} --output {output} --thread {threads} &> {log}
		"""


rule shapeit_concat_vcfs:
	"""
	Combine the phased per-chromosome VCFs into a single one.
	"""
	input:
		expand("{{results}}/phasing/{{callset}}_shapeit_{chrom}.bcf", chrom = [c for c in MAPS.keys()])
	output:
		"{results}/phasing/{callset}_shapeit.bcf"
	conda:
		"../envs/shapeit.yaml"
	benchmark:
		"{results}/phasing/{callset}_shapeit.benchmark.txt"
	threads: 24
        resources:
		mem_mb = 100000,
		walltime = "05:00:00"
	log:
		"{results}/phasing/{callset}_shapeit.log"
	shell:
		"""
		bcftools concat -o {output} --threads {threads} {input} &> {log}
		bcftools index {output}
		"""

