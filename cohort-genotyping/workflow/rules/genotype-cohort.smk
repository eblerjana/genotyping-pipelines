
rule genotyping_inputs:
	"""
	Helper rule.
	"""
	input:
		PANEL_MULTI,
		PANEL_BI
	output:
		"{results}/genotyping/input-panel-ready.txt"
	shell:
		"""
		echo "Panels ready" > {output}
		"""


rule genotyping_pangenie_prepare_panel:
	"""
	Prepare (uncompressed) PanGenie input panel.
	"""
	input:
		vcf = PANEL_MULTI
	output:
		temp("{results}/genotyping/panel-multi.vcf")
	shell:
		"""
		gunzip -c {input} > {output}
		"""


rule genotyping_pangenie_index:
	"""
	Create index for PanGenie.
	"""
	input:
		vcf = "{results}/genotyping/panel-multi.vcf",
		fasta = REFERENCE,
	output:
		directory("{results}/genotyping/pangenie/index/")
	log:
		"{results}/genotyping/pangenie/index.log"
	resources:
		mem_mb = 200000,
		walltime = "7:00:00"
	threads: 24
	params:
		out_prefix = "{results}/genotyping/pangenie/index/index"
	benchmark:
		"{results}/genotyping/pangenie/index.benchmark.txt"
	conda:
		"../envs/pangenie.yml"
	shell:
		"""
		mkdir {output}
		PanGenie-index -v {input.vcf} -r {input.fasta} -o {params.out_prefix} -t {threads}  &> {log}
		"""


rule genotyping_pangenie_genotype_sampling:
	"""
	Run genotyping using sampling.
	"""
	input:
		reads = lambda wildcards: ILLUMINA[wildcards.sample],
		index = "{results}/genotyping/pangenie/index/"
	output:
		temp("{results}/genotyping/pangenie/{sample}_pangenie_multi_genotyping.vcf")
	log:
		"{results}/genotyping/pangenie/{sample}_pangenie_multi_genotyping.log"
	resources:
		mem_mb = 60000,
		walltime = "3:00:00"
	params:
		index = "{results}/genotyping/pangenie/index/index",
		out_prefix = "{results}/genotyping/pangenie/{sample}_pangenie_multi",
	benchmark:
		"{results}/genotyping/pangenie/{sample}_pangenie.benchmark.txt"
	threads:
		24
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		PanGenie -f {params.index} -i <(gunzip -c {input.reads}) -o {params.out_prefix} -t {threads} -j {threads} -s {wildcards.sample}  &> {log}
		"""


rule genotyping_convert_genotypes_to_biallelic:
	"""
	Convert genotyped VCF to biallelic representation.
	"""
	input:
		vcf = "{results}/genotyping/pangenie/{sample}_pangenie_multi_genotyping.vcf",
		biallelic = PANEL_BI
	output:
		bi = "{results}/genotyping/pangenie/{sample}_pangenie_bi_genotyping.vcf.gz",
		bi_tbi = "{results}/genotyping/pangenie/{sample}_pangenie_bi_genotyping.vcf.gz.tbi",
		temp = temp("{results}/genotyping/pangenie/{sample}_pangenie_bi_genotyping_tmp.vcf.gz"),
		temp_tbi = temp("{results}/genotyping/pangenie/{sample}_pangenie_bi_genotyping_tmp.vcf.gz.tbi")
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_mb=30000
	priority: 1
	shell:
		"""
		cat {input.vcf} | python3 workflow/scripts/convert-to-biallelic.py {input.biallelic} | bgzip > {output.temp}
		tabix -p vcf {output.temp}
		bcftools sort {output.temp} -Oz -o {output.bi}
		tabix -p vcf {output.bi}
		"""
