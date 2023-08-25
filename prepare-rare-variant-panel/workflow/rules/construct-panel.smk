configfile: "config/config.yaml"

# remove overlapping records and prepare vcf
rule prepare_vcf:
	input:
		config["rare_variants_vcf"]
	output:
		temp("results/rare-alleles.vcf.gz")
	shell:
		"zcat {input} | workflow/scripts/prepare-vcf.py | bgzip > {output}"


# remove variants overlapping with panel variants and add remaining variants to multi-allelic panel
rule construct_panel:
	input:
		panel = config["pangenie_multi_vcf"],
		rare = "results/rare-alleles.vcf.gz"
	output:
		"results/panel-with-rare_multi.vcf"
	resources:
		mem_total_mb=20000
	conda:
		"../envs/rare.yml"
	shell:
		"""
		bedtools subtract -A -a {input.rare} -b {input.panel} | bcftools merge - {panel} > {output}
		"""
		
# remove variants overlapping with panel variants and add remaining variants to bi-allelic callset
rule construct_panel_biallelic:
	input:
		panel_bi = config["pangenie_bi_vcf"],
		rare = "results/rare-alleles.vcf.gz"
	output:
		"results/panel-with-rare_bi.vcf"
	resources:
		mem_total_mb=20000
	conda:
		"../envs/rare.yml"
	shell:
		"""
		bcftools norm -m- {rare} | bedtools subtract -A -a - -b {input.panel_bi} | bcftools merge - {panel_bi} > {output}
		"""


# write a report about the number of rare alleles added # TODO: write script
rule compute_statistics:
	input:
		"results/panel-with-rare_{mode}.vcf"
	output:
		"results/statistics_{mode}.txt"
	shell:
		"python3 workflow/scripts/compute-statistics.py {input} > {output}"
