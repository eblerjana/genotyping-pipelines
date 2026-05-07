
rule merge_create_file_list:
	"""
	Create file containing all filenames to be merged.
	"""
	input:
		vcfs = expand("{{results}}/genotyping/pangenie/{sample}_pangenie_bi_genotyping.vcf.gz", sample = ILLUMINA.keys()),
		tbis = expand("{{results}}/genotyping/pangenie/{sample}_pangenie_bi_genotyping.vcf.gz.tbi", sample = ILLUMINA.keys())
	output:
		"{results}/genotyping/pangenie_bi_filelist.tsv"
	run:
		f = open(output[0], 'w')
		for name in input.vcfs:
			print(name, file=f)
		f.close()


rule merge_by_region:
	"""
	Merge vcfs by region.
	"""
	input:
		filelist = "{results}/genotyping/pangenie_bi_filelist.tsv"
	output:
		vcf = "{results}/genotyping/merged-region-wise/pangenie_{region}_genotypes.vcf.gz",
		tbi = "{results}/genotyping/merged-region-wise/pangenie_{region}_genotypes.vcf.gz.tbi"
	log:
		"{results}/genotyping/merged-region-wise/pangenie_{region}_genotypes.log"
	resources:
		mem_mb = lambda wildcards, attempt: 100000 * attempt,
		walltime = "15:00:00"
	wildcard_constraints:
		region = "chr[0-9A-Z]+:[0-9]+-[0-9]+"
	threads: 12
	benchmark:
		"{results}/genotyping/merged-region-wise/pangenie_{region}_genotypes.benchmark.txt"
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		bcftools merge -r {wildcards.region} --regions-overlap 0 -l {input} --threads {threads} | bcftools +fill-tags -Oz -o {output.vcf} -- -t AN,AC,AF,AC_Hom,AC_Het &> {log}
		tabix -p vcf {output.vcf}
		"""


rule merge_create_regions_filelist:
	"""
	Create file containing all region-wise VCFs to merge.
	"""
	input:
		vcfs = expand("{{results}}/genotyping/merged-region-wise/pangenie_{region}_genotypes.vcf.gz", region = MERGING_REGIONS)
	output:
		"{results}/genotyping/pangenie_regions_filelist.tsv"
	run:
		f = open(output[0], 'w')
		for name in input.vcfs:
			print(name, file=f)
		f.close()



rule merge_concat_regions:
	"""
	Concat region-wise VCF files.
	"""
	input:
		filelist =  "{results}/genotyping/pangenie_regions_filelist.tsv"
	output:
		"{results}/genotyping/pangenie_all-samples_unfiltered.vcf.gz"
	log:
		"{results}/genotyping/pangenie_all-samples_unfiltered.log"
	benchmark:
		"{results}/genotyping/pangenie_all-samples_unfiltered.benchmark.txt"
	resources:
		mem_mb=10000,
		walltime="04:00:00"
	conda:
		"../envs/genotyping.yml"
	threads:
		24
	shell:
		"""
		bcftools concat -o {output} -O z -f {input.filelist} --threads {threads} &> {log}
		tabix -p vcf {output}
		"""


rule merge_clean_up:
	input:
		merged = "{results}/genotyping/pangenie_all-samples_unfiltered.vcf.gz",
		vcf = "{results}/genotyping/merged-region-wise/pangenie_{region}_genotypes.vcf.gz",
		tbi = "{results}/genotyping/merged-region-wise/pangenie_{region}_genotypes.vcf.gz.tbi"
	output:
		touch("{results}/genotyping/merged-region-wise/pangenie_{region}_genotypes_deleted.txt")
	shell:
		"""
		rm {input.vcf}
		rm {input.tbi}
		"""
