import gzip

configfile: "config.yaml"

nr_samples = {}

for callset in config['nr_samples'].keys():
	total = config['nr_samples'][callset]['total']
	males = config['nr_samples'][callset]['males']
	
	min_an = int(0.8 * total * 2)
	min_an_male = int(0.8 * ( ((total-males)*2) + males) )
	min_an_female = int(0.8 * males * 2)
	
	nr_samples[callset] = [min_an, min_an_male, min_an_female]

reference_to_ignore = "CHM13" if config['reference_version'] == "GRCh38" else "GRCh38"


chromosomes = [config['reference_prefix'] + str(i) for i in range(1,23)] + ['chrX']
vcfbub = config['vcfbub']

# threshold on bubble size used for vcfbub
threshold = 100000


# filter out LV>0 bubbles
rule filter_bubbles:
	input:
		lambda wildcards: config['vcf'][wildcards.caller]
	output:
		"results/vcf/{caller}/{caller}.vcf.gz"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=50000
	shell:
		"""
		{vcfbub} -l 0 {threshold} --input {input} | bgzip -c > {output}
		tabix -p vcf {output}
		"""

# remove sites that:
# - are covered by missing alleles in more than min_an haplotypes
# - are located outside of given chromosomes
# - contain Ns in sequence
rule filter_vcf:
	input:
		"results/vcf/{caller}/{caller}.vcf.gz"
	output:
		temp("results/vcf/{caller}/{caller}_filtered.vcf")
	conda:
		"../envs/genotyping.yml"
	log:
		"results/vcf/{caller}/{caller}_filtered.log"
	resources:
		mem_total_mb=20000,
		runtime_hrs=0,
		runtime_min=59
	params:
		min_an = lambda wildcards: nr_samples[wildcards.caller][0],
		min_an_male = lambda wildcards: nr_samples[wildcards.caller][1],
		min_an_female = lambda wildcards: nr_samples[wildcards.caller][2]
	shell:
		"bcftools view --samples ^{reference_to_ignore} {input} | bcftools view --min-ac 1 | python3 workflow/scripts/filter-vcf.py {min_an} {min_an_male} {min_an_female} --chromosomes {chromosomes} 2> {log} 1> {output}"


# remove alternative alleles that are not covered by any haplotype
rule trim_alt_alleles:
	input:
		"results/vcf/{caller}/{caller}_filtered.vcf"
	output:
		temp("results/vcf/{caller}/{caller}_filtered_trim.vcf")
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=20000,
		runtime_hrs=0,
		runtime_min=30
	shell:
		"bcftools view --trim-alt-alleles {input} > {output}"


# annotate the VCF and create an equivalent biallelic version of it
rule annotate_vcf:
	input:
		vcf="results/vcf/{caller}/{caller}_filtered_trim.vcf",
		gfa=lambda wildcards: config['gfa'][wildcards.caller]
	output:
		multi="results/vcf/{caller}/{caller}_filtered_ids.vcf",
		multi_tmp=temp("results/vcf/{caller}/{caller}_filtered_ids-tmp.vcf"),
		biallelic="results/vcf/{caller}/{caller}_filtered_ids_biallelic.vcf",
		bi_tmp=temp("results/vcf/{caller}/{caller}_filtered_ids-tmp_biallelic.vcf")
	log:
		"results/vcf/{caller}/{caller}_filtered_ids.log"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=200000,
		runtime_hrs=5,
		runtime_min=59
	params:
		outname = "results/vcf/{caller}/{caller}_filtered_ids-tmp"
	shell:
		"""
		python3 workflow/scripts/annotate_vcf.py -vcf {input.vcf} -gfa {input.gfa} -o {params.outname} &> {log}
		cat {output.multi_tmp} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n\"}}' > {output.multi}
		cat {output.bi_tmp} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n\"}}' > {output.biallelic}
		"""

# compress annotated
rule compress_annotated:
	input:
		"results/vcf/{caller}/{caller}_filtered_{which}.vcf"
	output:
		"results/vcf/{caller}/{caller}_filtered_{which}.vcf.gz"
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		bgzip -c {input} > {output}
		tabix -p vcf {output}
		"""
