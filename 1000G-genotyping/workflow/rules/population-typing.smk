from collections import defaultdict

configfile: "config/config.yaml"

pangenie = config['pangenie']
read_dir_unrelated = config['read_dir_unrelated']
read_dir_related = config['read_dir_related']
sample_index_unrelated = config['sample_index_unrelated']
sample_index_related = config['sample_index_related']
run_pilot = config['run_pilot']
ped_file = config['pilot_trios'] if run_pilot else config['trios']

populations = ['EUR', 'AFR', 'EAS', 'SAS', 'AMR', 'all']

# read all samples and create a dict: sample_name -> [read_file1, read_file2]
samples = defaultdict(lambda: dict())

for source in config['graph_vcf'].keys():
	samples_to_consider = set([s for s in config['panel_samples'][source]])

	if run_pilot:
		for s in open(config['pilot_samples'], 'r'):
			samples_to_consider.add(s.strip())
		for other in config['graph_vcf'].keys():
			if other == source:
				continue
			for ps in config['panel_samples'][other]:
				samples_to_consider.add(ps)
		print('Running pilot for ' + source + ' on ' + str(len(samples_to_consider)) + ' samples.')

	# read 2504 unrelated samples
	for line in open(sample_index_unrelated, 'r'):
		if line.startswith('study'):
			continue
		fields = line.split('\t')
		run_id = fields[3]
		sample_name = fields[2]
		if run_pilot and not sample_name in samples_to_consider:
			continue
		samples[source][sample_name] = read_dir_unrelated + sample_name + '_' + run_id + '.fasta.gz'


	# read 689 unrelated samples
	for line in open(sample_index_related, 'r'):
		if line.startswith('study'):
			continue
		fields = line.split('\t')
		run_id = fields[3]
		sample_name = fields[2]
		if run_pilot and not sample_name in samples_to_consider:
			continue
		samples[source][sample_name] = read_dir_related + sample_name + '_' + run_id + '.fasta.gz'


	# add missing panel samples
	for sample,path in config['additional_read_data'].items():
		samples[source][sample] = path


regions = defaultdict(list)
chromosomes = defaultdict(list)
chromosome_lengths = defaultdict(lambda: defaultdict(int))

for source in config['graph_vcf'].keys():
	for line in open(config['reference'][source] + '.fai', 'r'):
		fields = line.strip().split()
		chrom = fields[0]
		if not chrom in [str(i) for i in range(1,23)] + ['X', 'Y'] and not chrom in ['chr' + str(i) for i in range(1,23)] + ['chrX', 'chrY']:
			continue
		chromosomes[source].append(chrom)
		chromosome_lengths[source][chrom] = int(fields[1])

step_size = 1000000
for source in config['graph_vcf'].keys():
	for chromosome in sorted(chromosomes[source]):
		pos = 1
		while pos < chromosome_lengths[source][chromosome]:
			end_pos = min(pos+step_size-1, chromosome_lengths[source][chromosome])
			regions[source].append('{}:{}-{}'.format(chromosome, pos, end_pos))
			pos = pos + step_size


##################################################
#####   genotype variants using PanGenie    ######
##################################################


# if there are variants with AF=0, remove them
rule prepare_population_panel:
	input:
		lambda wildcards: config['graph_vcf'][wildcards.source]
	output:
		"results/population-typing/{source}/panel/panel.vcf"
	conda:
		"../envs/genotyping.yml"
	log:
		 "results/population-typing/{source}/panel/panel.log"
	resources:
		mem_total_mb=20000
	shell:
		"bcftools view --min-ac 1 {input} 2> {log} 1> {output}"


# genotype variants with pangenie
rule genotyping:
	input:
		reads=lambda wildcards: samples[wildcards.source][wildcards.sample],
		reference=lambda wildcards: config['reference'][wildcards.source],
		panel = "results/population-typing/{source}/panel/panel.vcf"
	output:
		reads = temp("results/population-typing/{source}/genotyping/reads-{sample}.fa"),
		path_segments = temp("results/population-typing/{source}/genotyping/pangenie-{sample}_path_segments.fasta"),
		genotyping_vcf = "results/population-typing/{source}/genotyping/pangenie-{sample}_genotyping.vcf"
	log:
		"results/population-typing/{source}/genotyping/pangenie-{sample}.log"
	threads: 24
	resources:
		mem_total_mb=100000,
		runtime_hrs=5,
		runtime_min=1
	params:
		out_prefix="results/population-typing/{source}/genotyping/pangenie-{sample}"
	benchmark:
		"results/population-typing/{source}/genotyping/pangenie-{sample}_benchmark.txt"
	priority: 1
	shell:
		"""
		gunzip -c {input.reads} > {output.reads}
		module load Singularity
		(/usr/bin/time -v {pangenie} -i {output.reads} -v {input.panel} -r {input.reference} -o {params.out_prefix} -s {wildcards.sample} -j {threads} -t {threads} -g ) &> {log}
		"""


# compress multiallelic file and create version without SNVs
rule postprocess_multi:
	input:
		genotypes = "results/population-typing/{source}/genotyping/pangenie-{sample}_genotyping.vcf"
	output:
		output_multi_all="results/population-typing/{source}/genotyping/{sample}_genotyping_multi_all.vcf.gz",
		output_multi_nosnvs="results/population-typing/{source}/genotyping/{sample}_genotyping_multi_nosnvs.vcf.gz"
	log:
		"results/population-typing/{source}/genotyping/{sample}_genotyping_postprocess.log"
	resources:
		mem_total_mb=10000,
		runtime_hrs=0,
		runtime_min=20
	conda:
		"../envs/genotyping.yml"
	benchmark:
		"results/population-typing/{source}/genotyping/{sample}_genotyping_postprocess_benchmark.txt"
	shell:
		"(/usr/bin/time -v workflow/scripts/postprocess-multi.sh {input.genotypes} {output.output_multi_all} {output.output_multi_nosnvs}) &> {log}"


# tabix vcf file
rule tabix_vcf:
	input:
		"results/population-typing/{filename}.vcf.gz"
	output:
		"results/population-typing/{filename}.vcf.gz.tbi"
	resources:
		mem_total_mb=20000,
		runtime_hrs=2,
		runtime_min=59
	conda:
		"../envs/genotyping.yml"
	shell:
		"tabix -p vcf {input}"
		

ruleorder: create_biallelic_vcf > tabix_vcf


# prepare list of multiallelic files to be merged
rule create_file_list_multi:
	input:
		vcfs=lambda wildcards: expand("results/population-typing/{{source}}/genotyping/{sample}_genotyping_multi_all.vcf.gz", sample=sorted(samples[wildcards.source].keys())),
#		tbi=lambda wildcards: expand("results/population-typing/{{source}}/genotyping/{sample}_genotyping_multi_all.vcf.gz.tbi", sample=sorted(samples[wildcards.source].keys()))
	output:
		"results/population-typing/{source}/multi_all_filelist.tsv"
	run:
		f = open(output[0], 'w')
		for name in input.vcfs:
			print(name, file=f)
		f.close()


# merge single-sample vcfs into multi-sample vcf for a specific region
rule merge_vcfs_by_region_multi:
	input:
		"results/population-typing/{source}/multi_all_filelist.tsv"
	output:
		"results/population-typing/{source}/merged-vcfs/region-wise/pangenie_merged_multi_all_{region}.vcf.gz"
	log:
		"results/population-typing/{source}/merged-vcfs/region-wise/pangenie_merged_multi_all_{region}.log"
	benchmark:
		"results/population-typing/{source}/merged-vcfs/region-wise/pangenie_merged_multi_all_{region}_benchmark.txt"
	resources:
		mem_total_mb=50000,
#		mem_total_mb=800000,
		runtime_hrs=2,
		runtime_min=59
	wildcard_constraints:
		region="chr[0-9A-Z]+:[0-9]+-[0-9]+"
	conda:
		"../envs/whatshap.yml" # NOTE: bcftools.yml has a problem using +fill-tags plugin
	shell:
		"(/usr/bin/time -v bcftools merge -r {wildcards.region} -l {input} | bcftools +fill-tags -Oz -o {output} -- -t AN,AC,AF) &> {log} "


# prepare list of biallelic files to be merged
rule create_file_list_bi:
	input:
		vcfs=lambda wildcards: expand("results/population-typing/{{source}}/genotyping/{sample}_genotyping_bi_all.vcf.gz", sample=sorted(samples[wildcards.source].keys())),
#		tbi=lambda wildcards: expand("results/population-typing/{{source}}/genotyping/{sample}_genotyping_bi_all.vcf.gz.tbi", sample=sorted(samples[wildcards.source].keys()))
	output:
		"results/population-typing/{source}/bi_all_filelist.tsv"
	run:
		f = open(output[0], 'w')
		for name in input.vcfs:
			print(name, file=f)
		f.close()


# merge all biallelic vcfs into a multisample vcf file
rule merge_vcfs_by_region_bi:
	input:
		filelist="results/population-typing/{source}/bi_all_filelist.tsv",
	output:
		vcf="results/population-typing/{source}/merged-vcfs/region-wise/pangenie_merged_bi_all_{region}.vcf.gz"
	log:
		"results/population-typing/{source}/merged-vcfs/region-wise/pangenie_merged_bi_all_{region}.log"
	resources:
		mem_total_mb=200000,
#		mem_total_mb=1000000,
		runtime_hrs=7,
		runtime_min=59
	wildcard_constraints:
		region="chr[0-9A-Z]+:[0-9]+-[0-9]+"
	threads: 10
	benchmark:
		"results/population-typing/{source}/merged-vcfs/region-wise/pangenie_merged_bi_all_{region}_benchmark.txt"
	conda:
		"../envs/whatshap.yml" # NOTE: bcftools.yml has a problem using +fill-tags plugin
	shell:
		"(/usr/bin/time -v bcftools merge -r {wildcards.region} -l {input.filelist} --threads {threads} | bcftools +fill-tags -Oz -o {output.vcf} -- -t AN,AC,AF) &> {log} "


# combine all region-specific vcfs into one whole genome one
rule concat_vcfs_filelist:
	input:
		vcfs= lambda wildcards: expand("results/population-typing/{{source}}/merged-vcfs/region-wise/pangenie_merged_{{what}}_all_{region}.vcf.gz", region=regions[wildcards.source]),
	#	tbi= lambda wildcards: expand("results/population-typing/{{source}}/merged-vcfs/region-wise/pangenie_merged_{{what}}_all_{region}.vcf.gz.tbi", region=regions[wildcards.source])
	output:
		"results/population-typing/{source}/merged-vcfs/region-wise/concat_{what}_all_filelist.tsv"
	run:
		f = open(output[0], 'w')
		for name in input.vcfs:
			print(name, file=f)
		f.close()


# concat per-region vcfs
rule concat_vcfs:
	input:
		"results/population-typing/{source}/merged-vcfs/region-wise/concat_{what}_all_filelist.tsv"
	output:
		vcf="results/population-typing/{source}/merged-vcfs/whole-genome/all-samples_{what}_all.vcf.gz"
	log:
		"results/population-typing/{source}/merged-vcfs/whole-genome/all-samples_{what}_all.log"
	benchmark:
		"results/population-typing/{source}/merged-vcfs/whole-genome/all-samples_{what}_all_benchmark.txt"
	resources:
		mem_total_mb=50000,
		runtime_hrs=29,
		runtime_min=59
	conda:
		"../envs/genotyping.yml"
	shell:
		"bcftools concat -o {output.vcf} -O z -f {input} &> {log} "


# remove SNVs from vcf file
rule filter_snvs:
	input:
		vcf="results/population-typing/{source}/merged-vcfs/whole-genome/all-samples_{what}_all.vcf.gz",
		tbi="results/population-typing/{source}/merged-vcfs/whole-genome/all-samples_{what}_all.vcf.gz.tbi"
	output:
		"results/population-typing/{source}/merged-vcfs/whole-genome/all-samples_{what}_nosnvs.vcf.gz"
	log:
		"results/population-typing/{source}/merged-vcfs/whole-genome/all-samples_{what}_nosnvs.log"
	benchmark:
		"results/population-typing/{source}/merged-vcfs/whole-genome/all-samples_{what}_nosnvs_benchmark.txt"
	resources:
		mem_total_mb=50000,
		runtime_hrs=27,
		runtime_min=59
	conda:
		"../envs/genotyping.yml"
	shell:
		"bcftools view --exclude-types snps -o {output} -O z {input.vcf} &> {log} "


# convert multi-allelic vcf to biallelic representation
rule create_biallelic_vcf:
	input:
		template= lambda wildcards: config['biallelic_vcf'][wildcards.source],
		vcf="results/population-typing/{source}/genotyping/{sample}_genotyping_multi_all.vcf.gz"
	output:
		vcf_bi_all="results/population-typing/{source}/genotyping/{sample}_genotyping_bi_all.vcf.gz",
		tbi_bi_all="results/population-typing/{source}/genotyping/{sample}_genotyping_bi_all.vcf.gz.tbi",
		vcf_bi_nosnvs="results/population-typing/{source}/genotyping/{sample}_genotyping_bi_nosnvs.vcf.gz",
		tbi_bi_nosnvs="results/population-typing/{source}/genotyping/{sample}_genotyping_bi_nosnvs.vcf.gz.tbi"
	log:
		"results/population-typing/{source}/genotyping/{sample}_genotyping_postprocess_bi.log"
	benchmark:
		"results/population-typing/{source}/genotyping/{sample}_genotyping_postprocess_bi_benchmark.txt"
	resources:
		mem_total_mb=30000,
		runtime_hrs=1,
		runtime_min=59
	conda:
		"../envs/genotyping.yml"
	shell:
		"(/usr/bin/time -v workflow/scripts/postprocess-bi.sh {input.template} {input.vcf} {output.vcf_bi_all} {output.vcf_bi_nosnvs}) &> {log}"
