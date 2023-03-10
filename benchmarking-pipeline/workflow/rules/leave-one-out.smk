configfile: "config/config.yaml"

# stores paths to reads
reads_leave_one_out = {}

for line in open(config['reads'], 'r'):
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	sample_name = fields[1]
	read_path = fields[7]
	reads_leave_one_out[sample_name] = read_path

allowed_variants = ['snp', 'indels', 'large-deletion', 'large-insertion', 'large-complex']
callsets_leave_one_out = [s for s in config['callsets'].keys()]
coverages_leave_one_out = ['full'] + [c for c in config['downsampling']]
versions_leave_one_out = [v for v  in config['pangenie'].keys()]


################################################################
######   prepare input panel and ground truth genotypes  #######
################################################################

# remove positions that are ".|." in the left out sample. These cannot be used for evaluation, as the true genotype is unknown
# for now, also remove CHM13, because PanGenie cannot handle haploid samples
# remove sample from the panel
rule remove_missing:
	input:
		lambda wildcards: config['callsets'][wildcards.callset]['multi'] if wildcards.representation == 'multi' else config['callsets'][wildcards.callset]['bi']
	output:
		temp("results/leave-one-out/{callset}/preprocessed-vcfs/{sample}_{callset}_{representation}_no-missing.vcf")
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=20000
	priority: 1
	wildcard_constraints:
		representation = "bi|multi"
	shell:
		"zcat {input} | python3 workflow/scripts/remove-missing.py {wildcards.sample} > {output}"


rule prepare_panel:
	input:
		"results/leave-one-out/{callset}/preprocessed-vcfs/{sample}_{callset}_multi_no-missing.vcf.gz"
	output:
		"results/leave-one-out/{callset}/input-panel/panel-{sample}_{callset}.vcf"
	conda:
		"../envs/genotyping.yml"
	priority: 1
	log:
		"results/leave-one-out/{callset}/input-panel/panel-{sample}_{callset}.log"
	resources:
		mem_total_mb=20000
	shell:
		"bcftools view --samples ^{wildcards.sample} {input} | bcftools view --min-ac 1 2> {log} 1> {output}"



# extract ground truth genotypes for sample
rule prepare_truth:
	input:
		"results/leave-one-out/{callset}/preprocessed-vcfs/{sample}_{callset}_bi_no-missing.vcf.gz"
	output:
		"results/leave-one-out/{callset}/truth/truth-{sample}_{callset}.vcf"
	conda:
		"../envs/genotyping.yml"
	priority: 1
	resources:
		mem_total_mb=20000
	log:
		"results/leave-one-out/{callset}/truth/truth-{sample}_{callset}.log"
	shell:
		"bcftools view --samples {wildcards.sample} {input} 2> {log} 1> {output}"



rule compress_vcf:
	input:
		"results/leave-one-out/{filename}.vcf"
	output:
		vcf = "results/leave-one-out/{filename}.vcf.gz",
		tbi = "results/leave-one-out/{filename}.vcf.gz.tbi"
	priority: 1
	shell:
		"""
		bgzip -c {input} > {output.vcf}
		tabix -p vcf {output.vcf}
		"""


########################################################
##################    run PanGenie    ##################
########################################################


# run pangenie
rule pangenie:
	input:
		reads = lambda wildcards: reads_leave_one_out[wildcards.sample] if wildcards.coverage == 'full' else "results/downsampling/{callset}/{coverage}/{sample}_{coverage}.fa.gz",
		fasta = lambda wildcards: config['callsets'][wildcards.callset]['reference'],
		vcf="results/leave-one-out/{callset}/input-panel/panel-{sample}_{callset}.vcf"
	output:
		reads = temp("results/leave-one-out/{callset}/{version}/{sample}/{coverage}/reads.fa"),
		path_segments = temp("results/leave-one-out/{callset}/{version}/{sample}/{coverage}/pangenie-{sample}_path_segments.fasta"),
		genotyping = temp("results/leave-one-out/{callset}/{version}/{sample}/{coverage}/pangenie-{sample}_genotyping.vcf")
	log:
		"results/leave-one-out/{callset}/{version}/{sample}/{coverage}/pangenie-{sample}.log"
	threads: 24
	resources:
		mem_total_mb=190000,
		runtime_hrs=5,
		runtime_min=1
	priority: 1
	params:
		out_prefix="results/leave-one-out/{callset}/{version}/{sample}/{coverage}/pangenie-{sample}",
		pangenie = lambda wildcards: config['pangenie'][wildcards.version]
	shell:
		"""
		gunzip -c {input.reads} > {output.reads}
		module load Singularity
		(/usr/bin/time -v {params.pangenie} -i {output.reads} -v {input.vcf} -r /hilbert{input.fasta} -o {params.out_prefix} -s {wildcards.sample} -j {threads} -t {threads} -g ) &> {log}
		"""


########################################################
##################    Evaluation      ##################
########################################################


rule alleles_per_bubble:
	input:
		lambda wildcards: config['callsets'][wildcards.callset]['multi']
	output:
		plot = "results/leave-one-out/{callset}/alleles-per-bubble.pdf",
		bed = "results/leave-one-out/{callset}/complex-bubbles.bed"
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
		bed = "results/leave-one-out/{callset}/complex-bubbles.bed",
		fai = lambda wildcards: config['callsets'][wildcards.callset]['reference'] + '.fai'
	output:
		bed = "results/leave-one-out/{callset}/biallelic-bubbles.bed",
		tmp = temp("results/leave-one-out/{callset}/biallelic-bubbles.fai")
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
		vcf = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/pangenie-{sample}_genotyping.vcf.gz",
		biallelic = lambda wildcards: config['callsets'][wildcards.callset]['bi']
	output:
		"results/leave-one-out/{callset}/{version}/{sample}/{coverage}/pangenie-{sample}_genotyping-biallelic.vcf"
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
		lambda wildcards: config['callsets'][wildcards.callset]['bi']
	output:
		lists = "results/leave-one-out/{callset}/untypable-ids/{sample}-untypable.tsv",
		summary = "results/leave-one-out/{callset}/untypable-ids-{sample}.tsv"
	params:
		out = "results/leave-one-out/{callset}/untypable-ids"
	conda:
		"../envs/genotyping.yml"
	shell:
		"zcat {input} | python3 workflow/scripts/untypable-ids-single.py {params.out} {wildcards.sample} > {output.summary}"


# determine untypable IDs
rule remove_untypable:
	input:
		vcf = "results/leave-one-out/{callset}/{path}{sample}{other}.vcf",
		ids = "results/leave-one-out/{callset}/untypable-ids/{sample}-untypable.tsv"
	output:
		vcf = "results/leave-one-out/{callset}/{path}{sample}{other}-typable-{vartype}.vcf.gz",
		tbi = "results/leave-one-out/{callset}/{path}{sample}{other}-typable-{vartype}.vcf.gz.tbi"
	wildcard_constraints:
		callset = "|".join([v for v in config['callsets'].keys()]),
		sample = "|".join([s for s in reads_leave_one_out.keys()]),
		vartype = "|".join(allowed_variants)
	resources:
		mem_total_mb = 20000,
		runtime_hrs = 1
	priority: 1
	shell:
		"""
		cat {input.vcf} | python3 workflow/scripts/skip-untypable.py {input.ids} | python3 workflow/scripts/extract-varianttype.py {wildcards.vartype} | bgzip -c > {output.vcf}
		tabix -p vcf {output.vcf}
		"""

rule rtg_format:
	input:
		lambda wildcards: config['callsets'][wildcards.callset]['reference']
	output:
		directory("results/leave-one-out/{callset}/SDF")
	resources:
		mem_total_mb=20000
	priority: 1
	shell:
		'rtg format -o {output} {input}'


def region_to_bed(wildcards):
	if wildcards.regions == "biallelic":
		return "results/leave-one-out/{callset}/biallelic-bubbles.bed".format(callset=wildcards.callset)
	if wildcards.regions == "multiallelic":
		return "results/leave-one-out/{callset}/complex-bubbles.bed".format(callset=wildcards.callset)
	assert(False)


# precision-recall
rule vcfeval:
	input:
		callset = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/pangenie-{sample}_genotyping-biallelic-typable-{vartype}.vcf.gz",
		callset_tbi = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/pangenie-{sample}_genotyping-biallelic-typable-{vartype}.vcf.gz.tbi",
		baseline = "results/leave-one-out/{callset}/truth/truth-{sample}_{callset}-typable-{vartype}.vcf.gz",
		baseline_tbi = "results/leave-one-out/{callset}/truth/truth-{sample}_{callset}-typable-{vartype}.vcf.gz.tbi",
		regions = region_to_bed,
		sdf = "results/leave-one-out/{callset}/SDF"
	output:
		summary = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/precision-recall-typable/{regions}_{vartype}/summary.txt"
	conda:
		"../envs/genotyping.yml"
	priority: 1
	wildcard_constraints:
		sample = "|".join([s for s in reads_leave_one_out.keys()]),
		regions = "biallelic|multiallelic",
		vartype = "|".join(allowed_variants)
	params:
		tmp = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/precision-recall-typable/{regions}_{vartype}_tmp",
		outname = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/precision-recall-typable/{regions}_{vartype}",
		which = "--all-records"
	resources:
		mem_total_mb = 20000,
		runtime_hrs = 1,
		runtime_min = 40
	shell:
		"""
		rtg vcfeval -b {input.baseline} -c {input.callset} -t {input.sdf} -o {params.tmp} --ref-overlap --evaluation-regions {input.regions} {params.which} --Xmax-length 30000 > {output.summary}.tmp
		mv {params.tmp}/* {params.outname}/
		mv {output.summary}.tmp {output.summary}
		rm -r {params.tmp}
		"""


# determine the variants that went into re-typing per category
rule collect_typed_variants:
	input:
		callset = "results/leave-one-out/{callset}/preprocessed-vcfs/{sample}_{callset}_bi_no-missing.vcf.gz",
		regions= region_to_bed,
		ids="results/leave-one-out/{callset}/untypable-ids/{sample}-untypable.tsv"
	output:
		"results/leave-one-out/{callset}/genotyped-ids/{sample}_{regions}_{vartype}.tsv"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		sample = "|".join([s for s in reads_leave_one_out.keys()]),
		regions = "biallelic|multiallelic",
		vartype = "|".join(allowed_variants)
	resources:
		mem_total_mb=50000
	priority: 1
	shell:
		"zcat {input.callset} | python3 workflow/scripts/skip-untypable.py {input.ids} | python3 workflow/scripts/extract-varianttype.py {wildcards.vartype} | bedtools intersect -header -a - -b {input.regions} -u -f 0.5 | python3 workflow/scripts/get_ids.py > {output}"


# compute concordances
rule genotype_concordances:
	input:
		callset = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/pangenie-{sample}_genotyping-biallelic-typable-{vartype}.vcf.gz",
		callset_tbi = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/pangenie-{sample}_genotyping-biallelic-typable-{vartype}.vcf.gz.tbi",
		baseline = "results/leave-one-out/{callset}/truth/truth-{sample}_{callset}-typable-{vartype}.vcf.gz",
		baseline_tbi = "results/leave-one-out/{callset}/truth/truth-{sample}_{callset}-typable-{vartype}.vcf.gz.tbi",
		regions = region_to_bed,
		typed_ids = "results/leave-one-out/{callset}/genotyped-ids/{sample}_{regions}_{vartype}.tsv"
	output:
		tmp_vcf1 = temp("results/leave-one-out/{callset}/{version}/{sample}/{coverage}/concordance/{regions}_{vartype}_base.vcf"),
		tmp_vcf2 = temp("results/leave-one-out/{callset}/{version}/{sample}/{coverage}/concordance/{regions}_{vartype}_call.vcf"),
		summary = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/concordance/{regions}_{vartype}/summary.txt"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		sample = "|".join([s for s in reads_leave_one_out.keys()]),
		regions = "biallelic|multiallelic",
		vartype = "|".join(allowed_variants)
	log:
		"results/leave-one-out/{callset}/{version}/{sample}/{coverage}/concordance/{regions}_{vartype}/summary.log"
	resources:
		mem_total_mb = 40000,
		runtime_hrs = 0,
		runtime_min = 40
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
		lambda wildcards: expand("results/leave-one-out/{{callset}}/{{version}}/{sample}/{{coverage}}/{{metric}}/{{regions}}_{{vartype}}/summary.txt", sample = config['callsets'][wildcards.callset]['leave_one_out_samples'])
	output:
		"results/leave-one-out/{callset}/{version}/plots/{coverage}/{metric}_{callset}-{version}-{coverage}_{regions}_{vartype}.tsv"
	params:
		samples = lambda wildcards: ','.join([c for c in config['callsets'][wildcards.callset]['leave_one_out_samples']]),
		outfile = "results/leave-one-out/{callset}/{version}/plots/{coverage}/{metric}_{callset}-{version}-{coverage}_{regions}_{vartype}",
		folder = "results/leave-one-out/{callset}/{version}"
	priority: 1
	shell:
		"python3 workflow/scripts/collect-results.py {wildcards.metric} {wildcards.coverage} {params.samples} {wildcards.regions} -variants {wildcards.vartype} -folder {params.folder} -outfile {params.outfile}"



# plot results of different subsampling runs
rule plotting_versions:
	input:
		lambda wildcards: expand("results/leave-one-out/{{callset}}/{version}/plots/{{coverage}}/{m}_{{callset}}-{version}-{{coverage}}_{{regions}}_{vartype}.tsv", m = wildcards.metric if wildcards.metric != 'untyped' else 'concordance', version=versions_leave_one_out, vartype=config['callsets'][wildcards.callset]['variants'])
	output:
		"results/leave-one-out/{callset}/plots/comparison-versions/{metric}/{metric}_{coverage}_{regions}.pdf"
	wildcard_constraints:
		metric="concordance|precision-recall-typable|untyped",
		regions="biallelic|multiallelic"
	priority: 1
	conda:
		"../envs/genotyping.yml"
	params:
		sources = lambda wildcards: ' '.join([wildcards.callset + '-' + v + '-' + wildcards.coverage + '_' + wildcards.regions for v in versions_leave_one_out])
	shell:
		"python3 workflow/scripts/plot-results.py -files {input} -outname {output} -sources {params.sources} -metric {wildcards.metric}"



# plot results of different coverages
rule plotting_coverages:
	input:
		lambda wildcards: expand("results/leave-one-out/{{callset}}/{{version}}/plots/{coverage}/{m}_{{callset}}-{{version}}-{coverage}_{{regions}}_{vartype}.tsv", m = wildcards.metric if wildcards.metric != 'untyped' else 'concordance', coverage=coverages_leave_one_out, vartype=config['callsets'][wildcards.callset]['variants'])
	output:
		"results/leave-one-out/{callset}/plots/comparison-coverages/{metric}/{metric}_{version}_{regions}.pdf"
	wildcard_constraints:
		metric="concordance|precision-recall-typable|untyped",
		regions="biallelic|multiallelic"
	priority: 1
	conda:
		"../envs/genotyping.yml"
	params:
		sources = lambda wildcards: ' '.join([wildcards.callset + '-' + wildcards.version + '-' + c + '_' + wildcards.regions for c in coverages_leave_one_out])
	shell:
		"python3 workflow/scripts/plot-results.py -files {input} -outname {output} -sources {params.sources} -metric {wildcards.metric}"



# plot resources (single core CPU time and max RSS) for different subsampling runs
rule plotting_resources:
	input:
		lambda wildcards: expand("results/leave-one-out/{{callset}}/{version}/{sample}/{{coverage}}/pangenie-{sample}.log", version = versions_leave_one_out, sample = config['callsets'][wildcards.callset]['leave_one_out_samples'])
	output:
		"results/leave-one-out/{callset}/plots/resources/resources_{callset}-{coverage}.pdf"
	conda:
		"../envs/genotyping.yml"
	params:
		outname = "results/leave-one-out/{callset}/plots/resources/resources_{callset}-{coverage}",
		samples	= lambda wildcards: " ".join(config['callsets'][wildcards.callset]['leave_one_out_samples']),
		versions = " ".join(versions_leave_one_out)
	shell:
		"python3 workflow/scripts/plot-resources.py -files {input} -outname {params.outname} -samples {params.samples} -sizes {params.versions}"
