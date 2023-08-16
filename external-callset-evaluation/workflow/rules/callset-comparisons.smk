configfile: "config/config.yaml"

truthsets_sv = [t for t in config["truthsets"].keys() if config["truthsets"][t]["vartype"] == "sv"]
truthsets_small = [t for t in config["truthsets"].keys() if config["truthsets"][t]["vartype"] == "snp-indel"]

truthsets = truthsets_sv + truthsets_small
callsets = []
for t in truthsets:
	for c in config["truthsets"][t]["callsets"].keys():
		callsets.append(c)
regions = ['all'] + [r for r in config['regions_to_bed'].keys()]


# determine list of samples present in panel vcf
assembly_samples = []
for line in open(config["panel_vcf"], 'r'):
	if line.startswith("##"):
		continue
	if line.startswith("#"):
		assembly_samples = line.strip().split()[9:]


####################################################################################################
#  find out which variants in the truth set are not contained in the pangenome graph (=untypables)
####################################################################################################


# assign each variant a unique ID (if a variant matches with panel, use the same ID as panel). Only keep
# variants on chromosomes 1-22, X and Y.
rule annotate_variants_callset:
	input:
		callset= lambda wildcards: config["truthsets"][wildcards.truthset]["path"],
		panel= config["panel_vcf"],
		ref_index=config["reference_fai"]
	output:
		"results/data/vcf/{truthset}/truthset-{truthset}.vcf.gz",
	wildcard_constraints:
		truthset = "|".join(truthsets)
	resources:
		mem_total_mb=20000,
		runtime_hrs=1,
		runtime_min=59
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		bcftools norm -m -any {input.callset} | python3 workflow/scripts/prepare-for-vcfeval.py {input.ref_index} | python3 workflow/scripts/annotate.py {input.panel} | bgzip -c > {output}
		tabix -p vcf {output}
		"""


# format reference (needed for running vcfeval)
rule rtg_format_callsets:
	input:
		reference=config["reference"]
	output:
		directory("results/data/SDF")
	resources:
		mem_total_mb=20000
	conda:
		"../envs/genotyping.yml"
	shell:
		"rtg format -o {output} {input}"


# for each panel sample determine its false negatives comparing to the truth set.
# these are variants only in the truth, but not detected in the sample itself.
# later, intersect all false negatives across the panel samples, to find out
# which variants are not present in the pangenome graph, but are in the truth.
# these variants are not accessible by re-genotyping methods, that are unable
# to detect variants themselves ("untypables").
rule determine_false_negatives_vcfeval:
	input:
		truthset="results/data/vcf/{truthset}/truthset-{truthset}.vcf.gz",
		panel=config["panel_vcf"],
		reference=config["reference"],
		ref_index=config["reference_fai"],
		sdf="results/data/SDF"
	output:
		sample_vcf="results/data/vcf/{truthset}/samples/{sample}-vcfeval.vcf.gz",
		fn="results/data/vcf/{truthset}/samples/{sample}/vcfeval/fn.vcf.gz"
	wildcard_constraints:
		truthset="|".join(truthsets_small)
	conda:
		"../envs/genotyping.yml"
	params:
		tmp="results/data/vcf/{truthset}/samples/{sample}/vcfeval_temp",
		outname="results/data/vcf/{truthset}/samples/{sample}/vcfeval"
	resources:
		mem_total_mb=30000,
		runtime_hrs=0,
		runtime_min=59
	log:
		"results/data/vcf/{truthset}/samples/{sample}/vcfeval.log"
	shell:
		"""
		bcftools view --samples {wildcards.sample} {input.panel} | bcftools view --min-ac 1 | python3 workflow/scripts/prepare-for-vcfeval.py {input.ref_index} | bgzip -c > {output.sample_vcf}
		tabix -p vcf {output.sample_vcf}
		rtg vcfeval -b {input.truthset} -c {output.sample_vcf} -t {input.sdf} -o {params.tmp} --squash-ploidy  &> {log}
		mv {params.tmp}/* {params.outname}/
		rm -r {params.tmp}
		"""


# same as previous rule, but for SVs using truvari instead of vcfeval
rule determine_false_negatives_truvari:
	input:
		truthset="results/data/vcf/{truthset}/truthset-{truthset}.vcf.gz",
		panel=config["panel_vcf"],
		reference=config["reference"],
		ref_index=config["reference_fai"],
		sdf="results/data/SDF"
	output:
		sample_vcf="results/data/vcf/{truthset}/samples/{sample}-truvari.vcf.gz",
		fn="results/data/vcf/{truthset}/samples/{sample}/truvari/fn.vcf.gz"
	wildcard_constraints:
		truthset="|".join(truthsets_sv)
	conda:
		"../envs/genotyping.yml"
	params:
		tmp="results/data/vcf/{truthset}/samples/{sample}/truvari_temp",
		outname="results/data/vcf/{truthset}/samples/{sample}/truvari"
	resources:
		mem_total_mb=30000,
		runtime_hrs=0,
		runtime_min=59
	log:
		"results/data/vcf/{truthset}/samples/{sample}/truvari.log"
	shell:
		"""
		bcftools view --samples {wildcards.sample} {input.panel} | bcftools view --min-ac 1 | python3 workflow/scripts/prepare-for-vcfeval.py {input.ref_index} | bgzip -c > {output.sample_vcf}
		tabix -p vcf {output.sample_vcf}
		truvari bench -b {input.truthset} -c {output.sample_vcf} -f {input.reference} -o {params.tmp} --multimatch -r 2000 --no-ref a -C 2000 --passonly &> {log}
		bgzip {params.tmp}/fn.vcf
		tabix -p vcf {params.tmp}/fn.vcf.gz
		mv {params.tmp}/* {params.outname}/
		rm -r {params.tmp}
		"""




# intersect the sets of FNs computed for each panel sample. The intersection then defines the set of unique/untypable variants
rule determine_unique:
	input:
		truthset="results/data/vcf/{truthset}/truthset-{truthset}.vcf.gz",
		samples = expand("results/data/vcf/{{truthset}}/samples/{sample}/{{method}}/fn.vcf.gz", sample=[s for s in assembly_samples if not s == "CHM13"])
	output:
		unique_tsv="results/data/vcf/{truthset}/{truthset}-unique_{method}.tsv",
		unique_vcf="results/data/vcf/{truthset}/{truthset}-unique_{method}.vcf"
	conda:
		"../envs/genotyping.yml"
	params:
		n_files = len([s for s in assembly_samples if not s == "CHM13"])
	resources:
		mem_total_mb=30000,
		runtime_hrs=0,
		runtime_min=30
	shell:
		"""
		bcftools isec -n={params.n_files} -w1 {input.samples}  > {output.unique_vcf}
		grep -v '#' {output.unique_vcf} | cut -f 3  > {output.unique_tsv}
		"""


rule count_unique_type:
	input:
		total="results/data/vcf/{truthset}/truthset-{truthset}.vcf.gz",
		unique="results/data/vcf/{truthset}/{truthset}-unique_{method}.vcf"
	output:
		total="results/data/vcf/{truthset}/{truthset}-total_{method}-{vartype}.vcf",
		unique="results/data/vcf/{truthset}/{truthset}-unique_{method}-{vartype}.vcf"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		vartype="snp-indel|sv"
	shell:
		"""
		bcftools view {input.total} | python3 workflow/scripts/extract-varianttype.py {wildcards.vartype} > {output.total}
		bcftools view {input.unique} | python3 workflow/scripts/extract-varianttype.py {wildcards.vartype} > {output.unique}
		"""




####################################################################################################
#  compare genotyping results to the truthset, excluding untypable variants which cannot be 
# genotyped correctly by a re-genotyper
####################################################################################################


# remove untypables from VCF (callset or groundtruth) and extract only variants of a specific type.
rule extract_variant_type_callset:
	input:
		vcf= lambda wildcards: "results/data/vcf/{truthset}/truthset-{truthset}.vcf.gz".format(truthset=wildcards.truthset) if wildcards.set in truthsets else config["truthsets"][wildcards.truthset]["callsets"][wildcards.set],
		untypable="results/data/vcf/{truthset}/{truthset}-unique_{method}.tsv",
		panel = config["panel_vcf"]
	output:
		vcf="results/callset-comparisons/{truthset}/{set}-typable-{vartype}_{method}.vcf.gz",
		tbi="results/callset-comparisons/{truthset}/{set}-typable-{vartype}_{method}.vcf.gz.tbi"
	wildcard_constraints:
		truthset = "|".join(truthsets),
		set= "|".join(truthsets + callsets),
		vartype="snp-indel|sv"
	resources:
		mem_total_mb=30000
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		bcftools norm -m -any {input.vcf} | python3 workflow/scripts/annotate.py {input.panel} | python3 workflow/scripts/skip-untypable.py {input.untypable} | python3 workflow/scripts/extract-varianttype.py {wildcards.vartype} | bgzip -c > {output.vcf}
		tabix -p vcf {output.vcf}
		"""

# keep untypables and only variants of a specific type
rule extract_variant_type_callset_all:
	input:
		vcf= lambda wildcards: "results/data/vcf/{truthset}/truthset-{truthset}.vcf.gz".format(truthset=wildcards.truthset) if wildcards.set in truthsets else config["truthsets"][wildcards.truthset]["callsets"][wildcards.set],
		panel = config["panel_vcf"]
	output:
		vcf="results/callset-comparisons/{truthset}/{set}-all-{vartype}_{method}.vcf.gz",
		tbi="results/callset-comparisons/{truthset}/{set}-all-{vartype}_{method}.vcf.gz.tbi"
	wildcard_constraints:
		truthset = "|".join(truthsets),
		set= "|".join(truthsets + callsets),
		vartype="snp-indel|sv",
	resources:
		mem_total_mb=30000
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		bcftools norm -m -any {input.vcf} | python3 workflow/scripts/annotate.py {input.panel} | python3 workflow/scripts/extract-varianttype.py {wildcards.vartype} | bgzip -c > {output.vcf}
		tabix -p vcf {output.vcf}
		"""


# prepare regions (stratification BED files intersected with callable regions as defined by the BED files that come with each ground truth)
rule prepare_evaluation_beds:
	input:
		callable_regions= lambda wildcards: config["truthsets"][wildcards.truthset]["callable_regions"],
		bed= lambda wildcards: config["regions_to_bed"][wildcards.region] if wildcards.region != 'all' else config["truthsets"][wildcards.truthset]["callable_regions"]
	output:
		"results/callset-comparisons/{truthset}/bed-files/{truthset}_{region}.bed"
	wildcard_constraints:
		truthset = "|".join(truthsets),
		region = "|".join(regions)
	resources:
		mem_total_mb=30000,
		runtime_hrs=1
	conda:
		"../envs/genotyping.yml"
	shell:
		"bedtools intersect -a {input.callable_regions} -b {input.bed} > {output}"



# compute precision/recall for small variants
rule vcfeval_callsets:
	input:
		callset="results/callset-comparisons/{truthset}/{callset}-{filter}-{vartype}_vcfeval.vcf.gz",
		callset_tbi="results/callset-comparisons/{truthset}/{callset}-{filter}-{vartype}_vcfeval.vcf.gz.tbi",
		baseline="results/callset-comparisons/{truthset}/{truthset}-{filter}-{vartype}_vcfeval.vcf.gz",
		baseline_tbi="results/callset-comparisons/{truthset}/{truthset}-{filter}-{vartype}_vcfeval.vcf.gz.tbi",
		regions= "results/callset-comparisons/{truthset}/bed-files/{truthset}_{region}.bed",
		reference=config["reference"],
		ref_index = config["reference_fai"],
		sdf="results/data/SDF"
	output:
		fixed_vcf=temp("results/callset-comparisons/{truthset}/vcfeval_{callset}_{vartype}_{filter}_region-{region}_fixed.vcf.gz"),
		fixed_vcf_tbi=temp("results/callset-comparisons/{truthset}/vcfeval_{callset}_{vartype}_{filter}_region-{region}_fixed.vcf.gz.tbi"),
		summary="results/callset-comparisons/{truthset}/vcfeval_{callset}_{vartype}_{filter}_region-{region}/summary.txt"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		vartype = "snp-indel",
		filter = "all|typable"
	log:
		"results/callset-comparisons/{truthset}/vcfeval_{callset}_{vartype}_{filter}_region-{region}.log"
	params:
		tmp = "results/callset-comparisons/{truthset}/vcfeval_{callset}_{vartype}_{filter}_region-{region}_tmp",
		outname = "results/callset-comparisons/{truthset}/vcfeval_{callset}_{vartype}_{filter}_region-{region}"
	resources:
		mem_total_mb=20000,
		runtime_hrs=0,
		runtime_min=40
	shell:
		"""
		bcftools view {input.callset} --min-ac 1 | python3 workflow/scripts/prepare-for-vcfeval.py {input.ref_index} | bgzip -c > {output.fixed_vcf}
		tabix -p vcf {output.fixed_vcf}

		rtg vcfeval -b {input.baseline} -c {output.fixed_vcf} -t {input.sdf} -o {params.tmp} --evaluation-regions {input.regions} > {output.summary}.tmp
		mv {params.tmp}/* {params.outname}/
		mv {output.summary}.tmp {output.summary}
		rm -r {params.tmp}
		"""


# compute precision/recall for SVs
rule truvari_callsets:
	input:
		callset="results/callset-comparisons/{truthset}/{callset}-{filter}-{vartype}_truvari.vcf.gz",
		callset_tbi="results/callset-comparisons/{truthset}/{callset}-{filter}-{vartype}_truvari.vcf.gz.tbi",
		baseline="results/callset-comparisons/{truthset}/{truthset}-{filter}-{vartype}_truvari.vcf.gz",
		baseline_tbi="results/callset-comparisons/{truthset}/{truthset}-{filter}-{vartype}_truvari.vcf.gz.tbi",
		regions= "results/callset-comparisons/{truthset}/bed-files/{truthset}_{region}.bed",
		reference=config["reference"],
		ref_index = config["reference_fai"]
	output:
		fixed_vcf=temp("results/callset-comparisons/{truthset}/truvari_{callset}_{vartype}_{filter}_region-{region}_fixed.vcf.gz"),
		fixed_vcf_tbi=temp("results/callset-comparisons/{truthset}/truvari_{callset}_{vartype}_{filter}_region-{region}_fixed.vcf.gz.tbi"),
		summary="results/callset-comparisons/{truthset}/truvari_{callset}_{vartype}_{filter}_region-{region}/summary.txt"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		vartype = "sv",
		filter = "all|typable"
	log:
		"results/callset-comparisons/{truthset}/truvari_{callset}_{vartype}_{filter}_region-{region}.log"
	params:
		tmp = "results/callset-comparisons/{truthset}/truvari_{callset}_{vartype}_{filter}_region-{region}_tmp",
		outname = "results/callset-comparisons/{truthset}/truvari_{callset}_{vartype}_{filter}_region-{region}"
	resources:
		mem_total_mb=20000,
		runtime_hrs=0,
		runtime_min=40
	shell:
		"""
		bcftools view {input.callset} --min-ac 1 | python3 workflow/scripts/prepare-for-vcfeval.py {input.ref_index} | bgzip -c > {output.fixed_vcf}
		tabix -p vcf {output.fixed_vcf}

		truvari bench -b {input.baseline} -c {output.fixed_vcf} -f {input.reference} -o {params.tmp} --multimatch --includebed {input.regions} -r 2000 --no-ref a -C 2000 --passonly &> {log}
		mv {params.tmp}/* {params.outname}/
		rm -r {params.tmp}
		"""
