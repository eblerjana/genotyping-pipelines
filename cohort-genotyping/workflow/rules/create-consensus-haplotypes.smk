

rule consensus_compute_consensus:
	"""
	Insert all variants into the reference genome to produce haplotypes.
	"""
	input:
		vcf = "{results}/phasing/{callset}_shapeit.bcf",
		reference = REFERENCE
	output:
		temp("{results}/unpolished-haplotypes/{callset}/{callset}_unpolished_{sample}_hap{haplotype}.fasta.gz")
	log:
		"{results}/unpolished-haplotypes/{callset}/{callset}_unpolished_{sample}_hap{haplotype}.log"
	benchmark:
		"{results}/unpolished-haplotypes/{callset}/{callset}_unpolished_{sample}_hap{haplotype}.benchmark.txt"
	conda:
		"../envs/shapeit.yaml"
	wildcard_constraints:
		haplotype = "1|2"
	resources:
		mem_mb = 20000,
		walltime = "04:00:00"
	shell:
		"""
		bcftools consensus --sample {wildcards.sample} --haplotype {wildcards.haplotype} -e 'ALT~\"<.*>\"' -f {input.reference} {input.vcf} 2> {log} | bgzip > {output}
		"""


rule consensus_compress_haplotypes:
	input:
		genomes = expand("{{results}}/unpolished-haplotypes/{{callset}}/{{callset}}_unpolished_{sample}_hap{haplotype}.fasta.gz", sample = ILLUMINA.keys() , haplotype = ["1", "2"]),
		reference = REFERENCE
	output:
		"{results}/unpolished-haplotypes/{callset}_unpolished.agc"
	log:
		"{results}/unpolished-haplotypes/{callset}_unpolished.log"
	benchmark:
		"{results}/unpolished-haplotypes/{callset}_unpolished.benchmark.txt"
	conda:
		"../envs/agc.yaml"
	resources:
		mem_mb = 80000,
		walltime = "20:00:00"
	threads: 32
	shell:
		"""
		agc create {input.reference} {input.genomes} -o {output} -t {threads} &> {log}
		"""


