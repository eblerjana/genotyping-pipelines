import os.path
import gzip

# set panel. If PanGenie-ready panels are given, use these. Otherwise,
# create them from provided MC graph

PANEL_BI = config["panel_bi"]
PANEL_MULTI = config["panel_multi"]
MC_GFA = config["mc_gfa"]
MC_VCF = config["mc_vcf"]
MIN_FRAC = 0.8
SAMPLES_TO_EXCLUDE = config["samples_to_exclude"]
PANEL_SAMPLES = []


if not PANEL_BI or not PANEL_MULTI:
	# make sure MC data is provided
	assert MC_GFA and MC_VCF, "No panels and no MC graph provided."
	PANEL_BI = "{results}/mc-vcf/mc_filtered_ids_biallelic.vcf.gz" # results of prepare-mc-vcf pipeline
	PANEL_MULTI = "{results}/mc-vcf/mc_filtered_ids.vcf.gz" # results of prepare-mc-vcf pipeline
	# determine panel samples
	for line in gzip.open(MC_VCF, 'rt'):
		if line.startswith('##'):
			continue
		if line.startswith('#'):
			fields = line.strip().split()
			PANEL_SAMPLES = fields[9:]
			break
else:
	for line in gzip.open(PANEL_BI, 'rt'):
		if line.startswith('##'):
			continue
		if line.startswith('#'):
			fields = line.strip().split()
			PANEL_SAMPLES = fields[9:]
			break


REFERENCE = config["reference"]
assert os.path.isfile(REFERENCE), "File " + REFERENCE + " does not exist."
assert os.path.isfile(REFERENCE + ".fai"), "File " + REFERENCE + ".fai (index) does not exist."

SAMPLE_SHEET = config["sample-sheet"]
ILLUMINA = {}
ONT = {}
READ_TECH = {}
CRAM_REF = config["cram_ref"]
WINDOW_WIDTHS = [25]

assert os.path.isfile(SAMPLE_SHEET), "File " + SAMPLE_SHEET + " does not exist."

for line in open(SAMPLE_SHEET, 'r'):
	if line.startswith("#"):
		continue
	fields = line.strip().split()
	if fields[7] == "nan":
		continue
	ILLUMINA[fields[1]] = fields[7]

	assert os.path.isfile(fields[7]), "File " + fields[7] + " does not exist."
	if fields[8] != "nan":
		ONT[fields[1]] = fields[8]
		assert fields[9] in ["ONT", "HIFI"], "Read technology must either be ONT or HIFI."
		READ_TECH[fields[1]] = fields[9]

AGC = {
	"pangenie_all-samples_unfiltered": "{results}/unpolished-haplotypes/pangenie_all-samples_unfiltered_unpolished.agc"
}

CONSENSUS = {}
CONSENSUS["pangenie_all-samples_unfiltered"] = {}

for sample in ONT:
	for haplotype in ["hap1", "hap2"]:
		CONSENSUS["pangenie_all-samples_unfiltered"][(sample, haplotype)] = "pangenie_all-samples_unfiltered_unpolished_" + sample + "_" + haplotype
SAMPLES = [s for s in ONT.keys()]


MAPS = config["maps"]
# parse chromosomes and define regions needed for merging
CHROMOSOMES = []
MERGING_REGIONS = []

step_size = 1000000
for line in open(REFERENCE + ".fai", 'r'):
	fields = line.strip().split()
	chrom = fields[0]
	CHROMOSOMES.append(chrom)
	pos = 1
	chrom_len = int(fields[1])
	while pos < chrom_len:
		end_pos = min(pos+step_size-1, chrom_len)
		MERGING_REGIONS.append('{}:{}-{}'.format(chrom, pos, end_pos))
		pos = pos + step_size
