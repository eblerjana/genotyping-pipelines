#!/bin/bash
set -e

if [ $# -ne 3 ] ; then
	echo "Usage: $0 <input.vcf.gz> <output_multi_all.vcf.gz> <output_multi_nosnvs.vcf.gz>"
	exit 1
fi

input="$1"
output_multi_all="$2"
output_multi_nosnvs="$3"

if [ ! -f $input ] ; then
	echo "File $input not found"
	exit 1
fi

# source ~/.bashrc
# conda activate bcftools
which bcftools
bcftools view -o ${output_multi_all} -O z ${input}
bcftools view --exclude-types snps -o ${output_multi_nosnvs} -O z ${input}
tabix ${output_multi_nosnvs}
tabix ${output_multi_all}
