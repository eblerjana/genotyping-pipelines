#!/bin/bash
set -e

if [ $# -ne 4 ] ; then
	echo "Usage: $0 <template.vcf> <input.vcf> <output_bi_all.vcf.gz> <output_bi_nosnvs.vcf.gz>"
	exit 1
fi

template_vcf="$1"
input="$2"
output_bi_all="$3"
output_bi_nosnvs="$4"

if [ ! -f $template_vcf ] ; then
	echo "File $template_vcf not found"
	exit 1
fi

if [ ! -f $input ] ; then
	echo "File $input not found"
	exit 1
fi

which bcftools

zcat ${input} | python3 workflow/scripts/convert-to-biallelic.py ${template_vcf} | bcftools view -o ${output_bi_all} -O z
bcftools view --exclude-types snps -o ${output_bi_nosnvs} -O z ${output_bi_all}
tabix ${output_bi_nosnvs}
tabix ${output_bi_all}
