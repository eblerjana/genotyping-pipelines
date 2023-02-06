#!/bin/sh

input_repeats=$1
input_complex=$2
input_fai=$3

output_tmp1=$4
output_tmp2=$5
output_tmp3=$6

output_rep_complex=$7
output_rep_simple=$8
output_nonrep_complex=$9
output_nonrep_simple=${10}

# non-repetitive regions
awk 'length($1) < 6' $input_repeats | bedtools sort -i - > $output_tmp1
sort -k1,1d -k 2,2n -k 3,3n $input_fai > $output_tmp2
bedtools complement -i $output_tmp1 -g $output_tmp2 -L > $output_tmp3

# repeats + complex
bedtools intersect -a $output_tmp1 -b $input_complex > $output_rep_complex

# repeats + simple
bedtools subtract -a $output_tmp1 -b $input_complex > $output_rep_simple

# non-repetitive + complex
bedtools intersect -a $output_tmp3 -b $input_complex > $output_nonrep_complex

# non-repetitive + simple
bedtools subtract -a $output_tmp3 -b $input_complex > $output_nonrep_simple

echo $output_tmp3
echo $output_nonrep_simple
