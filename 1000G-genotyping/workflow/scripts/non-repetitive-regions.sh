#!/bin/sh

input_repeats=$1
output_subset=$2
input_fai=$3
output_fai=$4
output_bed=$5

awk 'length($1) < 6' $input_repeats | bedtools sort -i - > $output_subset
sort -k1,1d -k 2,2n -k 3,3n $input_fai > $output_fai
bedtools complement -i $output_subset -g $output_fai -L > $output_bed
