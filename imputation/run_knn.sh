#!/bin/bash
k_param=$1
input_file=$2
output_file=$3

python3 knn_smooth.py -k ${k_param} -f ${input_file} -o intermediate_${k_param}.txt
Rscript postprocessing.R intermediate_${k_param}.txt ${output_file} knn_${k_param}