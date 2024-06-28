#!/bin/bash

mkdir ./results/
output_dir="./results/"

for file in *.fas; do

  output_file="${output_dir}${file%.*}.csv"
  output_file_sum="${output_dir}${file%.*}_sum.csv"

  blastn -db ~/Reference/blastdb/16S_ribosomal_RNA \
         -query "${file}" \
         -task blastn \
         -dust no \
         -outfmt "7 delim=, qacc sacc evalue bitscore qcovus pident sscinames" \
         -max_target_seqs 5 >  "$output_file"
  grep -v '#'  "$output_file" > "$output_file_sum"

  output_file_t1="${output_dir}${file%.*}_t1.csv"
  output_file_t1_sum="${output_dir}${file%.*}_t1_sum.csv"

  blastn -db ~/Reference/blastdb/16S_ribosomal_RNA \
         -query "${file}" \
         -task blastn \
         -dust no \
         -outfmt "7 delim=, qacc sacc evalue bitscore qcovus pident sscinames" \
         -max_target_seqs 1 >  "$output_file_t1"
  grep -v '#'  "$output_file_t1" > "$output_file_t1_sum"

done
