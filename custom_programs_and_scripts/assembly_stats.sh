#!/bin/bash

assembly_stats () {
   output="$1""_stats.csv"
   [ -f "$output" ] && echo "Output file exists" && exit
   mkdir -p logs
   echo "Geting Cleaned Reads statistics" >> logs/progress.txt
   touch "$output"
   csv_out_path=$(realpath "$output")
   echo samples,contigs,total_bp,mean_length,95_CI_length,min_length,max_length,median_legnth,contigs_bigger1kb >> $csv_out_path
   #cleaned_reads_folder=$(realpath "$1")
   cd "$1"
   for i in *
      do
        echo "Processing $i" >> ../logs/progress.txt
        phyluce_assembly_get_fasta_lengths --input $i --csv >> $csv_out_path
      done
}

assembly_stats "$1"



