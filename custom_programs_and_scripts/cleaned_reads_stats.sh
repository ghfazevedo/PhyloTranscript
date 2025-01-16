#!/bin/bash

cleaned_reads_stats () {
   [ -f cleaned_reads_stats_out.csv ] && echo "Output file exists" && exit
   mkdir -p logs
   echo "Geting Cleaned Reads statistics" >> logs/progress.txt
   touch cleaned_reads_stats_out.csv
   csv_out_path=$(realpath "cleaned_reads_stats_out.csv")
   echo sample,reads,total_bp,mean_length,95_CI_length,min,max,median >> $csv_out_path
   #cleaned_reads_folder=$(realpath "$1")
   cd "$1"
   for i in *
      do
        echo "Processing $i" >> ../logs/progress.txt
        phyluce_assembly_get_fastq_lengths --input $i/split-adapter-quality-trimmed --csv >> $csv_out_path
      done
}

cleaned_reads_stats "$1"



