#!/bin/bash

filter_by_min_n_seq () {
   output_folder="$1""_minseq""$2"
   [ -d $output_folder ] && echo "Output folder exists" && exit
   
   mkdir -p logs
   echo "Filtering for minimum number of sequences" >> logs/progress.txt
   input_folder=$(realpath "$1")
   mkdir $output_folder
   output_folder=$(realpath $output_folder)
   
   file_list=$(ls $input_folder/*.*)
   min_seq="$2"
   for i in $file_list
     do
       n_seq=$(grep -c ">" $i)
       if [[ $n_seq -ge $min_seq ]]; then
         cp $i $output_folder
       fi
     done
}

filter_by_min_n_seq "$1" "$2"

# $1 is the input folder with alignments to filter_by_min_n_seq
# $2 is the minimun number of sequences for the alignment to be kept.