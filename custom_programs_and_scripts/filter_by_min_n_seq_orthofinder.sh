#!/bin/bash

filter_by_min_n_seq_orthofinder () {
   input_folder=$(realpath "$1")
   output_folder=$input_folder"_minseq""$2"
   [ -d $output_folder ] && echo "Output folder exists" && exit
   
   mkdir -p logs
   echo "Filtering for minimum number of sequences: ""$2" | tee -a logs/progress.txt
   mkdir $output_folder
   cd $input_folder
   file_list=$(find -type f -name '*.*')
   min_seq="$2"
   for i in $file_list
     do
       n_seq=$(grep -c ">" $i)
       if [[ $n_seq -ge $min_seq ]]; then
         cp $i $output_folder
       fi
     done
}

filter_by_min_n_seq_orthofinder "$1" "$2"

# $1 is the input folder with orthofinder output. If more than one result folder, only one will be processed.
# $2 is the minimun number of sequences for the alignment to be kept.