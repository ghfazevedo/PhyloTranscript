#!/bin/bash

get_unique_orthofinder (){
   input_folder_with_orthofinder_results=$(realpath "$1")
   output_folder=$input_folder_with_orthofinder_results/Orthogroup_Sequences_Incomplete_Single_Copy_minseqfiltered_mafft_CIAlign_decontaminated_uniq
   [ -d $output_folder ] && echo "Output folder exists" && exit
   
   mkdir -p logs
   logs_path=$(realpath logs)
   echo "Getting Unique Orthofinder " "$(date)" | tee -a $logs_path/progress.txt
   
   wroking_directory=$PWD
   
   cd $input_folder_with_orthofinder_results
   
   busco -i Orthogroup_Sequences_Incomplete_Single_Copy_minseqfiltered_mafft_CIAlign_consensus -o arachnida -m transcriptome  -c 12  -l arachnida_odb10
   

   	  grep -s -h ">" arachnida/pseudo_reference.fasta/run_arachnida_odb10/busco_sequences/fragmented_busco_sequences/*.faa | cut -f1 -d":" | cut -f2 -d">" >> $input_folder_with_orthofinder_results/BUSCO_list.txt
      grep -s -h ">" arachnida/pseudo_reference.fasta/run_arachnida_odb10/busco_sequences/single_copy_busco_sequences/*.faa | cut -f1 -d":" | cut -f2 -d">" >> $input_folder_with_orthofinder_results/BUSCO_list.txt
   	  grep -s -h ">" arachnida/pseudo_reference.fasta/run_arachnida_odb10/busco_sequences/multi_copy_busco_sequences/*.faa | cut -f1 -d":" | cut -f2 -d">" >> $input_folder_with_orthofinder_results/BUSCO_list.txt

   
   BUSCO_list=$(cat BUSCO_list.txt)
   
   echo "Preparing Output folder " "$(date)" | tee -a $logs_path/progress.txt
   cp -r Orthogroup_Sequences_Incomplete_Single_Copy_minseqfiltered_mafft_CIAlign_decontaminated $output_folder
   
   echo "Deleting BUSCO alignments " "$(date)" | tee -a $logs_path/progress.txt
   for file in $BUSCO_list
     do
       rm $output_folder/$file
     done
   
   
   cd $wroking_directory
 }

get_unique_orthofinder "$1"