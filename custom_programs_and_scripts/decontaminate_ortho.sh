#!/bin/bash

decontaminate_ortho (){
   pseudoreference_folder=$(realpath "$1")
   folder_with_alignments=$(realpath "$2")
   output_folder=$folder_with_alignments"_decontaminated"
   [ -d $output_folder ] && echo "Output folder exists" && exit
   
   mkdir -p logs
   logs_path=$(realpath logs)
   echo "Start decontamination " "$(date)" | tee -a $logs_path/progress.txt
   
   wroking_directory=$PWD
   
   busco -i $pseudoreference_folder -o contaminants_bacteria -m transcriptome  -c 12  -l bacteria_odb10
   
   busco -i $pseudoreference_folder -o contaminants_alveolata -m transcriptome  -c 12  -l alveolata_odb10
   
   busco -i $pseudoreference_folder -o contaminants_euglenozoa -m transcriptome  -c 12  -l euglenozoa_odb10
   
   busco -i $pseudoreference_folder -o contaminants_fungi -m transcriptome  -c 12  -l fungi_odb10
   
   busco -i $pseudoreference_folder -o contaminants_nematoda -m transcriptome  -c 12  -l nematoda_odb10
   
   
   folders_busco_contaminants=$(ls -d contaminants_*)
   for i in $folders_busco_contaminants
      do
        grep -s -h ">" $i/pseudoreference.fasta/run_*_odb*/busco_sequences/fragmented_busco_sequences/*.faa | cut -f1 -d":" | cut -f2 -d">" >> possible_contaminats_list.txt
        grep -s -h ">" $i/pseudoreference.fasta/run_*_odb*/busco_sequences/single_copy_busco_sequences/*.faa | cut -f1 -d":" | cut -f2 -d">" >> possible_contaminats_list.txt
        grep -s -h ">" $i/pseudoreference.fasta/run_*_odb*/busco_sequences/multi_copy_busco_sequences/*.faa | cut -f1 -d":" | cut -f2 -d">" >> possible_contaminats_list.txt
    done
   
   contaminat_list=$(cat possible_contaminats_list.txt)
   
   echo "Preparing Output folder " "$(date)" | tee -a $logs_path/progress.txt
   cp -r $folder_with_alignments $output_folder
   
   echo "Deleting contaminat alignments " "$(date)" | tee -a $logs_path/progress.txt
   for file in $contaminat_list
     do
       rm $output_folder/$file
     done
   
 }

decontaminate_ortho "$1" "2"

#"$1" is the folder with pseudoreference
#"$2" is the folder with alignments