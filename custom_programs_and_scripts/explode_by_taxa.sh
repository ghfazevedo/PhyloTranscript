#!/bin/bash

print_message(){
echo '
##############################################################
##                Explode by taxa                            ##
###############################################################

 This program uses phyluce function phyluce_align_filter_alignments
 phyluce needs to be installed and the environment needs to be active
 The only argument is the folder with alignments

 Usage: explode_by_taxa.sh folder_with_alignments

See phyluce_align_explode_alignments --help and the phyluce webpage 
'
exit
}


explode_by_taxa(){
    input_folder=$(realpath "$1")
	output_folder=$input_folder"_bytaxa"

   [ -d "$output" ] && echo "Output folder exists:" "$output" && exit
   mkdir -p logs
   logs_path=$(realpath logs)
   echo "Exploding alignments" | tee -a $logs_path/progress.txt

   phyluce_align_explode_alignments \
           --alignments $input_folder \
           --input-format fasta \
           --output $output_folder \
           --by-taxon \
           --include-locus
   
}


if [ "$1" == "" ]
  then
     print_message
	 exit
  else
     explode_by_taxa "$1"
fi


