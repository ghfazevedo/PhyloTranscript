#!/bin/bash

busco_plot () {
   busco_output="$1"
   busco_output=$(realpath "$1")
   summary_folder_prefix=$(echo $busco_output | rev | cut -d'/' -f 1 | rev)
   summary_folder=$summary_folder_prefix"_summary"
   [ -d $summary_folder ] && echo "Output folder exists" && exit
   mkdir $summary_folder
   summary_folder=$(realpath $summary_folder)
   folder_list=$(ls -d $busco_output/*.fasta)
   for species_folder in $folder_list
      do
	    cp $species_folder/*.txt $summary_folder
	done
   generate_plot.py -wd $summary_folder
 }
 
busco_plot "$1"