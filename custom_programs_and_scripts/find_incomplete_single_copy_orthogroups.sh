#!/bin/bash
version=0.1
# Version Log

print_message () {
  echo '
##################################################
Find Incomplete Single Copy Orthogroups '$version'

Created by Guilherme Azevedo 2021 
##################################################

This program takes as input a folder with OrthoFinder result and retrieve all single copy orthologs with incomplete taxa occupancy. The folder must be organized as:
    orthofinder_results
		    |
            |_ Citation.txt
            |_ Comparative_Genomics_Statistics
            |_ Log.txt
            |_ Orthogroup_Sequences
            |_ Orthogroups
            |_ Single_Copy_Orthologue_Sequences
            |_ WorkingDirectory



Usage:

find_incomplete_single_copy_orthogroups.sh -I orthofinder_output_folder -N number_of_species [-h]

        -I orthofinder_output_folder    The folder with OrthoFinder results
        -N number_of_species            The total number of species 
                                          expected in a complete
                                          occupancy matrix.
        -h                              Print this message  



'
exit
}


#Set variables and default values
orthofinder_output_folder=""
number_of_species=""



while getopts "I:N:h" flag; do
    case "${flag}" in
            I)  orthofinder_output_folder=${OPTARG};;
            N)  number_of_species=${OPTARG};;
            h)  print_message;;
            ?)  print_message;;
    esac
done

if [ "$orthofinder_output_folder" == "" ]
  then
     echo "The input folder with orthofinder results must to be specified with the flag -I"
     echo
     print_message
     exit
  else
     orthofinder_output_folder=$(realpath $orthofinder_output_folder)
	 echo "Input folder is "$orthofinder_output_folder""
fi

if [ "$number_of_species" == "" ]
  then
     echo "The number of species expected in a comple occupancy matrix must to be specified with the flag -N"
     echo
     print_message
     exit
  else
     echo "Total number of species is "$number_of_species""
fi


output_folder=$(realpath $orthofinder_output_folder/Orthogroup_Sequences_Incomplete_Single_Copy)

[ -d $output_folder ] && echo "A folder with incomplete orthofinder orthogroups already exists in the output directory" && exit
echo "Output folder "$output_folder""

[ -f $orthofinder_output_folder/Orthogroups/not_single_copy_list.txt ] && echo "A file with not sigle copy genes exists at "$orthofinder_output_folder/Orthogroups/not_single_copy_list.txt"" && exit

mkdir -p logs
logs_path=$(realpath logs)

echo '
##################################################
Find Incomplete Single Copy Orthogroups '$version'

Created by Guilherme Azevedo 2021 
##################################################

This program takes as input a folder with OrthoFinder result and retrieve all single copy orthologs with incomplete taxa occupancy.


----------------------------------------------------------------------
'$(date)'

Start retrieving inconplete occupancy single copy orthogroups:

Input folder '$orthofinder_output_folder'
Output folder '$output_folder'
Total number of species: '$number_of_species'

' | tee -a $logs_path/progress.txt


touch $orthofinder_output_folder/Orthogroups/not_single_copy_list.txt
while IFS= read -r line
 do 
   list=$(echo $line | cut -f2- -d" " | cut -f1-$number_of_species -d" " )
   gene=$(echo $line | cut -f1 -d" ")
   echo "Screening orthogroup $gene" 
   for count in $list
     do
       if [[ $count -gt 1 ]]; then
         echo $gene >> $orthofinder_output_folder/Orthogroups/not_single_copy_list.txt
         break
       fi
     done
  done < $orthofinder_output_folder/Orthogroups/Orthogroups.GeneCount.tsv
  
# create folder with all single copy
echo "Preparing single copy orthogroup folder" | tee -a $logs_path/progress.txt

sem --jobs +0 cp -r $orthofinder_output_folder/Orthogroup_Sequences $output_folder
sem --wait

while IFS= read -r line
 do 
   sem --jobs +0 rm $output_folder/$line".fa"
 done < $orthofinder_output_folder/Orthogroups/not_single_copy_list.txt
sem --wait

working_directory=$(pwd)
cd $output_folder
sem --jobs +0 rename .fa .fasta *.fa
sem --wait

cd $working_directory


echo 'Finished retrieving inconplete occupancy single copy orthogroups.

'$(date)'
----------------------------------------------------------------------
######################################################################

' | tee -a $logs_path/progress.txt

