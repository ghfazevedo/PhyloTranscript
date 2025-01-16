#!/bin/bash
version=0.1
# Version Log

print_message () {
  echo '
###################################
Create Pseudoreference '$version'

Created by Guilherme Azevedo 2022 
###################################

This program runs CIAlign over all aligned fasta files in a folder to create a pseudoreference fasta file.
It requires CIAlign to be installed in your system. Please see https://github.com/KatyBrown/CIAlign and cite CIAlign:

Tumescheit C, Firth AE, Brown K. 2022. CIAlign: A highly customisable command line tool to clean, interpret and visualise multiple sequence alignments. PeerJ 10:e12983 https://doi.org/10.7717/peerj.12983

Usage:
create_pseudoreference.sh -I input_folder [-d min_diver -o output_folder -i min_insertion] 

      -A input_folder         The folder containing the fasta files 
                                with aligned sequences.
      -o output_folder        The name of the output folder. 
                                (Default: input_folder_CIAlign)
      -j n_jobs               Number of jobs to be parallelized.(Default:1)
                                 Note that j*n should not exceed total number of cores.
      -h                            Print this message and exit.
'
exit
}



#Set variables and default values
input_folder=""
output_folder=""
n_jobs=1


while getopts "A:o:j:h" flag; do
    case "${flag}" in
            A)  input_folder=${OPTARG};;
            o)  output_folder=${OPTARG};;
            h)  print_message;;
            j)  n_jobs=${OPTARG};;
            ?)  print_message;;
    esac
done

if [ "$input_folder" == "" ]
  then
     echo "The input folder with aligned fasta sequences needs to be specified with the flag -A"
     echo
     print_message
     exit
  else
     input_folder=$(realpath $input_folder)
	 echo "Input folder "$input_folder""
fi


if [ "$output_folder" == "" ]
  then
     output_folder=$input_folder"_pseudoreference"
	 [ -d $output_folder ] && echo "A folder with pseudoreference already exists in the output directory" && exit
	 mkdir $output_folder
	 output_folder=$(realpath $output_folder)
	 echo "Output folder "$output_folder""
  else
	 [ -d $output_folder ] && echo "A folder with pseudoreference already exists in the output directory" && exit
	 mkdir $output_folder
	 output_folder=$(realpath $output_folder)
	 echo "Output folder "$output_folder""
fi

mkdir -p logs
logs_path=$(realpath logs)

echo '
###################################
Create Pseudoreference '$version'

Created by Guilherme Azevedo 2021 
###################################

This program runs CIAlign over all aligned fasta files in a folder to create a pseudoreference fasta file.
cite CIAlign:

Tumescheit C, Firth AE, Brown K. 2022. CIAlign: A highly customisable command line tool to clean, interpret and visualise multiple sequence alignments. PeerJ 10:e12983 https://doi.org/10.7717/peerj.12983

----------------------------------------------------------------------
'$(date)'

Start CIAlign:
CIAlign --make_consensus --consensus_name $file 

Input folder '$input_folder'
Output folder '$output_folder'

' | tee -a  $logs_path/progress.txt

file_list=$(ls $input_folder)

for file in $file_list
    do
      sem --jobs $n_jobs CIAlign --infile $input_folder/$file --outfile_stem $output_folder/$file --make_consensus --consensus_name $file
    done	  
sem --wait	  

rm $output_folder/*.fasta_with_consensus.fasta
cat $output_folder/*_consensus.fasta >> $output_folder/pseudoreference.fasta
rm $output_folder/*_consensus.fasta

cat $output_folder/*.fasta_log.txt >> $logs_path/CIAlign_logs_consensus.txt

rm $output_folder/*.txt
rm $output_folder/*.fasta_cleaned.fasta

echo 'Finished creating pseudoreference.

'$(date)'
----------------------------------------------------------------------
######################################################################

' | tee -a $logs_path/progress.txt



