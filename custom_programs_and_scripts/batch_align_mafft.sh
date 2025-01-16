#!/bin/bash
version=0.1
# Version Log

print_message () {
  echo '
###################################
Batch Align Orthologs Version '$version'

Created by Guilherme Azevedo 2021 
###################################

This program runs MAFFT over all fasta files in a folder.
It requires MAFFT to be installed in your system. Please see https://mafft.cbrc.jp/alignment/software/ and cite MAFFT:

Katoh, Standley 2013 (Molecular Biology and Evolution 30:772-780)
MAFFT multiple sequence alignment software version 7: improvements in performance and usability.

batch_align_mafft.sh -I folder_with_orthologs [-m mafft_mode -o output_folder -n n_threads] 

      -I folder_with_orthologs      The folder containing the fasta files 
                                       to be aligned
      -m mafft_mode                 Currently only three options: E-INS-i
                                       or L-INS-i or FFT-NS-1 (see MAFFT manual).
                                       Default: E-INS-i
      -o output_folder                The name of the output file
      -n n_threads                  Number of threads to be used
                                       by job.(Default:12)
      -j n_jobs                  Number of jobs to be parallelized.(Default:1)
                                       Note that j*n should not exceed total cores.	  
      -h                            Print this message and exit.
'
exit
}



#Set variables and default values
folder_with_orthologs=""
mafft_mode="E-INS-i"
output_folder=""
n_threads=6
n_jobs=1


while getopts "I:m:o:n:j:h" flag; do
    case "${flag}" in
            I)  folder_with_orthologs=${OPTARG};;
            m)  mafft_mode=${OPTARG};;
            o)  output_folder=${OPTARG};;
            n)  n_threads=${OPTARG};;
            j)  n_jobs=${OPTARG};;
            h)  print_message;;
            ?)  print_message;;
    esac
done


if [ "$folder_with_orthologs" == "" ]
  then
     echo "The input folder with ortholog fasta sequences needs to be specified with the flag -I"
     echo
     print_message
     exit
  else
     folder_with_orthologs=$(realpath $folder_with_orthologs)
     echo "Input folder "$folder_with_orthologs""
fi

if [ "$mafft_mode" == "E-INS-i" ] || [ "$mafft_mode" == "L-INS-i" ] || [ "$mafft_mode" == "FFT-NS-1" ]
  then
     echo "MAFFT mode "$mafft_mode""
  else
     echo -e ' '$mafft_mode' Not a valid MAFFT mode option. Choose E-INS-i or L-INS-i.\n'
     print_message
     exit
fi


if [ "$mafft_mode" == "E-INS-i" ]
  then
     mafft_options="--genafpair --maxiterate 1000"
  elif [ "$mafft_mode" == "L-INS-i" ]
    then
     mafft_options="--localpair --maxiterate 1000"
  elif [ "$mafft_mode" == "FFT-NS-1" ]
    then
     mafft_options="--retree 1"    
fi

if [ "$output_folder" == "" ]
  then
     output_folder=$folder_with_orthologs"_mafft"
     [ -d $output_folder ] && echo "A folder with alignments already exists in the output directory" && exit
     mkdir $output_folder
     output_folder=$(realpath $output_folder)
     echo "Output folder "$output_folder""
  else
     [ -d $output_folder ] && echo "A folder with alignments already exists in the output directory" && exit
     mkdir $output_folder
     output_folder=$(realpath $output_folder)
     echo "Output folder "$output_folder""
fi

mkdir -p logs
logs_path=$(realpath logs)

echo '
###################################
Batch Align Orthologs Version '$version'

Created by Guilherme Azevedo 2021 
###################################

This program runs MAFFT over all fasta files in a folder.
Please cite MAFFT:

Katoh, Standley 2013 (Molecular Biology and Evolution 30:772-780)
MAFFT multiple sequence alignment software version 7: improvements in performance and usability.


----------------------------------------------------------------------
'$(date)'

Start alignning with MAFFT with options:
--thread '$n_threads' --adjustdirectionaccurately '$mafft_options'

Input folder '$folder_with_orthologs'
Output folder '$output_folder'

' >> $logs_path/progress.txt



files=$(ls $folder_with_orthologs)

for file in $files
  do
    echo "Aligning" $file >> $logs_path/progress.txt
    sem --id alignjob --jobs $n_jobs mafft --thread $n_threads --adjustdirectionaccurately $mafft_options $folder_with_orthologs/$file > $output_folder/$file
done
sem --wait --id alignjob

# remove the "_R_" from reversed sequences in the alignment
sed -i 's/>_R_/>/g' $output_folder/*.*

echo 'Finished aligning with MAFFT.

'$(date)'
----------------------------------------------------------------------
######################################################################

' >> $logs_path/progress.txt



