#!/bin/bash
version=0.2
# Version Log


print_message () {
  echo '
###################################
Batch Read Cleaner Version '$version'

Created by Guilherme Azevedo 2021 
###################################

This program runs Trimmomatic in batch mode and prepares the output for phyluce.

Usage:

batch_read_cleaner.sh -I InputFileList -R RawReads [-o] OutputFolder [-h] [... OTHER_OPTIONAL_ARGUMENTS] 

      -I  InputFileList       The comma delimited file with the name of the files to be processed folowed by the desired output name for the final assembly file
      -R  RawReads            The folder with the raw reads to be processed
      -o  OutputFolder        The name of the folder in whcih the final folder with the cleaned reads will be saved 
	                            (Optional. Deafault: clean_fastq in the working directory.)
      -a  adaptersfasta       Fastafile with the adapeter sequences.
      -m  minimum_length      Minimum length of the sequence to be kept. (Optional. Default 40 (same as in phyluce))
      -n  nthreads            Number of cores. (Optional. Default: 8)
      -h  help

See Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic) web page for more information.
Please cite:
Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.

Input file example:

sample1_R1.fastq.gz,sample1_R2.fastq.gz,Species_name_1
sample2_R1.fastq.gz,sample2_R2.fastq.gz,Species_name_2
sample3_R1.fastq.gz,sample3_R2.fastq.gz,Species_name_3
.
.
.

'
exit
}


#Set variables and default values
InputFileList=""
RawReads=""
OutputFolder=""
adaptersfasta=""
min_len=40
nthreads=8




while getopts "I:R:O:a:n:m:h" flag; do
    case "${flag}" in
            I)  InputFileList=${OPTARG};;
            R)  RawReads=${OPTARG};;
            o)  OutputFolder=${OPTARG};;
            a)  adaptersfasta=${OPTARG};;
            n)  nthreads=${OPTARG};;
            m)  min_len=${OPTARG};;
            h)  print_message;;
            ?)  printf '\nUsage: %s -I INPUT_FOLDER -O OutputFolder -R RAW_READS -a PATH_TO_ADAPTERS_FASTA [-h] [... OTHER_OPTIONAL_ARGUMENTS] \n \nRun this program with flag -h for full usage options\n' $0; exit 2;;                                                                                                                                                                                          
    esac
done


if [ "$InputFileList" == "" ]
  then
     echo "Input file needs to be specified with the flag -I"
     echo
     echo -e '\nUsage: trimmobatch -I INPUT_FOLDER  -R RAW_READS -o OutputFolder -a PATH_TO_ADAPTERS_FASTA [-h] [... OTHER_OPTIONAL_ARGUMENTS] \n \nRun this program with flag -h for full usage options\n'
     exit
  else
     InputFileList=$(realpath $InputFileList)
     echo "Input file "$InputFileList""
fi


if [ "$RawReads" == "" ]
  then
     echo "Path to folder with raw reads needs to be specified with the flag -R"
     echo
     echo -e '\nUsage: trimmobatch -I INPUT_FOLDER -O OutputFolder -R RAW_READS -a PATH_TO_ADAPTERS_FASTA [-h] [... OTHER_OPTIONAL_ARGUMENTS] \n \nRun this program with flag -h for full usage options\n'
     exit
  else
     RawReads=$(realpath $RawReads)
	 echo "Raw reads folder "$RawReads""
fi

if [ "$OutputFolder" == "" ]
  then
     OutputFolder=$(pwd)
	 echo "Output folder "$OutputFolder""
  else
     mkdir $OutputFolder
     OutputFolder=$(realpath $OutputFolder)
     echo "Output folder name "$OutputFolder""
fi

if [ "$adaptersfasta" == "" ]
  then
     echo "Path to fasta file with adapters sequence needs to be specified with the flag -a"
     echo "You can downlowad it from https://github.com/usadellab/Trimmomatic"
	 echo "Or you can download the Illumina TruSeq3-PE-2.fa the files using wget https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq3-PE-2.fa"
	 echo "You can also use a combined fats file with all paired adapters from 
	 "
     echo -e '\nUsage: trimmobatch -I INPUT_FOLDER -O OutputFolder -R RAW_READS -a PATH_TO_ADAPTERS_FASTA [-h] [... OTHER_OPTIONAL_ARGUMENTS] \n \nRun this program with flag -h for full usage options\n'
     exit
  else
     echo "Adapters file "$adaptersfasta""
fi

adaptersfasta=$(realpath $adaptersfasta)


mkdir -p $OutputFolder/logs

[ -d $OutputFolder/clean_fastq ] && echo "A folder with clean reads in the output directory already exists" && exit
mkdir $OutputFolder/clean_fastq
clean_fastq_dir=$(realpath $OutputFolder/clean_fastq)


echo '
###################################
Batch Read Cleaner Version '$version'

Created by Guilherme Azevedo 2021 
###################################

This program runs Trimmomatic in batch mode

Please cite Trimmomatic:
Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.

----------------------------------------------------------------------
'$(date)'

Output folder: '$OutputFolder'

Raw reads folder: '$RawReads'

Start cleaning raw reads with Trimmomatic '$(trimmomatic -version)'

Commands used:
PE -threads '$nthreads' ILLUMINACLIP:'$adaptersfasta':2:30:10:2:keepBothReads LEADING:5 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:'$min_len'

' >> $OutputFolder/logs/progress.txt

progressfile=$(realpath $OutputFolder/logs/progress.txt)


filelist=$(paste -d " " $InputFileList)


for line in $filelist
  do
    r1=$(echo $line | cut -f1 -d",")
    r2=$(echo $line | cut -f2 -d",")
    out=$(echo $line | cut -f3 -d",")
    
    echo "Cleaning "$out"" >> $progressfile
	trimmomatic PE -threads $nthreads $RawReads/$r1 $RawReads/$r2 $out"-READ1.fastq.gz" $out"-READ1-single.fastq.gz" $out"-READ2.fastq.gz" $out"-READ2-single.fastq.gz" ILLUMINACLIP:$adaptersfasta:2:30:10:2:keepBothReads LEADING:5 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:$min_len
	cat $out"-READ1-single.fastq.gz" $out"-READ2-single.fastq.gz" > $out"-READ-singleton.fastq.gz"
	rm $out"-READ1-single.fastq.gz" $out"-READ2-single.fastq.gz"
	mkdir $clean_fastq_dir/$out/split-adapter-quality-trimmed
	mv *READ*.fastq.gz $clean_fastq_dir/$out/split-adapter-quality-trimmed/
done


echo "Finished cleaning raw reads with Trimmomatic

'$(date)'
----------------------------------------------------------------------
######################################################################

" >> $progressfile

