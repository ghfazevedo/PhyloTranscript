#!/bin/bash
version=0.1
# Version Log

print_message () {
  echo '
###################################
Batch Quality Control Version '$version'

Created by Guilherme Azevedo 2021 
###################################

This program runs FastQC in batch mode taking as input a folder with cleaned reads. The folder needs to be organized as:
    clean_reads
	  |
	  |_ Species_name_CODE01
	  |  |
	  |  |_split-adapter-quality-trimmed
	  |	   |
	  |    |Species_name_CODE01-READ1.fastq.gz
	  |    |Species_name_CODE01-READ2.fastq.gz
	  |    |Species_name_CODE01-READ-singleton.fastq.gz
	  |
	  |_ Species_name_CODE02
	  |  |
	  |  |_split-adapter-quality-trimmed
	  |	   |
	  |    |Species_name_CODE02-READ1.fastq.gz
	  |    |Species_name_CODE02-READ2.fastq.gz
	  |    |Species_name_CODE02-READ-singleton.fastq.gz
	  |
	  ...

This is the same folder structure output by the illumiprocessor in the phyluce softwere for UCEs or output by the custom program batch_read_cleaner.sh.
  

Usage:

batch_quality_control -C clean_reads  [-o] OutputFolder [-h] [... OTHER_OPTIONAL_ARGUMENTS] 

      -I  InputFolder	Folder with the clean reads.
      -o  OutputFolder  Name of the output folder. (Optional. Default: input folder with sufix "_fastqc_output" 
      -h  help

See FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) web page for more information.
'
}


#Set variables and default values
InputFolder=""
OutputFolder=""


while getopts "I:o:h" flag; do
    case "${flag}" in
            I)  InputFolder=${OPTARG};;
            o)  OutputFolder=${OPTARG};;
            h)  print_message;exit;;
            ?)  printf '\nUsage: %s batch_quality_control -C clean_reads  [-o] OutputFolder [-h] [... OTHER_OPTIONAL_ARGUMENTS]  \n \nRun this program with flag -h for full usage options\n' $0; exit 2;;                                                                                                                                                                                          
    esac
done

if [ "$OutputFolder" == "" ]
  then
     OutputFolder=$InputFolder"_fastqc_output"
	 [ -d $OutputFolder ] && echo "Output directory exists" && exit
	 mkdir $OutputFolder
	 OutputFolder=$(realpath $OutputFolder)
	 echo "Output folder "$OutputFolder""
  else
     [ -d $OutputFolder ] && echo "Output directory exists" && exit
     mkdir $OutputFolder
	 OutputFolder=$(realpath $OutputFolder)
     echo "Output folder name "$OutputFolder""
fi

if [ "$InputFolder" == "" ]
  then
     echo "Input folder needs to be specified with the flag -I"
     echo
     print_message
     exit
  else
     InputFolder=$(realpath $InputFolder)
     echo "Input folder "$InputFolder""
fi

mkdir -p logs
log_dir=$(realpath logs)


echo "
###################################
Batch Quality Control Version '$version'

Created by Guilherme Azevedo 2021 
###################################

This program runs FastQC in batch mode taking as input a folder with cleaned reads.

----------------------------------------------------------------------
'$(date)'

Output folder: '$OutputFolder'

Start FastQC '$(fastqc -v)' with default configuration.

" >> $log_dir/progress.txt

progressfile=$(realpath $log_dir/progress.txt)


species_folders=$(ls $InputFolder)
for species in $species_folders
  do
    echo "Analyzing "$species >> $progressfile
    fastqc -o $OutputFolder $InputFolder/$species/split-adapter-quality-trimmed/*.fastq.gz
  done

rm -r $OutputFolder/*.zip


cd $OutputFolder
for i in *.*; 
  do 
    # Adapter content
    adapterW=$(grep -c 'alt\=\"\[WARN\]\"\/>Adapter\ Content' $i) ; 
      if [[ $adapterW -gt 0 ]]; 
        then 
          echo $i >> $log_dir/FastQC_CheckAdapterContent.txt; 
      fi; 
    adapterF=$(grep -c 'alt\=\"\[FAIL\]\"\/>Adapter\ Content' $i) ; 
      if [[ $adapterF -gt 0 ]]; 
        then 
          echo $i >> $log_dir/FastQC_CheckAdapterContent.txt; 
      fi; 

    #Per base sequence quality
    basequalW=$(grep -c 'alt\=\"\[WARN\]\"\/>Per\ base\ sequence\ quality' $i) ; 
      if [[ $basequalW -gt 0 ]]; 
        then 
          echo $i >> $log_dir/FastQC_CheckBaseQuality.txt; 
      fi; 
    basequalF=$(grep -c 'alt\=\"\[FAIL\]\"\/>Per\ base\ sequence\ quality' $i) ; 
      if [[ $basequalF -gt 0 ]]; 
        then 
          echo $i >> $log_dir/FastQC_CheckBaseQuality.txt; 
      fi;	  
	  
    #Overrepresented sequences
    overrepW=$(grep -c 'alt\=\"\[WARN\]\"\/>Overrepresented\ sequences' $i) ; 
      if [[ $overrepW -gt 0 ]]; 
        then 
          echo $i >> $log_dir/FastQC_CheckOverrepresented.txt; 
      fi; 
    overrepF=$(grep -c 'alt\=\"\[FAIL\]\"\/>Overrepresented\ sequences' $i) ; 
      if [[ $overrepF -gt 0 ]]; 
        then 
          echo $i >> $log_dir/FastQC_CheckOverrepresented.txt; 
      fi;	  
	  
done

echo "
Finished Quality Control

'$(date)'
----------------------------------------------------------------------
######################################################################

" >> $log_dir/progress.txt



