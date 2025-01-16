#!/bin/bash
version=0.1
# Version Log

print_message () {
  echo '
####################################
Remove Redundancy Version '$version'

Created by Guilherme Azevedo 2021 
####################################

This program uses a loop to run CD_HIT in a folder with fasta files containing contigs of different species.

Please cite CD-HIT:
Weizhong Li & Adam (2006). Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences.  Godzik Bioinformatics,  22:1658-9. Open access PDF Pubmed Citations

Limin Fu, Beifang Niu, Zhengwei Zhu, Sitao Wu and Weizhong Li (2012). CD-HIT: accelerated for clustering the next generation sequencing data. Bioinformatics, 28 (23): 3150-3152. doi: 10.1093/bioinformatics/bts565 


Usage:

remove_redundancy.sh -I ContigsFolder -O OutputPrefix [-h] [... OTHER_OPTIONAL_ARGUMENTS] 

      -I  ContigsFolder		  The folder with assembled the contigs
      -O  OutputPrefix        The prefix name for the folder to save the otputs (Optional. Default: input folder)
      -c  Similarity          Similarity threshold for CD-HIT. (Otional. Default 0.95)
      -t  nthreads            Number of cores. Default: 8
      -h  help

'
exit
}


#Set variables and default values
ContigsFolder=""
OutputPrefix=""
Similarity=0.95
nthreads=8


while getopts "I:O:c:t:h" flag; do
    case "${flag}" in
            I)  ContigsFolder=${OPTARG};;
            O)  OutputPrefix=${OPTARG};;
			c)  Similarity=${OPTARG};;
			t)  nthreads=${OPTARG};;
            h)  print_message;exit;;
            ?)  printf '\nremove_redundancy.sh -I ContigsFolder -O OutputPrefix [-h] [... OTHER_OPTIONAL_ARGUMENTS]  \n \nRun this program with flag -h for full usage options\n' $0; exit 2;;                                                                                                                                                                                          
    esac
done


if [ "$ContigsFolder" == "" ]
  then
     echo "Input file needs to be specified with the flag -I"
     echo
     echo -e '\nUsage: remove_redundancy.sh -I ContigsFolder -O OutputPrefix [-h] [... OTHER_OPTIONAL_ARGUMENTS] \n'
     exit
  else
     ContigsFolder=$(realpath $ContigsFolder)
     echo "Input file "$ContigsFolder""
fi

if [ "$OutputPrefix" == "" ]
  then
     OutputFolderClstr=$ContigsFolder"_cdhitclstr"
	 OutputFolderFasta=$ContigsFolder"_cdhitfasta"
	 [ -d $OutputFolderClstr ] && echo "A folder with CD_HIT output already exists in the output directory" && exit
	 [ -d $OutputFolderFasta ] && echo "A folder with CD_HIT output already exists in the output directory" && exit
	 mkdir $OutputFolderClstr $OutputFolderFasta
	 OutputFolderClstr=$(realpath $OutputFolderClstr)
	 OutputFolderFasta=$(realpath $OutputFolderFasta)
	 echo "Output folder "$OutputFolderClstr" and "$OutputFolderFasta""
  else
     OutputFolderClstr=$OutputPrefix"_cdhitclstr"
	 OutputFolderFasta=$OutputPrefix"_cdhitfasta"
	 [ -d $OutputFolderClstr ] && echo "A folder with CD_HIT output already exists in the output directory" && exit
	 [ -d $OutputFolderFasta ] && echo "A folder with CD_HIT output already exists in the output directory" && exit
	 mkdir $OutputFolderClstr $OutputFolderFasta
	 OutputFolderClstr=$(realpath $OutputFolderClstr)
	 OutputFolderFasta=$(realpath $OutputFolderFasta)
	 echo "Output folder "$OutputFolderClstr" and "$OutputFolderFasta""
fi



mkdir -p $PWD/logs

echo '
####################################
Remove Redundancy Version '$version'

Created by Guilherme Azevedo 2021 
####################################

This program uses a loop to run CD_HIT in a folder with fasta files containing contigs of different species.


Please cite CD-HIT:
Weizhong Li & Adam (2006). Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences.  Godzik Bioinformatics,  22:1658-9. Open access PDF Pubmed Citations

Limin Fu, Beifang Niu, Zhengwei Zhu, Sitao Wu and Weizhong Li (2012). CD-HIT: accelerated for clustering the next generation sequencing data. Bioinformatics, 28 (23): 3150-3152. doi: 10.1093/bioinformatics/bts565 


----------------------------------------------------------------------
'$(date)'

Input folder: '$ContigsFolder'
Output folder: '$OutputFolderClstr' and '$OutputFolderFasta'

Start removing redundant sequences with commands:

cd-hit-est -c $Similarity -n 10 -M 0 -T $nthreads -sf 1

' >> $PWD/logs/progress.txt



for file in $(ls $ContigsFolder)
  do
    echo "Processing "$file"" >> $PWD/logs/progress.txt
    cd-hit-est -i $ContigsFolder/$file -o $OutputFolderClstr/$file -c $Similarity -n 10 -M 0 -T $nthreads -sf 1
  done

mv $OutputFolderClstr/*.fasta $OutputFolderFasta/


echo 'Start removing redundant sequences


'$(date)'
----------------------------------------------------------------------
######################################################################

' >> $PWD/logs/progress.txt

