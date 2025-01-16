#!/bin/bash
version=0.4
# Version Log
# v0.5: added the reverse-forward orientation option
# v0.4: added the --careful option and --cov-cutoff auto as default
# v0.3: added the single cell/isolate option
# v0.2: corrected the number of threads option
print_message () {
  echo '
###################################
Batch Assembly with SPAdes Version '$version'

Created by Guilherme Azevedo 2021 
###################################

This program runs SPAdes in a folder containing clean reads and outputs a folder with all contigs. The input folder needs to have the following directory structure:

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

Usage:

spadesbatch -I CleanReadsFolder -O OutputFolder [-h] [... OTHER_OPTIONAL_ARGUMENTS] 

      -I  CleanReadsFolder    The folder with the clean reads (trimmomatic output)
      -O  OutputFolder        The name of the folder in whcih the final assemblies will be saved
      -s  spadesoptions       The optional arguments of spades: isolate or sc (single cell) (default: sc).
      -o  orientation         The orientation of reads: rf (reverse-forward) or fr (forward-reverse). Default: fr
      -n  nthreads            Number of cores. Default: 8
      -m  memory              Memory limit in Gb (Deffault 250) 
      -j n_jobs               Number of jobs to be parallelized.(Default:1)
                               Note that j*n should not exceed total number of cores.
      -h  help

See Spades (https://github.com/ablab/spades) web page for more information.


'
exit
}


#Set variables and default values
CleanReadsFolder=""
OutputFolder=""
nthreads=8
memory=250
spadesoptions="sc"
orientation="fr"
n_jobs=1


while getopts "I:o:s:n:m:o:j:h" flag; do
    case "${flag}" in
            I)  CleanReadsFolder=${OPTARG};;
            O)  OutputFolder=${OPTARG};;
            s)  spadesoptions=${OPTARG};;
            o)  orientation=${OPTARG};;
            n)  nthreads=${OPTARG};;
            m)  memory=${OPTARG};;
            j)  n_jobs=${OPTARG};;			
            h)  print_message;;
            ?)  printf '\nUsage: %s -I INPUT_FOLDER -O OUTPUT_PREFIX [-h] [... OTHER_OPTIONAL_ARGUMENTS] \n \nRun this program with flag -h for full usage options\n' $0; exit 2;;                                                                                                                                                                                          
    esac
done


if [ "$CleanReadsFolder" == "" ]
  then
     echo "Input file needs to be specified with the flag -I"
     echo
     echo -e '\nUsage: spadesbatch -I INPUT_FOLDER -O OUTPUT_PREFIX [-h] [... OTHER_OPTIONAL_ARGUMENTS] \n \nRun this program with flag -h for full usage options\n'
     exit
  else
     CleanReadsFolder=$(realpath $CleanReadsFolder)
     echo "Input file "$CleanReadsFolder""
fi

if [ "$OutputFolder" == "" ]
  then
     OutputFolder=$(pwd)
	 echo "Output folder "$OutputFolder""
  else
     [ -d $OutputFolder ] && echo "Output directory exists" && exit
     mkdir $OutputFolder
     OutputFolder=$(realpath $OutputFolder)
     echo "Output folder name "$OutputFolder""
fi


outfolder=$OutputFolder/spades_assemb
[ -d $outfolder ] && echo "A folder with spades assemblies already exists in the output directory" && exit
mkdir $outfolder &&


mkdir -p $OutputFolder/logs

echo '
###################################
Batch Assembly with SPAdes Version '$version'

Created by Guilherme Azevedo 2021 
###################################

This program runs SPAdes in batch mode.

Please cite SPAdes:
Prjibelski, A., Antipov, D., Meleshko, D., Lapidus, A., & Korobeynikov, A. (2020). Using SPAdes de novo assembler. Current Protocols in Bioinformatics, 70, e102. doi: 10.1002/cpbi.102

----------------------------------------------------------------------
'$(date)'

Input folder: '$CleanReadsFolder'
Output folder: '$OutputFolder'

Start assembling with SPAdes using the commands:
spades.py '$spadesflag' --careful --cov-cutoff auto --memory '$memory' -t '$nthreads'

' >> $OutputFolder/logs/progress.txt

progressfile=$OutputFolder/logs/progress.txt

speciesfolders=$(ls $CleanReadsFolder)

spadesflag="--"$spadesoptions


for i in $speciesfolders
do
echo "Processing "$i"" >> $progressfile
sem --jobs $n_jobs \
spades.py $spadesflag --careful --cov-cutoff auto --memory $memory -t $nthreads \
--pe1-1 $CleanReadsFolder/$i/split-adapter-quality-trimmed/*-READ1.fastq.gz \
--pe1-2 $CleanReadsFolder/$i/split-adapter-quality-trimmed/*-READ2.fastq.gz \
--pe1-s $CleanReadsFolder/$i/split-adapter-quality-trimmed/*-READ-singleton.fastq.gz \
--pe1-$orientation \
-o $outfolder/$i
done
sem --wait

resultsfolder=$(ls $outfolder)


finalcontigsfolder=$OutputFolder/spades_contigs
[ -d $finalcontigsfolder ] && echo "A folder with spades contigs already exists in the output directory. The contigs were not coppied to a single folder, but you can find them in the spades_assemblies folder" && exit
mkdir $finalcontigsfolder

for i in $resultsfolder
do
  sem --jobs $n_jobs cp $outfolder/$i/contigs.fasta $finalcontigsfolder/$i.fasta
done
sem --wait

# Check for warnings in SPAdes

for i in $(ls $resultsfolder)
  do
    warning=$(ls $i | grep -c "warnings.log")
	if [ "$warning" != 0 ]
	  then
	    echo $i >> $OutputFolder/logs/CheckSpadesWarnings.txt
		echo "There was warnings in SPAdes assemblies. Check the SPAdes logs" >> $progressfile
	fi
done


echo 'Finished assembling with spades


'$(date)'
----------------------------------------------------------------------
######################################################################

' >> $progressfile
