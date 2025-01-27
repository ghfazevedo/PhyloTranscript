#!/bin/bash
version=0.2
# Version Log
# v0.2 added a working directory variable so it can be run in any directory.

print_message () {
  echo '
###################################
call_snps_using_reference '$version'

Created by Guilherme Azevedo 2021 
###################################

This program runs phyluce workflow for mapping reads against a reference sequence and then runs a customized workflow for calling snps.
You only need to specify the folder with clean reads and the folder with the consensus_reference.
The script will automatically generate the configuration files needed.
The snp calling workflow is available from a fork of original phyluce at https://github.com/ghfazevedo/phyluce.git

Usage:

phase_trimmed_loci.sh -R CleanReadsFolder -O OutputPrefix -C ConsensusReference -e CondaEnv [-h] [... OTHER_OPTIONAL_ARGUMENTS] 

      -R  CleanReadsFolder    The folder with the clean reads (trimmomatic output)
      -O  OutputPrefix        The name of the folder in whcih the final assemblies will be saved
      -C  ConsensusReference  The folder with pseudoreference
      -e  CondaEnvironment    The conda environment with phyluce installed
      -n  n_cores             Number of cores. Default: 16
      -h  help

See phyluce workflow (https://phyluce.readthedocs.io/en/latest/daily-use/daily-use-4-workflows.html#) web page for more information and https://github.com/ghfazevedo/phyluce.git.
All default parameters are the same as used in phyluce v1.7.1.

If you use this program, cite phyluce and refer to the customized workflow fork.

'
exit
}



#path to the clean reads (relative or absolute)
cleanreads=""


#path to consensus_reference folder
consensus_reference=""


#choose output prefix
outputprefix=""

#choose the number of cores
n_cores=16

#chose conda enviroment
condaenv=""

# Set working directory
workingdirectory=$(pwd)

while getopts "R:O:C:e:n:h" flag; do
    case "${flag}" in
            R)  cleanreads=${OPTARG};;
            O)  outputprefix=${OPTARG};;
            C)  consensus_reference=${OPTARG};;
			e)  condaenv=${OPTARG};;
			n)  n_cores=${OPTARG};;
            h)  print_message;;
            ?)  print_message;;
    esac
done

if [ "$cleanreads" == "" ]
  then
     echo "Folder with clean reads file needs to be specified with the flag -R"
     echo
     print_message
     exit
  else
     echo "Clean reads folder "$cleanreads""
fi


if [ "$consensus_reference" == "" ]
  then
     echo "consensus_reference folder name needs to be specified wit flag -C"
     echo
     print_message
     exit
  else
     echo "consensus_reference folder "$consensus_reference""
fi

if [ "$condaenv" == "" ]
  then
     echo "Conda environment not specified. Using base"
  else
     source activate $condaenv
fi

if [ "$outputprefix" == "" ]
  then
     outputprefix=$consensus_reference
	 [ -d $outputprefix"_snps" ] && echo "A folder with snps already exists in the output directory" && exit
	 [ -d $outputprefix"_ref_mapped" ] && echo "A folder with mapped output already exists in the output directory" && exit
    
  else
	 [ -d $outputprefix"_snps" ] && echo "A folder with snps already exists in the output directory" && exit
 	 [ -d $outputprefix"_ref_mapped" ] && echo "A folder with mapped output already exists in the output directory" && exitfi
fi

mkdir -p logs
logs_path=$(realpath logs)

cleanreadspath=$(realpath $cleanreads)

consensus_referencepath=$(realpath $consensus_reference)

echo '
###################################
call_snps_using_reference '$version'

Created by Guilherme Azevedo 2021 
###################################

This program runs phyluce workflow for phasing UCE data.
You only need to specify the folder with clean reads and the folder with the consensus_reference.
The script will automatically generate the configuration files needed.

----------------------------------------------------------------------
'$(date)'

Start phasing.

' | tee -a $logs_path/progress.txt


### Create conf file

echo 'Creating configuration file for mapping' | tee -a $logs_path/progress.txt

cd $cleanreadspath
echo "reads:" > $workingdirectory/snps_wfA.conf
speciesfolder=$(ls)
for i in $speciesfolder
   do 
     echo "    "$i": "$cleanreadspath"/"$i"/split-adapter-quality-trimmed/"
   done >> $workingdirectory/snps_wfA.conf
echo -en "\ncontigs:\n" > $workingdirectory/snps_wfB.conf
for i in $speciesfolder
   do 
     echo "    "$i": "$consensus_referencepath
   done >> $workingdirectory/snps_wfB.conf
cat $workingdirectory/snps_wfA.conf $workingdirectory/snps_wfB.conf > $workingdirectory/"config_ref_mapping.yaml"
rm  $workingdirectory/snps_wfA.conf
rm $workingdirectory/snps_wfB.conf 
cd $workingdirectory

### Run workflow for mapping
echo 'Runnig workflow for mapping against reference' | tee -a $logs_path/progress.txt

phyluce_workflow \
    --config config_ref_mapping.yaml \
    --output $outputprefix"_ref_mapped" \
    --workflow mapping \
    --cores $n_cores

#prepare configuration file for snp calling
#path to bams

echo 'Creating configuration file for calling' | tee -a $logs_path/progress.txt

mappedpath=$(realpath $outputprefix"_ref_mapped")

cd $cleanreadspath
echo "bams:" > $workingdirectory/snps_wfC.conf
speciesfolder=$(ls)
for i in $speciesfolder
   do 
     echo "    "$i": "$mappedpath"/mapped_reads/"$i".fxm.sorted.md.bam"
   done >> $workingdirectory/snps_wfC.conf
echo -en "\nreference:\n" > $workingdirectory/snps_wfD.conf
echo "    "$consensus_referencepath >> $workingdirectory/snps_wfD.conf
cat $workingdirectory/snps_wfC.conf $workingdirectory/snps_wfD.conf > $workingdirectory/"config_snps.yaml"
rm  $workingdirectory/snps_wfC.conf; 
rm $workingdirectory/snps_wfD.conf 
cd $workingdirectory

#run workflow
echo 'Runnig workflow for calling snps' | tee -a $logs_path/progress.txt

phyluce_workflow \
	--config config_snps.yaml \
    --output $outputprefix"_snps" \
    --workflow vcfsnps \
    --cores $n_cores

bgzip $outputprefix"_snps"/calls/variants.vcf
bgzip $outputprefix"_snps"/calls/biallelic.vcf.recode.vcf



echo 'Finished calling snps.

'$(date)'
----------------------------------------------------------------------
######################################################################

' | tee -a $logs_path/progress.txt








