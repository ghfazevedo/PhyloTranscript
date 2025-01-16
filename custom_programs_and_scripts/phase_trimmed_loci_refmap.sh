#!/bin/bash
version=0.2
# Version Log
# v0.2 added a working directory variable so it can be run in any directory.

print_message () {
  echo '
###################################
phase_trimmed_loci Version '$version'

Created by Guilherme Azevedo 2021 
###################################

This program runs phyluce workflow for phasing UCE data mapping to a reference.
You only need to specify the folder with clean reads and the reference fasta file.
The script will automatically generate the configuration files needed.

Usage:

phase_trimmed_loci_refmap.s -R CleanReadsFolder -O OutputPrefix -C Reference -e CondaEnv [-h] [... OTHER_OPTIONAL_ARGUMENTS] 

      -R  CleanReadsFolder    The folder with the clean reads (trimmomatic output)
      -O  OutputPrefix        The name of the folder in whcih the final assemblies will be saved
      -C  Reference           The reference fasta file
      -e  CondaEnvironment    The conda environment with phyluce installed
      -n  n_cores             Number of cores. Default: 16
      -h  help

See phyluce workflow (https://phyluce.readthedocs.io/en/latest/daily-use/daily-use-4-workflows.html#) web page for more information

All default parameters are the same as used in phyluce v1.7.1.

'
exit
}



#path to the clean reads (relative or absolute)
cleanreads=""


#path to reference file
reference=""


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
            C)  reference=${OPTARG};;
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


if [ "$reference" == "" ]
  then
     echo "Reference file name needs to be specified wit flag -C"
     echo
     print_message
     exit
  else
     echo "Reference file "$reference""
fi

if [ "$condaenv" == "" ]
  then
     echo "Conda environment not specified. Using base"
  else
     source activate $condaenv
fi

if [ "$outputprefix" == "" ]
  then
     outputprefix=output
	 [ -d $outputprefix"_phased_refmap_contigs" ] && echo "A folder with phased contigs already exists in the output directory" && exit
	 [ -d $outputprefix"_phased_refmap_loci" ] && echo "A folder with phased loci already exists in the output directory" && exit
 	 [ -d $outputprefix"_phased_refmap_output" ] && echo "A folder with phased outputs already exists in the output directory" && exit
 	 [ -d $outputprefix"_mapped_refmap" ] && echo "A folder with mapped output already exists in the output directory" && exit
     [ -d $outputprefix"_phased_refmap_vcfs"] && echo "A folder with vcf files already exists in the output directory" && exit
     [ -d $outputprefix"_phased_refmap_bams"] && echo "A folder with bam files already exists in the output directory" && exit

  else
	 [ -d $outputprefix"_phased_refmap_contigs" ] && echo "A folder with phased contigs already exists in the output directory" && exit
	 [ -d $outputprefix"_phased_refmap_loci" ] && echo "A folder with phased loci already exists in the output directory" && exit
 	 [ -d $outputprefix"_phased_refmap_output" ] && echo "A folder with phased outputs already exists in the output directory" && exit
 	 [ -d $outputprefix"_mapped_refmap" ] && echo "A folder with mapped output already exists in the output directory" && exitfi
     [ -d $outputprefix"_phased_refmap_vcfs"] && echo "A folder with vcf files already exists in the output directory" && exit
     [ -d $outputprefix"_phased_refmap_bams"] && echo "A folder with bam files already exists in the output directory" && exit

fi

mkdir -p logs
logs_path=$(realpath logs)

cleanreadspath=$(realpath $cleanreads)

referencepath=$(realpath $reference)

echo '
###################################
phase_trimmed_loci Version '$version'

Created by Guilherme Azevedo 2021 
###################################

This program runs phyluce workflow for phasing UCE data.
You only need to specify the folder with clean reads and the folder with the reference.
The script will automatically generate the configuration files needed.

----------------------------------------------------------------------
'$(date)'

Start phasing.

' | tee -a $logs_path/progress.txt


### Create conf file

echo 'Creating configuration file for mapping' | tee -a $logs_path/progress.txt

cd $cleanreadspath
echo "reads:" > $workingdirectory/phase_wfA.conf
speciesfolder=$(ls)
for i in $speciesfolder
   do 
     echo "    "$i": "$cleanreadspath"/"$i"/split-adapter-quality-trimmed/"
   done >> $workingdirectory/phase_wfA.conf
echo -en "\ncontigs:\n" > $workingdirectory/phase_wfB.conf
for i in $speciesfolder
   do 
     echo "    "$i": "$referencepath
   done >> $workingdirectory/phase_wfB.conf
cat $workingdirectory/phase_wfA.conf $workingdirectory/phase_wfB.conf > $workingdirectory/"config_mapping_refmap.yaml"
rm  $workingdirectory/phase_wfA.conf
rm $workingdirectory/phase_wfB.conf 
cd $workingdirectory

### Run workflow for mapping
echo 'Runnig workflow for mapping' | tee -a $logs_path/progress.txt

phyluce_workflow \
    --config config_mapping_refmap.yaml \
    --output $outputprefix"_mapped_refmap" \
    --workflow mapping \
    --cores $n_cores

#prepare configuration file for phasing
#path to bams

echo 'Creating configuration file for phasing' | tee -a $logs_path/progress.txt

mappedpath=$(realpath $outputprefix"_mapped_refmap")

cd $cleanreadspath
echo "bams:" > $workingdirectory/phase_wfC.conf
speciesfolder=$(ls)
for i in $speciesfolder
   do 
     echo "    "$i": "$mappedpath"/mapped_reads/"$i".fxm.sorted.md.bam"
   done >> $workingdirectory/phase_wfC.conf
echo -en "\ncontigs:\n" > $workingdirectory/phase_wfD.conf
for i in *
   do 
     echo "    "$i": "$referencepath
   done >> $workingdirectory/phase_wfD.conf
cat $workingdirectory/phase_wfC.conf $workingdirectory/phase_wfD.conf > $workingdirectory/"config_phase_refmap.yaml"
rm  $workingdirectory/phase_wfC.conf; 
rm $workingdirectory/phase_wfD.conf 
cd $workingdirectory

#run workflow
echo 'Runnig modified workflow for phasing' | tee -a $logs_path/progress.txt

phyluce_workflow \
	--config config_phase_refmap.yaml \
    --output $outputprefix"_phased_refmap_output" \
    --workflow phasingtrimmed \
    --cores $n_cores

mkdir $outputprefix"_phased_refmap_contigs"
cd $outputprefix"_phased_refmap_contigs"
cp $workingdirectory/$outputprefix"_phased_refmap_output"/fastas/*.fasta .
cd $workingdirectory


#################################################################################################################
# This part renames the fasta file removing the orthogroup/BUSCO name from the seqname and adds it to the end:
# input:  >0000_Species_name_code_pilon
# output: >Species_name_code_phase |0000 
#################################################################################################################

echo 'Renaming fasta headers' | tee -a $logs_path/progress.txt

cd $outputprefix"_phased_refmap_contigs"
filelist=$(ls)

for file in $filelist
  do
    echo "Renaming "$file""
    seqnames=$(grep ">" $file)
    phase=$(echo $file | cut -f2 -d".")
    phase="p"$phase
    for i in $seqnames
        do
          lociname=$(echo $i | cut -f1 -d"_" | sed "s/>//g" )
          tobereplaced=$lociname"_"
          spname=$(echo $i | sed "s/$tobereplaced//g")
          newname=$(echo $spname" |"$lociname)
          sed -i "s/$i/$newname/g" $file
        done
    sed -i "s/pilon/$phase/g" $file
  done

cd $workingdirectory

# Make a folder with one file per locus

mkdir $outputprefix"_phased_refmap_loci"
loci_folder=$(realpath $outputprefix"_phased_refmap_loci")

cd $outputprefix"_phased_refmap_contigs"
filelist=$(ls)
for file in $filelist
  do
    modifiedseqnames=$(grep ">" $file | sed "s/ |/|/g")
	for i in $modifiedseqnames
        do
		  sample_name=$(echo $i | cut -f1 -d"|")
		  locus=$(echo $i | cut -f2 -d"|")
		  sed '$a>' $file | sed -n "/$sample_name |$locus/,/>/p" | head -n -1 >> $loci_folder/$locus".fasta"
        done
  done
cd $workingdirectory

# Make a folder with one vcf files per species and the bam files

mkdir $outputprefix"_phased_refmap_vcfs"
cd $outputprefix"_phased_refmap_vcfs"
cp $workingdirectory/$outputprefix"_phased_refmap_output"/fastas/*.vcf .
cd $workingdirectory

mkdir $outputprefix"_phased_refmap_bams"
cd $outputprefix"_phased_refmap_bams"
cp $workingdirectory/$outputprefix"_phased_refmap_output"/bams/*.bam .
cp $workingdirectory/$outputprefix"_phased_refmap_output"/bams/*.bai .
rm *.chimera.bam
cd $workingdirectory


echo 'Finished phasing.



'$(date)'
----------------------------------------------------------------------
######################################################################

' | tee -a $logs_path/progress.txt








