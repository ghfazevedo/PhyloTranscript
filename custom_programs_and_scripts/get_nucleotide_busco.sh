#!/bin/bash
version=0.2
# Version Log

print_message () {
  echo '
#################################################
Get Nucleotide BUSCO Sequences Version '$version'

Created by Guilherme Azevedo 2021 
#################################################

This program retrieves the single copy contigs from each species and combine them in one fasta file per gene in the output folder called *busco_nucleotide_seqs*. It takes as input a folder with the contigs in fasta format, a folder with the output generated from a BUSCO run and the BUSCO lineage dataset used in the BUSCO run.

It requires seqtk(https://github.com/lh3/seqtk) v1.3-r106 and util-linux(https://man7.org/linux/man-pages/man1/rename.1.html) v.2.36 to be installed in your system. Optionally you can provaide a prefix for the output folder (the default is busco_nucleotide_seqs).

Usage:
 get_nucleotide_busco -C contigs_folder -B busco_output -D busco_database [-o output_prefix] [-h] [... OTHER_OPTIONAL_ARGUMENTS]
     
     -C contigs_folder    Folder with the contigs, one per species,
                              in fasta format
     -B busco_output      Folder with BUSCO output contaning
                              one folder per species
     -D busco_dataset    The name of the dataset used 
                             (e.g. arachnida_odb10).
     -o ouput_prefix      The prefix for the output folder name 
                             (Optional. Default: busco_nucleotide_seqs
     -f [yes,no]         If "yes", it retrieves the BUSCO that were found 
                             to be fragmented together with the complete BUSCO. If "no", only the complete BUSCO are used.
                             (Default="yes")
     -h                   Print this message and exit
'
exit
}



#Set variables and default values
contigs_folder=""
busco_output=""
busco_dataset=""
ouput_prefix=""
fragmented_busco="yes"


while getopts "C:B:D:o:f:h" flag; do
    case "${flag}" in
            C)  contigs_folder=${OPTARG};;
            B)  busco_output=${OPTARG};;
            D)  busco_dataset=${OPTARG};;
            o)  ouput_prefix=${OPTARG};;
            f)  fragmented_busco=${OPTARG};;
            h)  print_message;;
            ?)  print_message;;
    esac
done



if [ "$contigs_folder" == "" ]
  then
     echo "Contigs folder needs to be specified with the flag -C"
     echo
     print_message
     exit
  else
     contigs_folder=$(realpath $contigs_folder)
     echo "Contigs folder "$contigs_folder""
fi

if [ "$busco_output" == "" ]
  then
     echo "The BUSCO output folder needs to be specified with the flag -B"
     echo
     print_message
     exit
  else
     busco_output=$(realpath $busco_output)
     echo "BUSCO folder "$busco_output""
fi

if [ "$busco_dataset" == "" ]
  then
     echo "The BUSCO dataset needs to be specified with the flag -D"
     echo
     print_message
     exit
  else
     echo "BUSCO dataset "$busco_dataset""
     busco_dataset="run_"$busco_dataset
fi

if [ "$ouput_prefix" == "" ]
  then
     ouput_prefix=busco_nucleotide_seqs
     [ -d $ouput_prefix ] && echo "A folder with busco nucleotides sequences already exists in the output directory" && exit
     mkdir $ouput_prefix
     ouput_prefix=$(realpath $ouput_prefix)
  else
     [ -d $ouput_prefix ] && echo "A folder with busco nucleotides sequences already exists in the output directory" && exit
     mkdir $ouput_prefix
     ouput_prefix==$(realpath $ouput_prefix)
fi


mkdir -p logs
logs_path=$(realpath logs)

echo '
#################################################
Get Nucleotide BUSCO Sequences Version '$version'

Created by Guilherme Azevedo 2021 
#################################################

This program retrieves the single copy contigs from each species and combine them in one fasta file per gene

----------------------------------------------------------------------
'$(date)'

Start retrieving nucleotide sequences from BUSCO outputs.

' >> $logs_path/progress.txt

sp_list=$(ls -d $busco_output/*.fasta | rev | cut -d'/' -f 1 | rev)


# Retrieve complete BUSCO function
retrieve_complete_BUSCO () {
for species in $sp_list
  do
    echo "Retrieving complete BUSCO for $species" | tee -a $logs_path/progress.txt
    sp_name=$(echo $species | cut -f1 -d".")
    cp -r $busco_output/$species/$busco_dataset/busco_sequences/single_copy_busco_sequences $busco_output/$species/$busco_dataset/busco_sequences/single_copy_busco_sequences_backup
    cd $busco_output/$species/$busco_dataset/busco_sequences/single_copy_busco_sequences_backup
    for file in *.*
      do
        contig_name=$(grep ">" $file | cut -f1 -d":" | cut -f2 -d">" )
        echo $contig_name > $ouput_prefix/contig_name_temp.txt
        seqtk subseq $contigs_folder/$species $ouput_prefix/contig_name_temp.txt >> $ouput_prefix/$file
        sed -i "s/$contig_name/$sp_name/g" $ouput_prefix/$file
      done
    cd $busco_output
  done
  
echo "Removing temporary files" | tee -a $logs_path/progress.txt
rm $ouput_prefix/contig_name_temp.txt
for species in $sp_list
  do
    rm -r $busco_output/$species/$busco_dataset/busco_sequences/single_copy_busco_sequences_backup
    
  done
}

  
# Retrieve fragmented BUSCO function
retrieve_fragmented_BUSCO () {
    for species in $sp_list
      do
        echo "Retrieving fragmented BUSCO for $species" | tee -a $logs_path/progress.txt
        sp_name=$(echo $species | cut -f1 -d".")
        cp -r $busco_output/$species/$busco_dataset/busco_sequences/fragmented_busco_sequences $busco_output/$species/$busco_dataset/busco_sequences/fragmented_busco_sequences_backup
        cd $busco_output/$species/$busco_dataset/busco_sequences/fragmented_busco_sequences_backup
        for file in *.*
          do
            contig_name=$(grep ">" $file | cut -f1 -d":" | cut -f2 -d">" )
            echo $contig_name > $ouput_prefix/contig_name_temp.txt
            seqtk subseq $contigs_folder/$species $ouput_prefix/contig_name_temp.txt >> $ouput_prefix/$file
            sed -i "s/$contig_name/$sp_name/g" $ouput_prefix/$file
          done
        cd $busco_output
      done
      
    echo "Removing temporary files" | tee -a $logs_path/progress.txt
    rm $ouput_prefix/contig_name_temp.txt
    for species in $sp_list
      do
        rm -r $busco_output/$species/$busco_dataset/busco_sequences/fragmented_busco_sequences_backup
      done 
      
}


# Retrieve multicopy BUSCO function
retrieve_multicopy_BUSCO () {
    echo "Locus,Species,Number_of_Contigs,Number_of_Unique_Contigs" > duplicated_genes_with_concat_copies_info.txt
    path_to_txt_with_concatenated_copies=$(realpath duplicated_genes_with_concat_copies_info.txt)
    touch duplicated_genes_with_concat_copies_unique_list.txt
    path_to_unique_list=$(realpath duplicated_genes_with_concat_copies_unique_list.txt)
    for species in $sp_list
      do
        echo "Retrieving multi copy BUSCO for $species" | tee -a $logs_path/progress.txt
        sp_name=$(echo $species | cut -f1 -d".")
        cp -r $busco_output/$species/$busco_dataset/busco_sequences/multi_copy_busco_sequences $busco_output/$species/$busco_dataset/busco_sequences/multi_copy_busco_sequences_backup
        cd $busco_output/$species/$busco_dataset/busco_sequences/multi_copy_busco_sequences_backup
        for file in *.*
          do
            grep ">" $file | cut -f1 -d":" | cut -f2 -d">" > $ouput_prefix/contig_name_temp.txt
            
            lines_before=$(cat $ouput_prefix/contig_name_temp.txt | wc -l )
            sort -u $ouput_prefix/contig_name_temp.txt -o $ouput_prefix/contig_name_temp.txt
            lines_after=$(cat $ouput_prefix/contig_name_temp.txt | wc -l )
            
            seqtk subseq $contigs_folder/$species $ouput_prefix/contig_name_temp.txt >> $ouput_prefix/$file
            line_number=1
            if [ "$lines_after" = "$lines_before" ]
              then
                while read -r line
                  do    
                    sed -i "s/$line/$sp_name\-genecopy0$line_number/g" $ouput_prefix/$file
                    line_number=$((line_number+1))
                done < $ouput_prefix/contig_name_temp.txt
              else
                echo $file,$sp_name,$lines_before,$lines_after >> $path_to_txt_with_concatenated_copies
                echo $file >> $path_to_unique_list
                while read -r line
                  do    
                    sed -i "s/$line/$sp_name\-genecopyconcat0$line_number/g" $ouput_prefix/$file
                    line_number=$((line_number+1))
                done < $ouput_prefix/contig_name_temp.txt
            fi    
          done
        cd $busco_output
      done
    sort -u $path_to_unique_list -o $path_to_unique_list  
    echo "Removing temporary files" | tee -a $logs_path/progress.txt
    rm $ouput_prefix/contig_name_temp.txt
    for species in $sp_list
      do
        rm -r $busco_output/$species/$busco_dataset/busco_sequences/multi_copy_busco_sequences_backup
      done 
      
}


# Separate multi cpoy from single copy in the dataset
separate_sinlge_and_multi_copy_BUSCO () {
    echo "Separating single copy from multi copy genes" | tee -a $logs_path/progress.txt
    cp -r $ouput_prefix $ouput_prefix"_singlecopy"
    mkdir $ouput_prefix"_multicopy"
    grep -l -d "recurse" "genecopy" $ouput_prefix"_singlecopy" > list_of_multi_copy_genes.txt
    while read -r line
      do
        mv $line $ouput_prefix"_multicopy"
    done < list_of_multi_copy_genes.txt
}


if [ "$fragmented_busco" == "yes" ]
  then
     echo "Both fragmented and complete BUSCO will be retrieved" | tee -a $logs_path/progress.txt
     retrieve_complete_BUSCO &&
     retrieve_fragmented_BUSCO &&
     retrieve_multicopy_BUSCO &&
     separate_sinlge_and_multi_copy_BUSCO
  else
     echo "Only complete BUSCO will be retrieved" | tee -a $logs_path/progress.txt
     retrieve_complete_BUSCO &&
     retrieve_multicopy_BUSCO &&
     separate_sinlge_and_multi_copy_BUSCO
fi

rename .faa .fasta $ouput_prefix/*.faa
rename .faa .fasta $ouput_prefix"_multicopy"/*.faa
rename .faa .fasta $ouput_prefix"_singlecopy"/*.faa

echo 'Finished retrieving nucleotide sequences from BUSCO outputs.

'$(date)'
----------------------------------------------------------------------
######################################################################

' | tee -a $logs_path/progress.txt
