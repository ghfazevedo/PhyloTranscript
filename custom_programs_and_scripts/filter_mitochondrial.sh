#!/bin/bash


filter_mitochondrial(){
    pseudoreference_file=$(realpath "$1")
    mitogenome_fasta_file=$(realpath "$2")
    folder_with_alignments=$(realpath "$3")
    output_folder_nuc=$folder_with_alignments"_nuclear"
    [ -d $output_folder_nuc ] && echo "Output folder exists" && exit
    output_folder_mito=$folder_with_alignments"_mito"
    [ -d $output_folder_mito ] && echo "Output folder exists" && exit
    
    
    mkdir -p logs
    logs_path=$(realpath logs)
    echo "Filtering Mitochondrial genome " "$(date)" | tee -a $logs_path/progress.txt
    
    wroking_directory=$PWD
    
    blastn -query $pseudoreference_file -subject $mitogenome_fasta_file  -outfmt 6 -max_target_seqs 1 -out mitochondria_blastresult.txt
    
    cp -r $folder_with_alignments $output_folder_nuc
    
    mkdir $output_folder_mito
    
    while read -r line
      do
        file=$(echo $line | cut -f1 -d" ")
        mv $output_folder_nuc/$file $output_folder_mito/$file
    done < mitochondria_blastresult.txt
}


filter_mitochondrial "$1" "$2" "$3"

# "$1" is the pseudoreference_file, "$2" is the mitochondrial genome fasta file and "$3"is the folder with all alignments
