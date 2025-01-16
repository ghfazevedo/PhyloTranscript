#!/bin/bash

get_fasta_from_bam () {
	input_folder=$(realpath "$1") 
    [ -d "fasta_from_bam_by_species" ] && echo "Output directory exists" && exit
    mkdir fasta_from_bam_by_species
    mkdir fasta_from_bam_by_gene

    for bamfile in $input_folder/*.bam
      do
         name=$(basename $bamfile .bam)
         samtools consensus --min-MQ 30 -d 2 -@ 20 -o fasta_from_bam_by_species/$name".fasta"  $bamfile
    done

    for fastafile in fasta_from_bam_by_species/*.fasta
      do
        species_name=$(basename $fastafile .fasta)
        seq_headers=$(grep ">" $fastafile)
        for header in $seq_headers
          do
            gene_name=$(echo $header | cut -f2 -d">" )
            sed '$a>' $fastafile | sed -n "/$header/,/>/p" | head -n -1 | sed "s/$header/>$species_name |$gene_name/g"  >> fasta_from_bam_by_gene/$gene_name".fasta"
        done
    done
}

get_fasta_from_bam "$1"

#"$1" is the folder with bam files