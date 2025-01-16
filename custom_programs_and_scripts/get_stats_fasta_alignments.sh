#!/bin/bash

get_stats_fasta_alignments () {
    mkdir -p logs
    phyluce_align_get_taxon_locus_counts_in_alignments \
        --alignments "$1" \
        --input-format fasta \
        --output "$1""_taxon_count".csv \
        --log-path logs
    
    phyluce_align_get_align_summary_data \
        --alignments "$1" \
        --input-format fasta \
        --log-path logs \
        --cores 12 \
        --output-stats "$1""_stats".csv 
}

get_stats_fasta_alignments "$1"

# The argument "$1" is a folder with alignments