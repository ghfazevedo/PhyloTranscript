#!/bin/bash

add_sp_names_to_fasta_headers () {
        contigs_folder_renamed="$1""_renamed"
        [ -d "$contigs_folder_renamed" ] && echo "A folder with renamed contigs already exists in the output directory" && exit
		echo "Copying transcript files"
        cp -r "$1" $contigs_folder_renamed
        cd $contigs_folder_renamed
        files=$(ls)
        
        echo "Renaming sequences"
        for file in $files
          do
            echo "Processing $file"
            name=$(echo $file | cut -f1 -d".")
            string_to_add=">"$name" | "
            sed -i "s/>/$string_to_add/g" $file
          done
}

add_sp_names_to_fasta_headers "$1"