#!/bin/bash

combine_busco_amino () {

   [ -d busco_aminoacids_seqs ] && echo "Output folder exists" && exit
   mkdir -p logs
   echo "Combining files with BUSCO aminoacid sequences" | tee -a logs/progress.txt
   mkdir busco_aminoacid_seqs
   out_folder=$(realpath busco_aminoacid_seqs)
   busco_results_folder=$(realpath "$1")
   folder_list=$(ls -d $busco_results_folder/*.fasta)
   n_species=$(ls -1 -d $busco_results_folder/*.fasta | wc -l)
   
# this script copies only the single copy and complete
   for species_folder in $folder_list
      do
        species_name=$(echo $species_folder | rev | cut -d'/' -f 1 | rev | cut -f1 -d".")
        cd $species_folder"/run_arachnida_odb10/busco_sequences"
        cp -r single_copy_busco_sequences single_copy_busco_sequences_renamed
        cd single_copy_busco_sequences_renamed
        file_list=$(ls *.faa)
        for i in $file_list;       
          do
            to_be_replaced=$(grep ">" $i)
            new_name=">"$species_name
            sed -i "s/$to_be_replaced/$new_name/g" $i
        done
        for i in $file_list
          do
            cat $i >> $out_folder/$i
        done
        cd ../../../..
    done
    #cp -r $out_folder $out_folder"_only_complete_BUSCO"
   
# this script copies fragmented single copy
   for species_folder in $folder_list
      do
        species_name=$(echo $species_folder | rev | cut -d'/' -f 1 | rev | cut -f1 -d".")
        cd $species_folder"/run_arachnida_odb10/busco_sequences"
        cp -r fragmented_busco_sequences fragmented_busco_sequences_renamed
        cd fragmented_busco_sequences_renamed
        file_list=$(ls *.faa)
        for i in $file_list;       
          do
            to_be_replaced=$(grep ">" $i)
            new_name=">"$species_name
            sed -i "s/$to_be_replaced/$new_name/g" $i
        done
        for i in $file_list
          do
            cat $i >> $out_folder/$i
        done
        cd ../../../..
    done



# this script copies multi copy
   for species_folder in $folder_list
      do
        species_name=$(echo $species_folder | rev | cut -d'/' -f 1 | rev | cut -f1 -d".")
        cd $species_folder"/run_arachnida_odb10/busco_sequences"
        cp -r multi_copy_busco_sequences multi_copy_busco_sequences_renamed
        cd multi_copy_busco_sequences_renamed
        file_list=$(ls *.faa)
        for i in $file_list;       
          do
            to_be_replaced=$(grep ">" $i)
            new_name=">"$species_name
			copy_number=1
			for name in $to_be_replaced
			  do
			    sed -i "s/$name/$new_name\-genecopy0$copy_number/g" $i
				copy_number=$((copy_number+1))
			done
        done
        for i in $file_list
          do
            cat $i >> $out_folder/$i
        done
        cd ../../../..
    done


# this part separates single and multi copy BUSCO
    echo "Separating single copy from multi copy genes" 
    cp -r $out_folder $out_folder"_singlecopy"
    mkdir $out_folder"_multicopy"
    grep -l -d "recurse" "genecopy" $out_folder"_singlecopy" > list_of_multi_copy_genes_aa.txt
    while read -r line
      do
        mv $line $out_folder"_multicopy"
    done < list_of_multi_copy_genes_aa.txt
}



combine_busco_amino "$1"

# "$1" is the folder with the busco outputs 