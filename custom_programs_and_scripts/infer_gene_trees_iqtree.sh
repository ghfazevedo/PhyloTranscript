#!/bin/bash

infer_gene_trees_iqtree () {
   if [ "$1" = "" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ] || [ "$1" = "help" ]
    then
       echo "Usage: infer_gene_trees_iqtree.sh FOLDER_WITH_ALIGNMENTS NUMBER_OF_THREADS NUMBER_OF_PARALLEL_JOBS"
       echo "Folder with alignments must be given as first argument."
       echo "NUMBER_OF_THREADS and NUMBER_OF_PARALLEL_JOBS must be given as second and third arguments, respectively"
       exit
    else inputfolder=$(realpath "$1")
   fi  
   outputprefix=$(basename "$1")
   outputfolder=$outputprefix"_genetrees"
   threads="$2"
   n_jobs="$3"

   [ -d $outputfolder ] && echo "Output folder exists" && exit
   
   
   mkdir -p logs
   echo 'Inferring gene trees with IQTREE '$(date)'' | tee -a logs/progress.txt
   
   mkdir $outputfolder
   
   
   for i in $inputfolder/*.*
     do
        sem --id iqtreejob --jobs $n_jobs iqtree -s $i -m MFP -mset mrbayes --mrate I,G,I+G -T $threads -B 1000 --prefix $outputfolder"/"$(basename $i)
    done
   sem --wait --id iqtreejob
   
   for i in $outputfolder/*.treefile
     do
       paste $i >> $outputprefix"_gene_trees.nwck"
       echo $(basename $i) >> $outputprefix"_gene_tree_order.txt"
    done

 
   rm -r $outputfolder/*.gz 
   rm -r $outputfolder/*.phy 
   rm -r $outputfolder/*.bionj 
   rm -r $outputfolder/*.mldist 
   rm -r $outputfolder/*.contree 
   rm -r $outputfolder/*.nex 
   rm -r $outputfolder/*.iqtree
   

   echo 'Finished inferring gene trees with IQTREE '$(date)' ' | tee -a logs/progress.txt
}

infer_gene_trees_iqtree "$1" "$2" "$3"

# "$1" is the folder with alignments
# "$2" is the number of threads for IQTREE
# "$3" is the number of parallel jobs for parallel