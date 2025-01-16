#!/bin/bash

infer_gene_trees_iqtree_quick (){
   inputfolder=$(realpath "$1")
   outputfolder=$inputfolder"_genetrees"
   threads="$2"
   if [ $threads = "" ]
    then
       threads="AUTO"
    fi
   [ -d $outputfolder ] && echo "Output folder exists" && exit
   mkdir -p logs
   echo "Inferring gene trees with IQTREE" | tee -a logs/progress.txt
   mkdir $outputfolder

   
   cp -r $inputfolder $outputfolder

   cd $outputfolder
   
   for i in *.fa*
     do
        iqtree -s $i -m LG -T $threads 
     done
   
   for i in *.treefile
     do
       paste $i >> gene_trees.nwck
       echo $i >> gene_tree_order.txt
     done
     
 
   rm *.fasta *.faa *.fan *.fas *.fa

   echo "Finished inferring gene trees with IQTREE" | tee -a logs/progress.txt
}

infer_gene_trees_iqtree_quick "$1"

# "$1" is the folder with alignments