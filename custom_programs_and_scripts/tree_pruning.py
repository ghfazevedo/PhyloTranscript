#! /usr/bin/env python
# -*- coding: utf-8 -*-


############################################################################
# Importing libraries                                                      #
############################################################################

import dendropy
from dendropy import TreeList
from dendropy import Tree
from pathlib import Path
import argparse

descr = """Prune tree based on a list of taxa to keep.
"""

parser = argparse.ArgumentParser(description=descr)

parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

required.add_argument("-t", "--tree", dest="tree_file", help="File with one or more trees in newick or nexus format", required=True)
required.add_argument("-r", "--retain", dest="taxa_to_retain", help="File with a list of taxa to retain, one per line.", required=True)
optional.add_argument("-o", "--output_name", dest="output_name", default="None", help="The base name for the output")
optional.add_argument("-f", "--tree_format", dest="tree_format", default="newick", help="The format of trees: newick or nexus. The output will be in the same format as the input. Default = newick")



#####################################
## Defining variables to be used   ##
#####################################


args = parser.parse_args()

tree_file = args.tree_file

taxa_to_retain = args.taxa_to_retain

output_name = args.output_name

tree_format = args.tree_format



########################### 
###### Script Begins ######
###########################

txt_file = open(taxa_to_retain, "r")

taxa_to_retain_list = txt_file.readlines()
for index, i in enumerate(taxa_to_retain_list):
    taxa_to_retain_list[index] = i.replace('\n', "")
    
for index, i in enumerate(taxa_to_retain_list):   
    taxa_to_retain_list[index] = i.replace('_', " ")

if tree_format == "nexus":
    extension = ".nex"
else:
    extension = ".nwck"
    
if output_name == "None":
    base_name = Path(tree_file).stem
    output_name = f"{base_name}_pruned{extension}"
    #alternative for other python versions
    #output_name = "{}_pruned.{}".format(base_name, tree_format)
else:
    output_name = f"{output_name}{extension}"


# Read tree
print("Reading Trees...")

trees=TreeList()
trees.read(path=tree_file,
    schema=tree_format)

print("Pruning Trees...")

for index, tree in enumerate(trees):
    percent=float((index/len(trees))*100)
    print(format(percent, '.2f'), "% done", sep = "", end = "\r" )
    try:
        tree.retain_taxa_with_labels(taxa_to_retain_list)
    except:
        trees[index]=Tree()

print("100% done")

pruned_trees=TreeList()

for index, tree in enumerate(trees):
    tree_string=tree.as_string(schema="newick")
    if tree_string != ';\n':
        pruned_trees.append(tree) 

pruned_trees_filtered=TreeList()

for index, tree in enumerate(pruned_trees):
    n_taxa = 0
    for leaf in tree.leaf_node_iter():
        n_taxa = n_taxa + 1
    if  n_taxa > 3 :
        pruned_trees_filtered.append(tree)


pruned_trees_filtered.write(
    path=output_name,
    schema=tree_format
    )

