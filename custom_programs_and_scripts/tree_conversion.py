#! /usr/bin/env python
# -*- coding: utf-8 -*-


############################################################################
## Importing libraries                                                     #
############################################################################

import dendropy
from dendropy import TreeList
from dendropy import Tree
from pathlib import Path
import argparse

descr = """Convert tree formats.
"""

parser = argparse.ArgumentParser(description=descr)

parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

required.add_argument("-t", "--tree", dest="tree_file", help="File with one or more trees in newick or nexus format", required=True)
required.add_argument("-of", "--output_format", dest="output_format", default="None", help="The output format: newick or nexus.")
required.add_argument("-if", "--input_format", dest="input_format", default="newick", help="The imput format: newick or nexus.")
optional.add_argument("--keep_annotation", default=False, action="store_true", help="Wheter to keep branch annotation")


#####################################
## Defining variables to be used   ##
#####################################


args = parser.parse_args()

tree_file = args.tree_file

output_format = args.output_format

input_format = args.input_format

########################### 
###### Script Begins ######
###########################


if output_format == "nexus":
    extension = ".nex"
else:
    extension = ".nwck"

base_name = Path(tree_file).stem
output_name = f"{base_name}{extension}"
#alternative for other python versions
#output_name = "{}_pruned.{}".format(base_name, extension)


# Read tree
print("Reading Trees...")

trees=TreeList()
trees.read(path=tree_file,
    schema=input_format)

print("Converting Trees...")

if args.keep_annotation:
    trees.write(path=output_name,
                schema=output_format
    )
else:
    for tree in trees:
        for node in tree.preorder_node_iter():
            node.label=None
    
    trees.write(path=output_name,
                schema=output_format
    )