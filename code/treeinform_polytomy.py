import sys, dendropy
import os
import itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ete3 import Tree
from collections import defaultdict

# usage: python ../ms_treeinform/code/treeinform_polytomy.py -t newick -b 0.05 -n 5301 -o poly.nw

# annotates each node in a tree with sum of branch lengths for subtree
def annotate_node(tree, threshold):
    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf():
            node.add_feature("branchlength", 0)
            node.add_feature("under", True)
        if not node.is_leaf():
            children = node.get_children()
            branchlength = children[0].get_distance(children[1]) + children[0].branchlength + children[1].branchlength
            node.add_feature("branchlength", branchlength)
            if branchlength < threshold: # adds flag for if node.bl under threshold
                node.add_feature("under", True)
            else: node.add_feature("under", False)

# sets candidate branch lengths to 0
# dictionary of fusing candidates for subtrees < threshold
def prune(tree):
    for node in tree.traverse(strategy='levelorder', is_leaf_fn=under):
        if node.branchlength != 0 and node.under == True:
            candidates = defaultdict(list)
            # get set of fusing candidates for subtree < threshold
            for leaf in node.get_leaves():
                species = leaf.name.split('@')[0]
                species_id = leaf.name.split('@')[1]
                candidates[species].append(leaf)
            # fuse candidates
            for k,v in candidates.items():
                if len(v) == 1: # single candidate, keep branch length
                    d = v[0].get_distance(node)
                    v[0].dist = d
                else: # multiple candidates, fuse
                    for i in v:
                        #i.delete()
                        node.remove_child(i)
                        node.add_child(i, dist=0)
    #tree.resolve_polytomy(recursive=True)
    # potential bug in original version: 2 clades with same species and node branchlengths ID
    return tree

# boolean for whether a subtree can be skipped
def under(node):
    return node.under

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Annotates each node with total length of subtree, then outputs a text file where each row is a set of sequences from the same species in a subtree below specified subtree length threshold")
    p.add_argument('-t', '--tree', help='tree file extension')
    p.add_argument('-b', '--threshold', help='branch length threshold')
    p.add_argument('-o', '--out', help='out extension')
    p.add_argument('-n', '--num_trees', help='number of gene trees')
    opts = p.parse_args()

    file_ext = opts.tree
    threshold = float(opts.threshold)
    out = opts.out
    num_trees = int(opts.num_trees)
    #treelist = dendropy.TreeList()
    #treelist.read_from_path(file_name, schema="newick")

    # overwrite file
    with open(out, 'w') as f:
        f.write('')

    for i in range(num_trees):
    #for i in range(1,2):
    #for tree in treelist:
        file_name = str(i) + "." + file_ext
        try:
            ete_tree = Tree(file_name)
            R = ete_tree.get_midpoint_outgroup()
            ete_tree.set_outgroup(R)
            annotate_node(ete_tree, threshold)
            tree=prune(ete_tree)
            out_file = str(i) + out
            tree.write(format=5, outfile=out_file)
        except:
            pass
        #ete_tree = Tree(tree.as_string("newick")[5:])
        
