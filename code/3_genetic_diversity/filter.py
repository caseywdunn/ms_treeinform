#!/usr/bin/env python
from ete3 import Tree
from collections import defaultdict

# usage:
# python filter.py -o test.tab -i test_coverage.txt

# return number of genes for each species in a gene tree
def filter(infile, outfile):
    gene_family = defaultdict(int)
    for leaf in tree.iter_leaves():
        species = leaf.name.split('@')[0]
        gene_family[species] = gene_family[species]+1
    return gene_family

# write tab-delimited file for CAFE
def write(postfix, directory, nmin, nmax, taxon, out):
    # process taxon file
    t = ''
    with open(taxon, 'r') as f:
        t = f.read().split(',')

    # process trees
    sizes = []
    for i in range(nmin,nmax+1):
        tree_str = directory + str(i) + "." + postfix
        # error catching
        try:
            print(tree_str)
            tree = Tree(tree_str)
        except:
            print("couldn't read")
            continue
        gene_family = num_genes(tree, t)
        sizes.append(gene_family)

    # write file
    # does not check for if taxon names in file match those from gene trees
    with open(out, 'w') as f:
        f.write('Description\tID')
        for species in t:
            f.write('\t%s' % (species))
        f.write('\n')
        for i in range(len(sizes)):
            f.write('%d\t%d' % (i, i))
            gene_family = sizes[i]
            for species in t:
                f.write('\t%d' % (gene_family[species]))
            f.write('\n')

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='build gene family size input for CAFE from set of gene trees')

    ds = ' [%(default)s]'
    parser.add_argument('-i', '--infile', help='infile')
    parser.add_argument('-o', '--out', help='output file')
    opts = parser.parse_args()

    infile = opts.infile
    out = opts.out

    write(postfix, directory, nmin, nmax, taxon, out)
