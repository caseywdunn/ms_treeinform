#!/bin/bash
#SBATCH --mem=16G
#SBATCH -t 24:00:00

source activate agalma

# NOTE: this script doesn't actually work, although the sqlite commands are correct.
# This is owing to the fact that I don't know how to non-interactively input the
# dot-commands, i.e. .mode, .separator, etc.

DATADIR=/users/aguang/scratch/treeinform/ms_treeinform/data/revisions/trinity
CORSETDIR=/users/aguang/data/aguang/corset
cd $CORSETDIR

# copy to agalma-corset.sqlite
"""
cp $DATADIR/agalma_trinity.sqlite $CORSETDIR/agalma_corset.sqlite

# transcriptome SRX IDs
srx=( SRX288276 SRX288285 SRX288430 SRX288431 SRX288432 )
#ids=( 20 21 23 24 22 )

sqlite3 agalma_corset.sqlite <<EOF

CREATE TABLE corset_clusters(
    trinity_gene TEXT NOT NULL,
    trinity_isoform TEXT NOT NULL,
    corset_cluster VARCHAR(255),
    run_id INT
);
EOF

#corset=( SRX288276-clusters SRX288285_cdhit-clusters SRX288430_cdhit-clusters SRX288431-clusters SRX288432_cdhit-clusters )
for i in "${srx[@]}"
do
    sed 's/TRINITY_//' $CORSETDIR/$i-clusters.txt > $CORSETDIR/$i-clusters.mod.txt
    sed -i 's/_i/\t/' $CORSETDIR/$i-clusters.mod.txt # add '.original' if on OSX
    sqlite3 agalma_corset.sqlite <<EOF

    CREATE TABLE tmp(
        trinity_gene TEXT NOT NULL,
        trinity_isoform TEXT NOT NULL,
        corset_cluster TEXT NOT NULL
        );
   .mode csv
   .separator \t
   .import $CORSETDIR/$i-clusters.mod.txt tmp
   ALTER TABLE tmp ADD COLUMN run_id int;
   UPDATE tmp SET run_id = (SELECT id from runs where catalog_id='$i' AND done=1);
   UPDATE tmp SET version = 0;
   UPDATE tmp SET model_id = (SELECT model_id from agalma_genes...
   UPDATE tmp SET model_id = (SELECT model_id from agalma_genes WHERE agalma_genes.gene = trinity_gene AND agalma_genes.isoform = trinity_isoform AND agalma_genes.run_id = run_id);
corset_cluster AS gene
   INSERT INTO corset_clusters SELECT trinity_gene,trinity_isoform,corset_cluster,run_id FROM tmp;
   DROP TABLE tmp;
EOF
done
"""


# need to alter these commands from tirnity_genes AS agalma_genes to just agalma_genes
# join corset_clusters as a new table in database with agalma_genes structure
sqlite3 agalma_corset.sqlite "CREATE TABLE corset_genes AS SELECT agalma_genes.run_id,agalma_genes.version,agalma_genes.model_id,corset_cluster AS gene,agalma_genes.isoform FROM trinity_genes AS agalma_genes JOIN corset_clusters ON agalma_genes.gene=corset_clusters.trinity_gene AND agalma_genes.isoform=corset_clusters.trinity_isoform AND agalma_genes.run_id=corset_clusters.run_id AND agalma_genes.version=0;"

#JGI_NEMVEC NCBI_HYDMAG
imports=( 5 6 )
for i in "${imports[@]}"
do
    sqlite3 agalma_corset.sqlite "INSERT INTO corset_genes SELECT * FROM trinity_genes WHERE run_id=$i AND version=0;"
done

# make corset_genes the agalma_genes table and make agalma_genes the trinity_genes table (for homologize)
#sqlite3 agalma_corset.sqlite "ALTER TABLE agalma_genes RENAME TO trinity_genes;"
sqlite3 agalma_corset.sqlite "ALTER TABLE corset_genes RENAME TO agalma_genes;"