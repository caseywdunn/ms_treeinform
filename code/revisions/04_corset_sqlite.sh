#!/bin/bash
#SBATCH -t 1:00:00

source activate agalma

DATADIR=/users/aguang/scratch/treeinform/ms_treeinform/data/revisions/trinity
CORSETDIR=/users/aguang/data/aguang/corset
cd $CORSETDIR

# copy to agalma-corset.sqlite
#cp $DATADIR/agalma_trinity.sqlite $CORSETDIR/agalma_corset.sqlite

# transcriptome SRX IDs
srx=( SRX288276 SRX288285 SRX288430 SRX288431 SRX288432 )
ids=( 20 21 23 24 22 )

# create a corset clusters table in database
#sqlite3 agalma_corset.sqlite "CREATE TABLE corset_clusters(
#    trinity_gene TEXT NOT NULL,
#    trinity_isoform TEXT NOT NULL,
#    corset_cluster VARCHAR(255),
#    run_id INT
#);"

for i in "${srx[@]}"
do
#    sed 's/TRINITY_//' $CORSETDIR/$i-clusters.txt > $CORSETDIR/$i-clusters.mod.txt
#    sed -i 's/_i/\t/' $CORSETDIR/$i-cluster.mod.txt # add '.original' if on OSX
    sqlite3 agalma_corset.sqlite "CREATE TABLE tmp(
        trinity_gene TEXT NOT NULL,
        trinity_isoform TEXT NOT NULL,
        corset_cluster TEXT NOT NULL
        );"
   sqlite3 agalma_corset.sqlite ".mode csv \n .separator \t \n .import $CORSETDIR/$i-clusters.mod.txt tmp"
   sqlite3 agalma_corset.sqlite "ALTER TABLE tmp ADD COLUMN run_id int;"
   sqlite3 agalma_corset.sqlite "UPDATE tmp SET run_id = (SELECT id from runs where catalog_id='$i' AND done=1);"
   sqlite3 agalma_corset.sqlite "INSERT INTO corset_clusters SELECT trinity_gene,trinity_isoform,corset_cluster,run_id FROM tmp;"
   sqlite3 agalma_corset.sqlite "DROP TABLE tmp;"
done

# join corset_clusters as a new table in database with agalma_genes structure
sqlite3 agalma_corset.sqlite "CREATE TABLE corset_genes AS SELECT agalma_genes.run_id,agalma_genes.version,agalma_genes.model_id,corset_cluster AS gene,agalma_genes.isoform FROM agalma_genes JOIN corset_clusters ON agalma_genes.gene=corset_clusters.trinity_gene AND agalma_genes.isoform=corset_clusters.trinity_isoform AND agalma_genes.run_id=corset_clusters.run_id AND agalma_genes.version=0;"

# make corset_genes the agalma_genes table and make agalma_genes the trinity_genes table (for homologize)
sqlite3 agalma_corset.sqlite "ALTER TABLE agalma_genes RENAME TO trinity_genes;"
sqlite3 agalma_corset.sqlite "ALTER TABLE corset_genes RENAME TO agalma_genes;"