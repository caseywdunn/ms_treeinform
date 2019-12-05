#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -N 8
#SBATCH --mem=64G
#SBATCH -C intel
#SBATCH --array 5

export PATH=/gpfs/runtime/cbc_conda/bin/:$PATH
source activate_cbc_conda
conda activate agalma

DIR=/gpfs/data/cbc/aguang/treeinform/drosophila
DATADIR=$DIR/data
SCRATCHDIR=$DIR/scratch
srx=( SRX246999 SRX247001 SRX247003 SRX054483 SRX054470 SRX054462 SRX054487 )

#blastn -db $DATADIR/data/dmel_db -query $DATADIR/scratch/transcriptome-${SLURM_ARRAY_TASK_ID+1}/trinity_out_dir/Trinity.fasta -outfmt 6 -out $DATADIR/scratch/blast/${srx[$SLURM_ARRAY_TASK_ID]}.tsv -num_threads 8 -max_target_seqs 20 -evalue 1e-10
blastn -db $DATADIR/${srx[$SLURM_ARRAY_TASK_ID]} -query $DIR/scratch/transcriptome-${SLURM_ARRAY_TASK_ID+1}/trinity_out_dir/Trinity.fasta -outfmt 6 -out $SCRATCHDIR/blast/${srx[$SLURM_ARRAY_TASK_ID]}.tsv -num_threads 8 -max_target_seqs 20 -evalue 1e-10
