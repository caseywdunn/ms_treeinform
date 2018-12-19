#!/bin/bash
#SBATCH -t 4:00:00
#SBATCH --mem=24G
#SBATCH -n 8
#SBATCH --array=0-13

source activate agalma

DROSOPHILA_DIR=/gpfs/data/datasci/aguang/treeinform/drosophila
drosophila_fastqs=( SRR768436 SRR768438 SRR768440 SRR166831 SRR166818 SRR166810 SRR166835 )
drosophila_ids=( 1 2 3 4 5 6 7 ) # transcriptome IDs
drosophila_srx=( SRX246999 SRX247001 SRX247003 SRX054483 SRX054470 SRX054462 SRX054487 )

AGALMA_DIR=/gpfs/data/datasci/aguang/treeinform/agalma
agalma_ids=( 1 2 3 4 5 6 7 )
agalma_srx=( SRX288276 SRX288285 SRX288430 SRX288431 SRX288432 )
agalma_fastqs=( SRR871525 SRR871526 SRR871527 SRR871528 SRR871529 )

if [ $SLURM_ARRAY_TASK_ID -lt 7 ]
then
    cd $DROSOPHILA_DIR/scratch/genetic_diversity

    ID=${drosophila_ids[$SLURM_ARRAY_TASK_ID]}
    fastq=${drosophila_fastqs[$SLURM_ARRAY_TASK_ID]}
    srx_id=${drosophila_srx[$SLURM_ARRAY_TASK_ID]}

    echo "${srx_id}.sorted.bam ID" > list.txt
    reads2snps -bam_list_file list.txt -bam_ref_file $DROSOPHILA_DIR/scratch/transcriptome-${ID}/trinity_out_dir/Trinity.fasta -out ${srx_id} -nbth 8
else
    cd $AGALMA_DIR/scratch/genetic_diversity

    ID=${agalma_ids[$SLURM_ARRAY_TASK_ID]}
    fastq=${agalma_fastqs[$SLURM_ARRAY_TASK_ID]}
    srx_id=${agalma_srx[$SLURM_ARRAY_TASK_ID]}

    echo "${srx_id}.sorted.bam ID" > list.txt
    reads2snps -bam_list_file list.txt -bam_ref_file $AGALMA_DIR/scratch/transcriptome-${ID}/trinity_out_dir/Trinity.fasta -out ${srx_id} -nbth 8
fi

