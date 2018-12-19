#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH --mem=64G
#SBATCH -n 8
#SBATCH --array=1-13

# Assumes you have reads and contigs from running code in drosophila_assembly
# and agalma_assembly already or through an alternative method.

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
    mkdir -p $DROSOPHILA_DIR/scratch/genetic_diversity
    cd $DROSOPHILA_DIR/scratch/genetic_diversity

    ID=${drosophila_ids[$SLURM_ARRAY_TASK_ID]}
    fastq=${drosophila_fastqs[$SLURM_ARRAY_TASK_ID]}
    srx_id=${drosophila_srx[$SLURM_ARRAY_TASK_ID]}

    # map reads to contigs with BWA
    R1=$DROSOPHILA_DIR/data/${fastq}_1.fastq
    R2=$DROSOPHILA_DIR/data/${fastq}_2.fastq
    bwa index $DROSOPHILA_DIR/scratch/transcriptome-${ID}/trinity_out_dir/Trinity.fasta
    bwa mem -t 8 $DROSOPHILA_DIR/scratch/transcriptome-${ID}/trinity_out_dir/Trinity.fasta $R1 $R2 > ${srx_id}.sam

    samtools view -bS ${srx_id}.sam > ${srx_id}.bam
    rm ${srx_id}.sam # you're free to keep the sam files if you want but they are very large
    samtools sort -m 8G -@ 8 -o ${srx_id}.sorted.bam -T temp -O bam ${srx_id}.bam
    samtools depth ${srx_id}.bam > ${srx_id}_coverage.txt
else
    mkdir -p $AGALMA_DIR/scratch/genetic_diversity
    cd $AGALMA_DIR/scratch/genetic_diversity

    ID=${agalma_ids[$SLURM_ARRAY_TASK_ID]}
    fastq=${agalma_fastqs[$SLURM_ARRAY_TASK_ID]}
    srx_id=${agalma_srx[$SLURM_ARRAY_TASK_ID]}

    # map reads to contigs with BWA
    R1=$AGALMA_DIR/data/${fastq}_1.fastq
    R2=$AGALMA_DIR/data/${fastq}_2.fastq
    bwa index $AGALMA_DIR/scratch/transcriptome-${ID}/trinity_out_dir/Trinity.fasta
    bwa mem -t 8 $AGALMA_DIR/scratch/transcriptome-${ID}/trinity_out_dir/Trinity.fasta $R1 $R2 > ${srx_id}.sam

    samtools view -bS ${srx_id}.sam > ${srx_id}.bam
    rm ${srx_id}.sam # you're free to keep the sam files if you want but they are very large
    samtools sort -m 8G -@ 8 -o ${srx_id}.sorted.bam -T temp -O bam ${srx_id}.bam
    samtools depth ${srx_id}.bam > ${srx_id}_coverage.txt
fi
