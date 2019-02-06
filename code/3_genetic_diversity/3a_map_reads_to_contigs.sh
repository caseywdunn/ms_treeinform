#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH --mem=8G
#// SBATCH --mem=64G
#// SBATCH -n 8
#SBATCH --array=7-11

# Assumes you have reads and contigs from running code in drosophila_assembly
# and agalma_assembly already or through an alternative method.

source activate agalma

DROSOPHILA_DIR=/gpfs/data/datasci/aguang/treeinform/drosophila
AGALMA_DIR=/gpfs/data/datasci/aguang/treeinform/agalma/trinity
fastqs=( SRR768436 SRR768438 SRR768440 SRR166831 SRR166818 SRR166810 SRR166835  SRR871525 SRR871526 SRR871527 SRR871528 SRR871529 )
ids=( 1 2 3 4 5 6 7 5 6 7 8 9 ) # transcriptome IDs; 1-7 is drosopohila and 5-9 is agalma
srx=( SRX246999 SRX247001 SRX247003 SRX054483 SRX054470 SRX054462 SRX054487 SRX288276 SRX288285 SRX288430 SRX288431 SRX288432 )

if [ $SLURM_ARRAY_TASK_ID -lt 7 ]
then
    mkdir -p $DROSOPHILA_DIR/scratch/genetic_diversity
    cd $DROSOPHILA_DIR/scratch/genetic_diversity

    ID=${ids[$SLURM_ARRAY_TASK_ID]}
    fastq=${fastqs[$SLURM_ARRAY_TASK_ID]}
    srx_id=${srx[$SLURM_ARRAY_TASK_ID]}

    # map reads to contigs with BWA
    R1=$DROSOPHILA_DIR/data/${fastq}_1.fastq
    R2=$DROSOPHILA_DIR/data/${fastq}_2.fastq
    bwa index $DROSOPHILA_DIR/scratch/transcriptome-${ID}/trinity_out_dir/Trinity.fasta
    bwa mem -t 8 $DROSOPHILA_DIR/scratch/transcriptome-${ID}/trinity_out_dir/Trinity.fasta $R1 $R2 > ${srx_id}.sam

    samtools view -bS ${srx_id}.sam > ${srx_id}.bam
    rm ${srx_id}.sam # you're free to keep the sam files if you want but they are very large
    samtools sort -m 8G -@ 8 -o ${srx_id}.sorted.bam -T temp -O bam ${srx_id}.bam
    samtools depth ${srx_id}.sorted.bam > ${srx_id}_coverage.txt
else
#    mkdir -p $AGALMA_DIR/scratch/genetic_diversity
    cd $AGALMA_DIR/scratch/genetic_diversity

    ID=${ids[$SLURM_ARRAY_TASK_ID]}
    fastq=${fastqs[$SLURM_ARRAY_TASK_ID]}
    srx_id=${srx[$SLURM_ARRAY_TASK_ID]}

    # map reads to contigs with BWA
    R1=$AGALMA_DIR/data/${fastq}_1.fastq
    R2=$AGALMA_DIR/data/${fastq}_2.fastq
#    bwa index $AGALMA_DIR/scratch/transcriptome-${ID}/trinity_out_dir/Trinity.fasta
#    bwa mem -t 8 $AGALMA_DIR/scratch/transcriptome-${ID}/trinity_out_dir/Trinity.fasta $R1 $R2 > ${srx_id}.sam

#    samtools view -bS ${srx_id}.sam > ${srx_id}.bam
#    rm ${srx_id}.sam # you're free to keep the sam files if you want but they are very large
#    samtools sort -m 8G -@ 8 -o ${srx_id}.sorted.bam -T temp -O bam ${srx_id}.bam
    samtools depth ${srx_id}.sorted.bam > ${srx_id}_coverage.txt
fi
