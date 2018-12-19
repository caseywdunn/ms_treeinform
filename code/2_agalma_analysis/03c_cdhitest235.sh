#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH -n 8
#SBATCH --mem=60G
#SBATCH --array=2,3,5

source activate corset

DATADIR=/gpfs/scratch/aguang/treeinform/ms_treeinform/data/revisions/trinity
mkdir -p /users/aguang/scratch/treeinform/ms_treeinform/data/revisions/corset
CORSETDIR=/users/aguang/data/aguang/corset
cd $CORSETDIR

# transcriptome SRX IDs
srx=( SRX288276 SRX288285 SRX288430 SRX288431 SRX288432 )
fastqs=( SRR871525 SRR871526 SRR871527 SRR871528 SRR871529 )

# transcriptome ids
ids=( 20 21 23 24 22 )

#SLURM_ARRAY_TASK_ID=2
ID=${ids[$SLURM_ARRAY_TASK_ID-1]}
fastq=${fastqs[$SLURM_ARRAY_TASK_ID-1]}
srx_id=${srx[$SLURM_ARRAY_TASK_ID-1]}


#cd-hit-est -i $DATADIR/scratch/transcriptome-$ID/trinity_out_dir/Trinity.fasta -o cdhit_Trinity.fasta  -c 1 -M 60000 -T 8
#bowtie2-build cdhit_Trinity.fasta $CORSETDIR/cdhit_transcriptome-$ID
#R1=$DATADIR/data/${fastq}_1.fastq
#R2=$DATADIR/data/${fastq}_2.fastq
#bowtie2 -p 16 --all -x $CORSETDIR/cdhit_transcriptome-$ID -1 $R1 -2 $R2 > ${srx_id}_cdhit.sam  
#samtools view -S -b ${srx_id}_cdhit.sam > ${srx_id}_cdhit.bam
#rm ${srx_id}_cdhit.sam # if we don't do this we will be way over quota
#corset ${srx_id}_cdhit.bam -r true-stop
corset -i corset ${srx_id}_cdhit.bam.corset-reads -p ${srx_id}_cdhit -f true -x 10 -l 15