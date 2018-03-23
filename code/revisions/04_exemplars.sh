#!/bin/bash
#SBATCH -t 3-00:00:00
#SBATCH --mem=24G

# run RSEMeval and select exemplars based on max(confidence)

source activate agalma

DATADIR=/gpfs/scratch/aguang/treeinform/revisions/trinity
CORSETDIR=/users/aguang/scratch/treeinform/ms_treeinform/data/revisions/corset
cd $CORSETDIR

# transcriptome SRX IDs
fastq=( SRX288276 SRX288285 SRX288430 SRX288431 SRX288432 )

# insert size (manually entered)

#RSEMeval
rsem-eval-calculate-score --num-threads 8 --paired-end --seed 807633215 $DATADIR/scratch/transcriptome-$[$SLURM_ARRAY_TASK_ID-1]/filtered.0.fq $DATADIR/scratch/transcriptome-$SLURM/filtered.1.fq $DATADIR/scratch/transcriptome-$SLURM/${fastq[$i]}.fa rsem_eval

--bam input assembly_fasta_file sample_name L
/gpfs_home/aguang/miniconda3/envs/agalma/opt/rsem-eval-1.9/bin/rsem-eval-calculate-score --num-threads 8 --paired-end --seed 807633215 /gpfs/scratch/aguang/treeinform/ms_treeinform/data/revisions/trinity/scratch/transcriptome-1/filtered.0.fq /gpfs/scratch/aguang/treeinform/ms_treeinform/data/revisions/trinity/scratch/transcriptome-1/filtered.1.fq /gpfs/scratch/aguang/treeinform/ms_treeinform/data/revisions/trinity/scratch/transcriptome-1/SRX288276.fa rsem_eval 288.988 2>>/gpfs/scratch/aguang/treeinform/ms_treeinform/data/revisions/trinity/scratch/transcriptome-1/rsem-eval-calculate-score.log 1>>/gpfs/scratch/aguang/treeinform/ms_treeinform/data/revisions/trinity/scratch/transcriptome-1/rsem-eval-calculate-score.log

#bowtie
for i in {1..5}
do
    bowtie2-build $DATADIR/scratch/transcriptome-${i}/trinity_out_dir/Trinity.fasta $CORSETDIR/transcriptome-${i}
    FILES=`ls $DATADIR/data/${fastq[$i]}_1.fq | sed 's/_1.fq//g'` # need to work on this, also it is fastq data
    for F in $FILES ; do
        R1=${F}_1.fq
        R2=${F}_2.fq
        bowtie2 --all -S $CORSETDIR/transcriptome-${i} -1 $R1 -2 $R2 > ${F}.sam  
        samtools view -S -b ${F}.sam > ${F}.bam
    done
done

#corset
corset $CORSETDIR/*.bam