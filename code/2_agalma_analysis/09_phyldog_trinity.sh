#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH -n 12
#SBATCH -N 8
#SBATCH --mem=110G
#SBATCH -J phyldog
#SBATCH --array=1-9

source activate agalma
module load phyldog/Aug2016

treeinform=( 28 29 30 31 32 33 34 35 36 )
multalign=( 40 43 46 49 52 56 60 64 68 )

CODEDIR=/users/aguang/scratch/treeinform/ms_treeinform/code
WORKDIR=/users/aguang/scratch/treeinform/ms_treeinform/data/revisions/trinity/scratch/treeinform-${treeinform[$SLURM_ARRAY_TASK_ID-1]}
PHYLDOGDIR=/users/aguang/scratch/treeinform/ms_treeinform/data/revisions/trinity/phyldog/threshold-${SLURM_ARRAY_TASK_ID}

mkdir -p $PHYLDOGDIR

cd $PHYLDOGDIR
mkdir -p links
mkdir -p OptionFiles
mkdir -p ResultFiles

cd $WORKDIR/multalign-${multalign[$SLURM_ARRAY_TASK_ID-1]}
mkdir -p fa-gb # move fa-gb files for phyldog into here
cd alignments
for f in *.fa-gb
do
  NAME=${f%.fa-gb}
  OUT="$PHYLDOGDIR/links/$f"
  cat $f | python $CODEDIR/phyldog/parse_links.py >$OUT
  mv $f ../fa-gb
done
echo "links parsed"

cd ../..
cut -d: -f1 $PHYLDOGDIR/links/* | sort -u >$PHYLDOGDIR/SpeciesNames.txt
echo "Species names printed"

cd $PHYLDOGDIR
# inputFile.txt for prepareData.py
echo $WORKDIR/multalign-${multalign[$SLURM_ARRAY_TASK_ID-1]}/fa-gb > inputFile.txt
echo Protein >> inputFile.txt
echo fasta >> inputFile.txt
echo $PHYLDOGDIR/links >> inputFile.txt
echo $PHYLDOGDIR/OptionFiles >> inputFile.txt
echo $PHYLDOGDIR/ResultFiles >> inputFile.txt
echo yes >> inputFile.txt # Give starting species tree?
echo $CODEDIR/species_tree.tre >> inputFile.txt #absolute path to starting species tree
echo no >> inputFile.txt # Optimize species tree?
echo yes >> inputFile.txt # Optimize DL parameters?
echo branchwise >> inputFile.txt # How to optimize DL parameters?
echo no >> inputFile.txt # Want to assume all species have same number of genes?
echo yes >> inputFile.txt # Optimize gene trees?
echo 72 >> inputFile.txt # Time limit in hours

python $CODEDIR/phyldog/prepareData.py < inputFile.txt

cd OptionFiles
for f in *.opt
do
  sed -i ':a;N;$!ba;s/input.sequence.sites_to_use=all/output.events.file=$(RESULT)$(DATA)_Events.txt \ninput.sequence.sites_to_use=all/' $f
done

srun phyldog param=$PHYLDOGDIR/OptionFiles/GeneralOptions.txt