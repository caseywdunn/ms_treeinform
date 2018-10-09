#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -N 8
#SBATCH -c 16
#SBATCH --mem=60G
#SBATCH -C intel
#SBATCH -J phyldog

#source activate agalma
# need to install phyldog...
module load phyldog/Aug2016

WORKDIR=/gpfs/scratch/aguang/ms_treeinform/data/drosophila_analysis/trinity/scratch/
CODEDIR=/gpfs/scratch/aguang/ms_treeinform/code

cd $WORKDIR
mkdir -p links
mkdir -p OptionFiles
mkdir -p ResultFiles

cd multalign-9
mkdir -p fa-gb # move fa-gb files for phyldog into here
cd alignments
for f in *.fa-gb
do
  NAME=${f%.fa-gb}
  OUT="$WORKDIR/links/$f"
  cat $f | python $CODEDIR/parse_links.py >$OUT
  mv $f ../fa-gb
done

cd ../..
cut -d: -f1 $WORKDIR/links/* | sort -u >SpeciesNames.txt

# inputFile.txt for prepareData.py
echo $WORKDIR/multalign-9/fa-gb > inputFile.txt
echo Protein >> inputFile.txt
echo fasta >> inputFile.txt
echo $WORKDIR/links >> inputFile.txt
echo $WORKDIR/OptionFiles >> inputFile.txt
echo $WORKDIR/ResultFiles >> inputFile.txt
echo yes >> inputFile.txt # Give starting species tree?
echo $CODEDIR/species_tree.tre >> inputFile.txt #absolute path to starting species tree
echo no >> inputFile.txt # Optimize species tree?
echo yes >> inputFile.txt # Optimize DL parameters?
echo branchwise >> inputFile.txt # How to optimize DL parameters?
echo no >> inputFile.txt # Want to assume all species have same number of genes?
echo yes >> inputFile.txt # Optimize gene trees?
echo 72 >> inputFile.txt # Time limit in hours

python $CODEDIR/prepareData.py < inputFile.txt

cd OptionFiles
for f in *.opt
do
  sed -i ':a;N;$!ba;s/input.sequence.sites_to_use=all/output.events.file=$(RESULT)$(DATA)_Events.txt \ninput.sequence.sites_to_use=all/' $f
done

srun phyldog param=$WORKDIR/OptionFiles/GeneralOptions.txt
