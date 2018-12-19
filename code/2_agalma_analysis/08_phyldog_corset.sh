#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH -n 16
#SBATCH -N 8
#SBATCH --mem=60G
#SBATCH -J phyldog
#SBATCH -C intel

source activate agalma
module load phyldog/Aug2016

CODEDIR=/users/aguang/scratch/treeinform/ms_treeinform/code
WORKDIR=/users/aguang/scratch/treeinform/ms_treeinform/data/revisions/corset/scratch
PHYLDOGDIR=/users/aguang/scratch/treeinform/ms_treeinform/data/revisions/corset/phyldog

mkdir -p $PHYLDOGDIR

cd $PHYLDOGDIR
mkdir -p links
mkdir -p OptionFiles
mkdir -p ResultFiles


cd $WORKDIR/genetree-101/alignments
for f in *.fa
do
  NAME=${f%.fa}
  OUT="$PHYLDOGDIR/links/$f"
  cat $f | python $CODEDIR/phyldog/parse_links.py >$OUT
done
echo "links parsed"

cd ../..
cut -d: -f1 $PHYLDOGDIR/links/* | sort -u >$PHYLDOGDIR/SpeciesNames.txt
echo "Species names printed"

cd $PHYLDOGDIR
# inputFile.txt for prepareData.py
echo $WORKDIR/genetree-101/alignments > inputFile.txt
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
sed -i 's/branchProbabilities.optimization=average_then_branchwise/branch.expected.numbers.optimization=branchwise/' GeneralOptions.txt


GENETREEDIR=$WORKDIR/genetree-101/trees/
for f in *.opt
do
  NEWICK=${f/.opt/.newick}
  sed -i ':a;N;$!ba;s/input.sequence.sites_to_use=all/output.events.file=$(RESULT)$(DATA)_Events.txt \ninput.sequence.sites_to_use=all/' $f
  sed -i 's/######## Second, model options ########/rearrangement.gene.tree=spr \n######## Second, model options ########/' $
  sed -i "s+init.gene.tree=bionj+init.gene.tree=user\ngene.tree.file=${GENETREEDIR}${NEWICK}+" $f
  sed -i 's/.fa.reduced/.fa/' $f
done

srun phyldog param=$PHYLDOGDIR/OptionFiles/GeneralOptions.txt