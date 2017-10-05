#!/bin/bash
#SBATCH -t 3-00:00:00
#SBATCH --mem=16G

module load agalma/1.0.0

# make a new database & new folders for each threshold
mkdir -p threshold-mix2
cd threshold-mix
cp ~/repos/ms_treeinform/data/agalma-7taxa-20170207.sqlite ./threshold-mix.sqlite

export AGALMA_DB=threshold-mix.sqlite

agalma treeinform -p 12 --id threshold-mix -t 0.02
agalma homologize --id threshold-mix
agalma multalign --id threshold-mix
agalma genetree --id threshold-mix