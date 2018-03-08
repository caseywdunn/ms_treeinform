#!/bin/bash
#SBATCH -t 3-00:00:00
#SBATCH --mem=24G

module load conda
conda create -n corset
conda config --add bioconda
conda install -n corset samtools
conda install -n corset corset
conda install -n corset bowtie2

source activate corset