#!/bin/bash
#SBATCH -t 3-00:00:00
#SBATCH --mem=24G

conda create -n corset
conda config --add channels conda-forge
conda config --add channels bioconda
conda install -n corset corset
conda install -n corset bowtie2
conda install -n corset samtools

source activate corset