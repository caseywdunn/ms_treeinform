#!/bin/bash

# this script replicates setting up a conda environment of agalma1.0.1, with trinity2.5.1 loaded,
# on Brown University's CCV computing cluster Oscar. At this time (February 19, 2018) this is the
# latest version of trinity.

module load conda # on your own computer you will need to have conda already installed. you can comment 
# out this line if that is the case

conda create -n agalma

conda install -c conda-forge -n agalma pycairo
conda install -c conda-forge -n agalma libgcc-ng
conda install -c conda-forge -n agalma libstdcxx-ng
conda install -c dunnlab -n agalma agalma=1.0.1
conda remove -n agalma agalma trinity
conda install -c dunnlab -n agalma trinity=2.6.5
conda install -f -c dunnlab -n agalma agalma=1.0.1

source activate agalma