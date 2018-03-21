#!/bin/bash

# this script replicates setting up a conda environment of agalma1.1.0 with trinity/2.5.1 loaded,
# on Brown University's CCV computing cluster Oscar. At this time (February 19, 2018) this is the
# latest version of trinity.

conda create -n agalma -c dunnlab agalma