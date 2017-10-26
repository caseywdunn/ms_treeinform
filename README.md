# Revising assemblies with phylogenetic information

This repository contains the code for our manuscript on using phylogenetic information to revise and improve transcriptome assemblies in [Agalma1.0](http://bitbucket.org/caseywdunn/agalma). It is described in this preprint:

Guang A, Howison M, Zapata F, Lawrence CE, Dunn CW. Revising transcriptome assemblies with phylogenetic information in Agalma 1.0. bioRxiv 202416; doi: [https://doi.org/10.1101/202416](https://doi.org/10.1101/202416)

# Organization

The repository is organized thusly:

* **treeinform_preprint.pdf:** pdf of main manuscript preprint.
* **supplementary.Rmd:** R notebook complete with code to generate all the figures and analyses for the Supplementary Information.
* **figure1.Rmd:** R notebook with code to generate Figure 1 in the main manuscript.
* **data/:** folder with all the data that is used in supplementary.Rmd and for Figure 1.
* **code/:** folder with all the code that was used to generate **data/**. This is organized into:
	* **agalma/:** bash scripts used to run Agalma1.0 on output from a previous genetree run from Agalma1.0, starting from treeinform with different thresholds up to the 2nd genetree run. The sqlite database storing the previous genetree run is not in this repository as it is 1.5GB; if you would like to access the database please email us.
	* **phyldog/:** bash scripts used to run phyldog on the 2nd genetree run's output from Agalma1.0.
	* **functions.R:** functions accessed by supplementary.Rmd.
	* **with_internal.tre** and **species_chron.tre:** non-ultrametric and ultrametric tree files used to generate the calibration nodes and times used in supplementary.Rmd. Species_chron.tre was generated with ape::chronos on with_internal.tre.
	* **runjags.txt:** mixture model structure and parameters specified in JAGS language, for use in runjags in supplementary.Rmd.

# Bash scripts for Agalma1.0 and phyldog

All bash scripts were run on the epscor condo in the Center for Computing and Visualization at Brown inside a conda env. They are written to work with SLURM. If you would like help with running the scripts yourselves or to access the sqlite database the scripts were run on, please email us.