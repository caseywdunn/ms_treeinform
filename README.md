# Revising assemblies with phylogenetic information

This repository contains the code for our work on using phylogenetic information to revise and improve transcriptome assemblies in [Agalma1.0](http://bitbucket.org/caseywdunn/agalma). It is described in this manuscript:

Guang A, Howison M, Zapata F, Lawrence C, Dunn CW (2021) Revising transcriptome assemblies with phylogenetic information. PLOS ONE 16(1): e0244202. https://doi.org/10.1371/journal.pone.0244202

# Organization

The repository is organized thusly:

* **treeinform_preprint.pdf:** pdf of main manuscript preprint.
* **treeinform_PLoS_ONE.pdf:** pdf of manuscript submission to PLOS ONE.
* **plos_one.Rmd:** R notebook to generate figures for PLOS ONE submission.
* **bioinformatics_submission/:** folder with R notebooks and code for a prior submission to Bioinformatics
* **data/:** folder with all the data that is used in plos_one.Rmd.
* **code/:** folder with all the code that was used to generate **data/**. This is organized into:
	* **1_drosophila_analysis/:** bash scripts used to run Agalma1.0 on Drosophila reads from SRA (accession numbers in the scripts) with different treeinform thresholds.
	* **2_agalma_analysis/:** bash scripts used to run Agalma1.0 on Siphonophora reads from SRA (accession numbers in the scripts), JGI, NCBI-EST with different treeinform thresholds.
	* **cafe/:** bash and python scripts used to run CAFE for initial estimates of lambda and mu in our mixture model.
	* **phyldog/:** bash scripts used to run phyldog on the 2nd genetree run's output from Agalma1.0.
	* **functions.R:** functions accessed by plos_one.Rmd
	* **summary_stats.R:** additional functions accessed by plos_one.Rmd
	* **with_internal.tre** and **species_chron.tre:** non-ultrametric and ultrametric tree files used to generate the calibration nodes and times used in plos_one.Rmd. Species_chron.tre was generated with ape::chronos on with_internal.tre.
	* **runjags.txt:** mixture model structure and parameters specified in JAGS language, for use in runjags in plos_one.Rmd.

# Bash scripts for Agalma1.0 and phyldog

All bash scripts were run on the epscor and datasci condos in the Center for Computing and Visualization at Brown inside a conda env. They are written to work with SLURM. If you would like help with running the scripts yourselves or to access the sqlite database the scripts were run on, please email us.
