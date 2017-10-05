# Revising assemblies with phylogenetic information

This repository contains the code for our manuscript on using phylogenetic information to revise and improve transcriptome assemblies. It is organized thusly:

* supplementary.Rmd is an R notebook complete with code to generate all the figures and analyses for the Supplementary Information.
* data/ is a folder with all the data that is used in supplementary.Rmd.
* code/ is a folder with all the code that was used to generate data/. This is organized into:
	* agalma/; bash scripts used to run Agalma1.0. This includes (1) a script to run Agalma1.0 from homologize up to the first genetree (i.e. before running treeinform), and (2) several scripts to run Agalma1.0 on the output from the first genetree run, starting from treeinform with different thresholds up to the 2nd genetree run.
	* phyldog/; bash scripts used to run phyldog on the 2nd genetree run's output from Agalma1.0.
	* functions.R; functions accessed by supplementary.Rmd.
	* with_internal.tre and species_chron.tre; non-ultrametric and ultrametric tree files used to generate the calibration nodes and times used in supplementary.Rmd. Species_chron.tre was generated with ape::chronos on with_internal.tre.
	* runjags.txt; mixture model structure and parameters specified in JAGS language, for use in runjags in supplementary.Rmd.