Reads2snp v2.0 : a program for calling genotypes and SNPs based on NGS transcriptome reads.
========


Introduction
--------

reads2snp is a SNP and genotype caller: it predicts the genotype of distinct individuals at distinct positions of a set of sequences based on read mapping / read counts. Its typical input is a bam file. Its typical output is a vcf file. It is written in C++, based on the bio++ libraries, multi-threaded with openMP, available under Linux and MacOS.

reads2snp is developed in the context of the PopPhyl project (http://kimura.univ-montp2.fr/PopPhyl/), in which transcriptome data are used to characterize the population genomics of non-model organisms.

The specificities of reads2snp are to (i) make no prior hypothesis regarding the population frequency of alleles, (ii) estimate the sequencing error rate from read counts only, (iii) account for unequal expression bias between alleles, (iv) take specific care of spurious SNPs due to hidden paralogy.


Installation

The distributed binaries for GNU/Linux are statically linked and ready to run. The MacOSX bundle is also ready to run on any system but is dynamically linked and provides the shared libraries.


Compilation
--------

If, for any reason, a source compilation is needed, here are the dependencies :

*   gcc >= 4.4 with openMP support
*   Bio++ >= 2.0 (bpp-core, bpp-phyl, bpp-seq)
*   libgsl >= 1.0
*   zlib >= 1.1

The makefile should be adapted to fit with Bio++ path if needed.

The “make all” command produces a statically linked binary under GNU/Linux. The “make mac” command produces a dynamically linked binary under MacOSX and the script  utils/build_mac_bundle.sh  produces an autonomous binary with the shared libraries in the root source folder.


Execution
--------

reads2snp [options] (-alr alr_file | -bam bam_file -bamref bam_ref_file)


Input
--------

reads2snp needs information on the mapping of reads from distinct individuals to a common set of reference sequences (contigs). This is typically passed through bam files:

reads2snp -bam_list_file my_list.txt -bam_ref_file my_ref.fasta

The bam list file must be a two-colon file including one line per analyzed individual, such as:

file1.bam   indiv1_name
file2.bam   indiv2_name
…

reads2snp will analyze the listed bam files, and call individuals according to the listed names. Each bam file is supposed to contain the results of the mapping of reads from the corresponding individual to a reference that is shared across individuals. The reference must be passed to reads2snp too using the -bam_ref_file option. The reference file is a fasta file containing one or several sequences (contigs).

bam files can be produced by your favorite mapper, or generated from sam files using samtools :

samtools view -@ 4 -uhSt exple.fasta exple.sam | samtools sort -@ 4 - exple

This will produce a sorted bam file “exple.bam” using 4 threads.

Alternatively, reads2snp can take as input an alr file, which only contains read counts, not read alignments – and is thus much more compact. A description of the alr file format is provided at the bottom of this file. An alr input file is passed using the -alr option (no reference file needed):

reads2snp -alr exple.alr


Analysis
--------

reads2snp analyzes contigs sequentially. For each contig, the sequencing error rate and, if required, the allelic expression bias are estimated in the maximum-likelihood framework using a multinomial model for read counts. Then, for each position and each individual, the posterior probabilities of the 10 possible diploid genotypes are calculated, and the most likely genotype is called if its posterior probability is above some user-specificied threshold (Tsagkogeorga et al 2012). Only sufficiently covered positions are genotyped.

A second pass on the data is subsequently performed, in search for dubious SNPs potentially due to hidden paralogy. This second step (paraclean) only concerns positions at which at least one heterozygous genotype was called. A likelihood-ratio test is applied that compares a one-locus and a two-locus model. If the former is rejected in favor of the latter, the considered SNP is filtered out. A Dirichlet-multinomial model is used in this analysis to account for the overdispersion of read counts (Gayral et al 2013).


Output
--------

reads2snp writes as output the predicted diploid genotype for each individual at each position of each contig. The two main output files are:

*   a variant calling (.vcf) file containing the called genotypes and associated information regarding the filtering process; uncalled/filtered genotypes are coded as ".|.";

*   a fasta file (.fas) containing essentially the same information in the form of aligned DNA sequences – two per contig per diploid individual. Uncalled/filtered genotypes are coded as "NN"; beware: no phasing is done by reads2snp; SNPs in the .fas output files are arbitrarily phased;

reads2snp writes two additional files (for expert users only):
*   a genotype (.gen) file containing the posterior probability of maximally likely genotypes;
*   a parameter (.par) file containing error rate and allelic expression bias estimates for each contig;

All four output files share the same basename, which can be set using the -out option:

reads2snp -alr exple.alr -out my_outfile


Options
--------

[I/O options]

-bam <string> -bamref <string>  : bam and reference input file names

-alr <string>   : alr input file name

-out <string>   : output file basename

[main options]

-min <int>  : sets the minimal number of reads required to call a genotype (default = 10; increase to get less numerous, safer calls)

-th1 <float>    : sets the minimal posterior probability required to call a genotype (default = 0.95; increase to get less numerous, safer calls)

-par <0|1>  : options for hidden paralogy filtering [0: no filtering; default=1: do filtering]

-th2 <float>    : sets the maximal p-value required to reject a paralogous SNP (default = 0.01; increase to get more stringent filtering)

[advanced options]

-nbth <int> : number of openMP threads to be used during the execution (for multi-cpu computers only; default = 0; increase to speed up the exectution)

-aeb        : account for allelic expression bias (model 2 in Tsagkogeorga et al 2012)

-rand       : genotypes will be randomly drawn according to their posterior probabilities, so that not always the most probable genotype is called, thus avoiding the bias (threshold effect) identified by Kim et al 2011 BMC Bioinform 12: 231 (see Gayral et al 2013).

-acc|-tol|-for  : options regarding how to deal with instances of complex allele counts; a complex allele count occurs when, at a given position in a given individual, the third allele is seen in at least 3 reads and at least 5% of the reads. -acc: a genotype will be called anyway; -tol: missing data will be called in case of complex allele counts [default]; -for: complex allele counts forbidden, execution halted.

-fis <float>    : set the Fis (expected heterozygote depletion); default = 0 (Hardy-Weinberg equilibrium); increase it to get more homozygous genotypes; the same Fis is assumed for genotype calling and paralog filtering; useful for selfing species (e.g. see Nabholz et al 2014 Mol Ecol in press).

-opt <bfgs|newton>  : likelihood optimizer used in paralogue filtering; default = bfgs;

-pre <float>    : precision on parameters in likelihood optimization; default = 0.001

-spa        : remove contigs containing at least one paralogous SNPs

[bam files options]

-rlg <int>  : only consider reads longer than threshold

-bqt <int>  : only consider positions of quality above threshold

-rqt <int>  : only consider reads of mapping quality above threshold


Developed by
--------

Vincent Cahais, Julien Veyssier, Khalid Belkhir, Sylvain Glémin, Nicolas Galtier


Supported by
--------

[Montpellier Bioinformatics & Biodiversity platform] ( http://mbb.univ-montp2.fr "MBB")


References
--------

Tsagkogeorga G., Cahais V. & Galtier N. 2012. The population genomics of a fast evolver: high levels of diversity, functional constraint and molecular adaptation in the tunicate Ciona intestinalis. Genome Biology and Evolution 4:740-749.

Gayral P., Melo-Ferreira J., Glémin S., Bierne N., Carneiro M., Nabholz B., Lourenço J.M., Alves P.C., Ballenghien M., Faivre N., Belkhir K., Cahais V., Loire E., Bernard A. & Galtier N. 2013. Reference-free population genomics from Next-Generation transcriptome data and the vertebrate-invertebrate gap. PLoS Genetics 9:e10003457.




- - - -

.alr file format

An alr file consists of several blocks, each corresponding to a contig. Each block includes a contig name line, a header line, and a series of data lines – one per nucleotide position.

The contig name line is like:

>contig_name

The header line is like:

maj M/P indiv1  indiv2  indiv3  …

where indiv1, indiv2, … are the identifiers of distinct individuals.

"maj" is for majoritary. It corresponds to the base most frequently read across individuals at the current position. M/P is for "monomorphic/polymorphic". It says whether, at the current position, only reads of the "maj" kind were observed (monomorphic position), or at least one read different from "maj" was found (polymorphic position).

Data lines contain the information of read counts per individual per position. They should be either like (monomorphic position):

G   M   13  4   7   ...

or like (polymorphic position):

T   P   23[0/0/0/23]    31[1/0/0/30]    18[0/11/0/7]    …

The first entry of a data line is the majoritary base. The second entry should be either M (monomorphic) or P (polymorphic). Each of the following entries corresponds to one individual. In the case of monomorphic positions, just the total number of reads is given. In the case of polymorphic positions, the total number of reads is given too, followed by the detailed read counts in A/C/G/T order.

The first data line above corresponds to a position at which individuals 1, 2 and 3 have 13, 4 and 7 reads of G, respectively, and no read of A, C or T. The second data line above corresponds to a position at which individual 1 has 23 reads of T, individual 2 has 1 read of A and 30 reads of T and individual 3 has 11 reads of C and 7 reads of T.


