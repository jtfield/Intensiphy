# Intensiphy

Intensiphy is an expansion of function wrapper for the [Extensiphy pipeline](https://github.com/McTavishLab/extensiphy/tree/main). With Intensiphy, you can specify the organism you wish to add to your starting phylogeny and raw genomic data will be found in the [NCBI Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra). Those sequences will be assembled with Extensiphy and added to your phylogeny automatically. If your query organism is used for too much broad research that you don't want included in your phylogeny, you can input a csv file of SRA IDs that you would specifically like to include and Intensiphy will do the rest!

[Setup and Use](#setup-and-use)

[Building and testing your own Intensiphy Docker image](#building-and-testing-your-own-intensiphy-docker-image)

[Intensiphy Controls and Flags](#intensiphy-controls-and-flags)

[Output Files](#output-files)

[Phylogenetic Estimation](#phylogenetic-estimation)

[Additional Software](#additional-software)

[Dependencies](#dependnecies)

[Reporting Problems](#reporting_problems)

## Setup and Use

### Mamba
The simplest and most hassle free way to install the dependencies of Intensiphy is using [Mamba](https://github.com/mamba-org/mamba). The [Mamba Installation](https://github.com/McTavishLab/extensiphy/blob/main/tutorial/alternative_installation_methods.md#installing-dependencies-with-anaconda) section of this repository will walk through this process in more detail.

### Advanced
If you're comfortable installing programs by hand, the [Advanced Installation Methods](https://github.com/McTavishLab/extensiphy/blob/main/tutorial/alternative_installation_methods.md#installing-depedencies-by-hand) section is for you. This is largely only tested on Linux (Ubuntu) operating systems.

### Intensiphy Tutorial
We recommend you run through the [Intensiphy tutorial](https://github.com/McTavishLab/extensiphy/blob/main/tutorial/extensiphy_tutoria.md) for a more in-depth walkthrough of Intensiphy's features. The tutorial will walk through how to run Intensiphy using different data types and options. You can copy code snippets into your terminal window.

### Additional Tutorials
To help explain some of the jargon (technical words and terms) that goes along with bioinformatics programs, we've also included some other tutorials.
* The [command line tutorial](https://github.com/McTavishLab/extensiphy/blob/main/tutorial/command_line_tutorial.md) will help you get a grasp on how to find files in your computer using the shell/terminal/command line (you'll be a hacker in no time!).
* The [suffix tutorial](https://github.com/McTavishLab/extensiphy/blob/main/tutorial/suffix_tutorial.md) will help clarify the read suffix arguments.


### Quick test run
If you have followed one of the install approaches above, you are now ready to try a test run!  
We'll use the `combo.fas` alignment file as our starting alignment. `combo.fas` can be found in:
```
/extensiphy/testdata/combo.fas
```

Now, either from Mamba env, or from the directory where you installed Extensiphy, run:

```bash
./extensiphy.sh -a ./testdata/combo.fas -d ./testdata -1 _R1.fq -2 _R2.fq -u PHYLO -o EP_output
```
This is a simple run on three paired end read samples, which are found in the directory `extensiphy/testdata`
* The `-a` flag provides the path to the existing alignment to update.
* The `-d` flag provides the path to your directory of fastq files.
* The `-1` and `-2` flags specify the filename endings for each of the readfiles. (defaults are `_R1.fq` and `_R2.fq`, more info at https://github.com/McTavishLab/extensiphy/blob/main/tutorial/suffix_tutorial.md)
* The `-u` flag specfies what analysis to run. Here we are building a phylogeny. (default is `ALIGN`, building an alignment only.)
* The `-o` flag specifies the output directory. (default is `EP_output`)

Once Extensiphy has finished running on the test data, you should see a lines saying:
```
Alignment file is: /project/extensiphy/EP_output/RESULTS/extended.aln

Tree file is: /project/extensiphy/EP_output/RESULTS/RAxML_bestTree.consensusFULL

```
* If you did not get this message, you'll have to check output log `ep_dev_log.txt`
to learn more about the issue before proceeding.  

We just added 3 new taxa to a starting multiple sequence alignment and obtained a tree that includes these new taxa.


* If you are using docker - exit the container by typing
```
exit
```

* You can copy the extended tree to your local directory using:

```
docker cp ep_container:/project/extensiphy/EP_output/RESULTS/RAxML_bestTree.consensusFULL .
```

* For a deeper walk through, take a look through the [tutorial](https://github.com/McTavishLab/extensiphy/blob/main/tutorial/extensiphy_tutoria.md).

* To get right down to business and update your own alignment, continue to the next section.


### Using Extensiphy on your own data.

We'll use brackets `[]` to indicate variables you should replace with your own files or paths.
Replace the `[stuff inside the brackets]` with the appropriate paths and folder names you've used so far.


If you have installed Extensiphy locally, you can just pass in the paths to your data, and run the analysis.

````
./extensiphy.sh -a [path to your_input_alignment] -d [path to your_directory_of_reads] -1 [r1_suffix] -2 [r2_suffix] -u [either PHYLO or ALIGN, depending on if you want a phylogeny or just and alignment] -o [your_output_dir]
````

If you are using docker, it is simplest to link your data directory to a new container.

Put the input alignment and raw reads you want to align in a directory. e.g. [my_data_dir]

We'll build a new Extensiphy Docker container and connect the directory containing your data to the container.

```bash
docker run --name ep_container_link -i -t -v [/path/to/my_data_dir]:/project/linked_data mctavishlab/extensiphy bash
```

This shares the 'my_data_dir' folder between your operating system and the docker container. (In this example it is named "my_data_dir" locally and "linked_data" in your docker container, but you can name them the same thing in both places if you prefer.)

Now you can run `extensiphy.sh` as earlier but we'll specify the directory where your data is located.

```bash
./extensiphy.sh -a /project/linked_data/[alignment_file] -d /project/linked_data -1 [suffix_1] -2 [suffix_2] -o linked_data/[output_dir_name]
```

By putting the outputs into the linked directory, you can access them directly through your operating system without having to copy them.


## Extensiphy Controls and Flags:

### Required flags
```
- (-a) alignment in fasta format,
- (-d) directory of paired end fastq read files for all query taxa,
- (-u) produce only an updated alignment or perform full phylogenetic estimation (ALIGN or PHYLO) (DEFAULT: ALIGN),
```

### Optional flags
```
- (-t) tree in Newick format produced from the input alignment that you wish to update with new sequences or specify NONE to perform new inference (DEFAULT: NONE),
- (-1, -2) suffix (ex: R1.fastq or R2.fastq) for both sets of paired end files. Required if suffix is different than default (DEFAULTS: R1.fq and R2.fq),
- (-m) alignment type (SINGLE_LOCUS_FILES, PARSNP_XMFA or CONCAT_MSA) (DEFAULT: CONCAT_MSA),
- (-o) directory name to hold results (DEFAULT: creates EP_output),
- (-r) Selected a reference sequence from the alignment file for read mapping or leave as default and the first sequence in the alignment will be chosen (DEFAULT: RANDOM),
- (-p) number of taxa to process in parallel,
- (-c) number of threads per taxon being processed,
- (-e) set read-type as single end (SE) or pair-end (PE) (DEFAULT: PE),
- (-g) output format (CONCAT_MSA or SINGLE_LOCUS_FILES) (DEFAULT: CONCAT_MSA),
- (-s) specify the suffix (.fa, .fasta, etc) (DEFAULT: .fasta),
- (-b) bootstrapping tree ON or OFF (DEFAULT: OFF)
- (-i) set whether to clean up intermediate output files to save disk space)(KEEP, CLEAN)(DEFAULT: KEEP)

 if using single locus MSA files as input,
- (-f) csv file name to keep track of individual loci when concatenated (DEFAULT: loci_positions.csv),
- (-n) Set size of locus minimum size cutoff used as input or output (Options: int number)(DEFAULT: 700)     
```

## Output Files
* Concatenated alignment file: found in your output folder
```
[OUTDIR]/RESULTS/extended.aln
```

* Phylogeny in newick file format (if you selected to output a phylogeny): found in your output folder
 ```
 [OUTDIR]/RESULTS/RAxML_bestTree.consensusFULL
 ```

* Taxon specific intermediate files (if you kept intermediate files): found in your output folder
 ```
 [OUTDIR]/[TAXON_NAME]
 ```
  .sam, .bam and .vcf files can be found in here for any additional analyses.


## Phylogenetic Estimation
Extensiphy is targeted towards producing an updated sequence alignment and allowing users to use the alignment with any phylogenetic estimation method they choose. We provide a phylogenetic estimation as a convenience but you are in no way locked into using this estimation method. You can simply take the [alignment output by Extensiphy](#output-files) and use that alignment as input for your favorite estimation method. A few notes:

* Currently, phylogenetic estimation with Extensiphy is performed by RAxML using the GTR model. This setting cannot be changed at this time.

* To avoid estimating a phylogeny using the packaged RAxML program and settings, use the `-u ALIGN` option when running Extensiphy.

* To use an alternative method of phylogenetic estimation, when an Extensiphy run is complete, the `[OUTDIR]/RESULTS/extended.aln` file should be used as input for your chosen estimation method.

* If you wish to use multiple single locus alignment files as input to another estimation method, please see the [tutorial](https://github.com/McTavishLab/extensiphy/blob/main/tutorial/extensiphy_tutoria.md#outputting-updated-single-locus-alignment-files) for more information on updating single locus alignments.

## Additional Software
Extensiphy is the primary program of this software package.
However, another piece of software is included: Gon\_phyling.
Gon\_phyling is a piece of software for building starting alignments and phylogenies when you only have raw-read fastq files.
Gon\_phyling isn't the focus software but we provide it in case you might find it useful.
Checkout the [program and README](https://github.com/McTavishLab/extensiphy/tree/main/gon_phyling) in the `gon_phyling` directory.


## Dependencies

Dependencies (Separate programs you'll need to install):

1. [Python 3](https://www.python.org/)
2. [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)
3. [RAxMLHPC](https://github.com/stamatak/standard-RAxML)
4. [Seqtk](https://github.com/lh3/seqtk)
5. [Samtools](http://www.htslib.org/)
6. [Bcftools](http://www.htslib.org/)
7. [Fastx toolkit](http://hannonlab.cshl.edu/fastx_toolkit/download.html)
8. [Dendropy](https://dendropy.org/)

## Reporting Problems
Software will have bugs. We try to address issues with Extensiphy as they arise.
If you run into an issue, please report it using Extensiphy's [Issue Tracker](https://github.com/McTavishLab/extensiphy/issues).
You can also search the Issue Tracker for solved fixes for previously identified issues.
Finally, you can contact us at jtoscanifield@ucmerced.edu to discuss any problems with installing or running Extensiphy.