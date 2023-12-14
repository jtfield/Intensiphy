# Intensiphy

Intensiphy is an expansion of function wrapper for the [Extensiphy pipeline](https://github.com/McTavishLab/extensiphy/tree/main). With Intensiphy, you can specify the organism you wish to add to your starting phylogeny and raw genomic data will be found in the [NCBI Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra). Those sequences will be assembled with Extensiphy and added to your phylogeny automatically. If your query organism is used for too much broad research that you don't want included in your phylogeny, you can input a csv file of SRA IDs that you would specifically like to include and Intensiphy will do the rest!

[Setup and Use](#setup-and-use)

[Intensiphy Controls and Flags](#intensiphy-controls-and-flags)

[Output Files](#output-files)

[Phylogenetic Estimation](#phylogenetic-estimation)

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
To help explain some of the jargon (technical words and terms) that goes along with bioinformatics programs, we've written some tutorials, packaged with our Extensiphy program.
* The [command line tutorial](https://github.com/McTavishLab/extensiphy/blob/main/tutorial/command_line_tutorial.md) will help you get a grasp on how to find files in your computer using the shell/terminal/command line (you'll be a hacker in no time!).
* The [suffix tutorial](https://github.com/McTavishLab/extensiphy/blob/main/tutorial/suffix_tutorial.md) will help clarify the read suffix arguments.


### Quick test run
If you have followed one of the install approaches above, you are now ready to try a test run!  
We'll use the `ip_combo.fas` alignment file as our starting alignment. `ip_combo.fas` can be found in:
```
/intensiphy/testdata/ip_combo.fas
```

Now from the directory where you installed Intensiphy, run:

```bash
./intensiphy.py --align_file ./testdata/combo.fas --accession_method USER_INPUT --accs_file ./testdata/sra.csv 
```
This is a simple run, which downloads and assembles sequences found in the `sra.csv` file.
* The `--align_file` flag provides the path to the existing alignment to update.
* The `--accession_method` flag specfies what method accession numbers will be provided or collected by the program. Here we have collected and curated a file of NCBI SRA numbers.
* The `--accs_file` flag provides the path to the accession numbers file. The file must be in `.csv` format

Once Intensiphy has finished running on the test data, you should see a line saying:
```
Assembled sequences are found in /project/Intensiphy/ip_output/sequence_storage/

Number of sequences assembled during this run: 5

Number of sequences in the current sequence library: 10

```
* If you did not get this message, you'll have to check output log `ip_dev_log.txt`
to learn more about the issue before proceeding.  

We just constructed a sequence library using our original alignment and added 5 new sequences to the library.


## Intensiphy Controls and Flags:

### Required flags
```
- (--align_file) alignment in fasta format,
```

#### Either
```
- (--organism) the name or taxon ID of the organism you wish to collect sequences for from the NCBI SRA (Example: Neisseria gonorrhoeae[Organism] or txid482),
```

#### Or
```
- (--accession_method) Dictates how collecting and inputting accession numbers will be handled. (OPTIONS: USER_INPUT and AUTO_DL), (DEFAULT:AUTO_DL),
- (--accs_file) accession file. Used if you are using the [--accession_method USER_INPUT] flag to pass in a curated file of SRA numbers,
```

### Optional flags
```
- (--cores) number of cores allocated to Intensiphy for the alignment and assembly steps (DEFAULT: 2),
- (--ref) reference sequence label (without suffix or file ending information). (Example: SRR1500345) (DEFAULT: random selection)
- (--placement) Toggles the phylogenetic placement function once all accession files have been downloaded and assembled. Toggles automatic phylogenetic estimation of starting tree if one was not input using the [--starting_tree] flag. (OPTIONS: ON, OFF) (DEFAULT: OFF)
- (--starting_tree) a phylogeny produced from the input alignment file. Used for phylogenetic placement.
- (--ip_out_dir) path and folder name you would like to use to store the outputs of this program (DEFAULT: ip_output).
```

## Output Files
* Single sequence files: found in your sequence storage folder
```
[OUTDIR]/sequence_storage/[SEQUENCE_FOLDERS]
```

## Phylogenetic Estimation
Explanation of Intensiphy's interaction and utility with phylogenetic estimation coming soon!


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