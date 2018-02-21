# VSG-Seq Analysis Pipeline
A pipeline for analyzing VSG-seq data. More information about the approach can be found in our [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4514441/). 

## Basic Usage 
```
python VSGSeqPipeline.py -i input_file_data.txt 
```
You can adjust various parameters for your run. Read about all of the options using `python VSGSeqPipeline.py -h`.
`input_file_data.txt` is a tab-delimited file describing your input files and the samples they came from. More info in [Input Files](#input-files).

## Required Software

You'll need biopython (we use [anaconda](https://anaconda.org/anaconda/python)) and the following software installed and in your PATH:  

* [trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)  
* [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/releases)  (settings are currently optimized for v2.4.0)
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (required for Trinity)
* [blastn](https://www.ncbi.nlm.nih.gov/books/NBK279671/)  
* [cd-hit](https://github.com/weizhongli/cdhit)
* [bowtie](http://bowtie-bio.sourceforge.net/index.shtml)
* [MULTo](http://sandberg.cmb.ki.se/multo/) (in your current working directory *or* indicate where it lives with `-p`)

## Input Files

 In addition to your sequencing files in FASTQ format (this pipeline uses single-end sequencing reads), you'll need a tab-delimited file containing information about your samples. The first line is a header containing whatever attributes of your sample you'd like to incorporate into downstream analysis. These attributes will be incorporated into the final file describing the results. The first column must be the name of each input FASTQ. [Here's](Examples/filesToSubmit.txt) an example. The program expects those files to be in the current working directory.

## BLAST Databases for Identifying VSGs

You can use any reference you want to identify VSGs. We have a few options available in [VSG_blastdbs](VSG_blastdbs).

There are three different VSG databases:
* EATRO1125 VSGs (`EATRO1125_vsgs`)
* Lister427 VSGs (`tb427_vsgs`)
* Combined database of BOTH Lister427 and EATRO1125 (`concatAntattb427`)

There is one 'NonVSG' database (`NOTvsgs`). This database has been cobbled together after multiple iterations of assembling expressed VSGs, inspecting them by hand, and identifying common false positives (e.g., certain ESAGs assemble frequently and this will filter those out). If you run the pipeline using this filter (the default), you'll need this database available for BLAST.  

These files need to be in your working directory *or* your `blastn` installation needs to be [configured](http://telliott99.blogspot.com/2009/12/blast-ncbirc-file.html) such that it can find them whereever they live on your machine.

The fasta files these were created from are available in [VSGdb_fasta](VSGdb_fasta).

## Output Files

All intermediate files produced in the pipeline are saved in one folder. A summary file shows the expression of each VSG in each sample, both in terms of RPKM (calculated using MULTo) and percentage of the population (RPKM for that VSG/total RPKM). If you assembled VSGs from your reads, it will also contain information on how similar those VSGs are to VSGs in your reference database. See an example [here](Examples).
