# VSG-Seq Analysis Pipeline

## Basic Usage
Here's how you use it:  
```
python VSG-Seq_Pipeline.py -i input_file_data.txt 
```
You can adjust various parameters

## Required Software

You'll need biopython (we use [anaconda](https://anaconda.org/anaconda/python)) and the following software installed and in your PATH:  

* [trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)  
* [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/releases)  (settings are currently optimized for v2.4.0)
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (required for Trinity)
* [blastn](https://www.ncbi.nlm.nih.gov/books/NBK279671/)  
* [cd-hit](https://github.com/weizhongli/cdhit)
* [bowtie](http://bowtie-bio.sourceforge.net/index.shtml)
* [MULTo](http://sandberg.cmb.ki.se/multo/)  

## Input Files

 In addition to your sequencing files in FASTQ format (this pipeline uses single-end sequencing reads), you'll need a tab-delimited file containing information about your samples. The first line is a header line, containing whatever attributes of your sample you'd like to incorporate into downstream analysis. These attributes will be incorporated into the final file describing the results. The first column is the filename for each file. See xxx for an example of one such file.

## BLAST Databases for Identifying VSGs

You can use any reference you want to identify VSGs. We have a few options available in [VSG_blastdbs/](VSG_blastdbs/).

There are three different VSG databases:
* EATRO1125 VSGs ('EATRO1125_vsgs')
* Lister427 VSGs ('tb427_vsgs')
* Combined database of BOTH Lister427 and EATRO1125 ('concatAntattb427')

There is one 'NonVSG' database ('NOTvsgs'). This database has been cobbled together after multiple iterations of assembling expressed VSGs, inspecting them by hand and identifying common false positives (e.g., sometimes ESAGs assemble and this can filter those out).\
If you run the pipeline using this filter (the default), you'll need this database available for BLAST.

The fasta files from which these were created are available in [VSGdb_fasta](VSGdb_fasta).

## Output Files

