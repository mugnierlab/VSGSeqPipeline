# VSG-Seq Analysis Pipeline

## Basic Usage
Here's how you use it:  
```
python VSG-Seq_Pipeline.py -i input_file_data.txt 
```
You can adjust various parameters

## Required Software

You'll need biopython (we use [anaconda](https://anaconda.org/anaconda/python)) and the following software installed and in your PATH:  

* [MULTo](http://sandberg.cmb.ki.se/multo/)  
* [blastn](https://www.ncbi.nlm.nih.gov/books/NBK279671/)  
* [trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)  
* [cd-hit](https://github.com/weizhongli/cdhit)  
* [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/releases)  (settings are currently optimized for v2.4.0)  

## BLAST Databases for Identifying VSGs

You can use any reference you want to identify VSGs. We have a few options available in [VSG_blastdbs/](VSG_blastdbs/).

There are three different VSG databases:
	* EATRO1125 VSGs ( )
	* Lister427 VSGs ( )
	* Combined database of BOTH Lister427 and EATRO1125

There is one 'NonVSG' database. This database has been cobbled together after multiple iterations of assembling expressed VSGs, inspecting them by hand and identifying common false positives (e.g., sometimes ESAGs assemble and this can filter those out).\
If you run the pipeline using this filter (the default), you'll need this database available for BLAST.

The fasta files from which these were created are available in [VSG_blastdbs/source_fasta](VSG_blastdbs/source_fasta)

## Output Files
