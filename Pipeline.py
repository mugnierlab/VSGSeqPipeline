from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
from Bio import SeqIO 
from Bio import Seq
from sys import argv
import VSGFunctions as vsgf
import subprocess
import argparse
import time
import os

# usage is ORF_finder.py infile minimum protein length outfile
# this goes through a FASTA file and finds open reading frames - as soon as 
# an ORF is found, it writes that record to a new FASTA file and stops looking in
# that record.
#if no start is found, beginning -> stop is an ORF
#if no stop is found, start -> end is an ORF

# finds all open reading frames within a specified length

parser = argparse.ArgumentParser()

# arguments about input files and naming
parser.add_argument('-s', nargs= '+',help='listed sequence files to be put through the pipeline, fq', action="store", dest="s")
parser.add_argument('-st',help='text file with .fastq names of sequence files to be put through the pipeline, same for entering any step of pipeline, fq', action="store", dest="st", default='')
parser.add_argument('-d', help='additional descriptive terms to name your run', action="store", dest='d', default='')
parser.add_argument('-stderr', help='0, direct standard error to an output file. 1, output to console', action="store", dest='stderr', default=0)
# trimming setting
parser.add_argument('-g', help='trim_galore, stringency for trim galore', action ="store", dest = "g", default="3") 
parser.add_argument('-trim', help='cutadapt/trim_galore, minimum length of read for trimming', action ="store", dest = "trim", default="50") 
# trinity settings
parser.add_argument('-minp', help='Trinity, minimum protein length you are filtering for', action ="store", type=int, dest = "minp", default=300) 
parser.add_argument('-mem', help='Trinity, max memory allocation for trinity, G', action ="store", dest = "mem", default="10") 
parser.add_argument('-cpu', help='Trinity/cd-hit-est/MULTo/bowtie, number of processors to use', action="store", dest='cpu', default='2')
# Blast settings
parser.add_argument('-vsgdb', help='BLAST, name of the vsg database to blast against', action ="store", dest = "vsgdb", default="tb427_vsgs")
parser.add_argument('-NoNonVSGdb', help='BLAST, adding this flag means you dont want to BLAST ORFS againat the non-VSG database' , dest='NoNonVSGdb', action='store_true', default=False)
#parser.set_defaults(NoNonVSGdb=False)
# cd-hit-est parameters
parser.add_argument('-sit', help='cd-hit-est, sequence identiy threshold - how much the alignment has to match, percentage. value is 0.0 through 1.0 ', action ="store", dest = "sit", default=".98")
parser.add_argument('-t', help='cd-hit-est, number of threads for cd-hit-est to use, 0 default(all CPU will be used)', action ="store", dest = "t", default="0")
# MULTO settings
parser.add_argument('-p', help='MULTo, path to MULTo1.0 folder. default is /Users/mugnierlab/, please dont use "~/", python doesnt like this in the path', action="store", dest='p', default='/Users/mugnierlab/') # default assumes MULTo is in your home dir
parser.add_argument('-v', help='MUTLo(bowtie step), number of mismatches allowed', action="store", dest='v', default='2')
parser.add_argument('-reuseMulto', help='MUTLo, name of the multo files to be resued, typically same as header, default will make new MULTo files', action="store", dest='rmulto', default='')
#where the pipeline will start
parser.add_argument('-start', help='the step you want the pipeline to start at. 0 = input is raw untrimmed data. 1 = start after trimming, input is already trimmed. 2 = Start after Trinity, input are trinity files. 3 = start after finding ORF, input are ORFs. 4= Start after BLAST/cd-hit-est, input are VSG/ORFs. 5 = continues till the end to MULTo', type=int, action="store", dest='start', default=0)
parser.add_argument('-header', help='input header variable in the format of Y-m-d-H_M , include descriptive headers too ', action="store", dest='head', default='')
#where the pipeline will stop
parser.add_argument('-stop', help='the step you want the pipeline to stop at. 1 = stop after trimming. 2 = Stop after Trinity. 3= stop after finding ORF. 4= Stop after BLAST, 5 = continues till the end to MULTo', type=int, action="store", dest='stop', default=5)

# add segment to tkae sequence data, trim with trim_galore and cutadapt

# run trimmed sequences trough trinity

# then take the trimmed trinity sequences and run through this pipeline

arguments = parser.parse_args()

header = time.strftime("%Y-%m-%d-%H_%M")  # Year/Month/Day-Hour:Minute , names the folder for output files
if arguments.d != '':
	for d in arguments.d:
		header = header + "-"+ str(d)
if header != "":
	header = arguments.head
if arguments.stderr == 0:
	if not os.path.exists(header+"/StandardError"):
		os.makedirs(header+"/StandardError")

# start from the very beguinning 
if arguments.start == 0:
	if not os.path.exists(header):
		os.makedirs(header) # creates the folder
	
	trinityfiles =  vsgf.makeFilesList(arguments.s, arguments.st)
	vsgf.trimSequences(header, trinityfiles, arguments)
	if arguments.stop >1:
		vsgf.trinity(header, trinityfiles, arguments)
		if arguments.stop > 2:
			vsgf.findORFs(header, trinityfiles, arguments)
			if arguments.stop > 3:
				vsgf.blastCDHIT(header, trinityfiles, arguments)
				if arguments.stop > 4:
					vsgf.makeMulto(header, trinityfiles, arguments)
# start with trimmmed files
elif arguments.start == 1:
	trinityfiles =  vsgf.makeFilesList(arguments.s, arguments.st)
	if arguments.stop >1:
		vsgf.trinity(header, trinityfiles, arguments)
		if arguments.stop > 2:
			vsgf.findORFs(header, trinityfiles, arguments)
			if arguments.stop > 3:
				vsgf.blastCDHIT(header, trinityfiles, arguments)
				if arguments.stop > 4:
					vsgf.makeMulto(header, trinityfiles, arguments)
# start with trinity files
elif arguments.start == 2:
	trinityfiles =  vsgf.makeFilesList(arguments.s, arguments.st)
	if arguments.stop > 2:
		vsgf.findORFs(header, trinityfiles, arguments)
		if arguments.stop > 3:
			vsgf.blastCDHIT(header, trinityfiles, arguments)
			if arguments.stop > 4:
				vsgf.makeMulto(header, trinityfiles, arguments)
# start with ORFs
elif arguments.start == 3:
	trinityfiles =  vsgf.makeFilesList(arguments.s, arguments.st)
	if arguments.stop > 3:
		vsgf.blastCDHIT(header, trinityfiles, arguments)
		if arguments.stop > 4:
			vsgf.makeMulto(header, trinityfiles, arguments)
# start with BLAST-ed VSG ORFs
elif arguments.start == 4:
	trinityfiles =  vsgf.makeFilesList(arguments.s, arguments.st)
	if arguments.stop > 4:
		vsgf.makeMulto(header, trinityfiles, arguments)
else:
	print("Bad start argument, double check the start number you entered")













