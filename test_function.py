#!/usr/bin/env python

from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
from Bio import SeqIO 
from Bio import Seq
from sys import argv
import subprocess
import argparse
import time
import os

def makeFilesList(files, filesList):
	if filesList != "": # take in a text file with the sequence file dir/name in first column. tab-separated, with headers line describing all variables
		file = open(filesList,"r")
		seqfiles = []
		header_list = []
		firstline = 0
		sampledict = {}
                filebasenames = []
                
		for line in file.readlines():
			if firstline == 0:
				line = line.strip('\n')
                                header_list = line.split('\t')[0:]
				firstline = 1
			else:
				line = line.strip('\n')
                                linesplit = line.split('\t')
				sampledict[linesplit[0].split('.')[0]] = line.split('\t')[1:]
				seqfiles.append(linesplit[0])
                                filebasenames.append(linesplit[0].split('.')[0])
	
	
	
	
	return (filebasenames, seqfiles, header_list, sampledict)
    
file_data = makeFilesList("", "filesToSubmit.txt")

filebasenames = file_data[0]
seqfiles = file_data[1]
headers = file_data [2]
sampledict = file_data[3]
print filebasenames
print seqfiles
print headers
print sampledict