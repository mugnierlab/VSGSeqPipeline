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
	if filesList != "": # take in a text file with the sequence file dir/name
		file = open(filesList,"r")
		seqfiles = []
		for line in file.readlines():
			seqfiles.append(line)
		files = seqfiles
	
	trinityfiles = []
	
	for file in files:
		filenameS = file.split('/')# temp change for current file run
		if len(filenameS) == 2:
			filenameS = filenameS[1].split('.')[0] 
		else:
			filenameS = file.split('.')[0]
		trinityfiles.append(filenameS)
	
	return trinityfiles

file_data = makeFilesList("", "filesToSubmit_orig.txt")

print file_data