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

def makeFilesList(filesList):
	#if filesList != "": # take in a text file with the sequence file dir/name in first column. tab-separated, with headers line describing all variables - MUST HAVE FQ suffix(for now?)no periods in middle
	file = open(filesList,"r")
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
			sampledict[linesplit[0].split('.')[0]] = line.split('\t')[0:]
			filebasenames.append(linesplit[0].split('.')[0])
	return (filebasenames, header_list, sampledict)

def trimSequences(header, filebasenames, arguments):
	stringency = arguments.g
	for file in filebasenames:
		stderr_tg = ""
		stderr_ca = ""
		if arguments.stderr == 0 :
			stderr_tg = " 2> " + header + "/StandardError/trim_galore-"+file+".txt"
			stderr_ca = " 2> " + header + "/StandardError/cutadapt-"+file+".txt"
		subprocess.call(['trim_galore --stringency '+stringency+' --length '+arguments.trim+' --dont_gzip --output_dir ' + header + "/ " +str(file) +".fq"+stderr_tg], shell=True) # trim off sequenceing adapters
		subprocess.call(['cutadapt -m '+arguments.trim+' -b ATTTAGGTGACACTATAG -b CTATAGTGTCACCTAAAT '+ header + "/"+str(file)+'_trimmed.fq > '+ header + "/"+str(file)+'_trimmed2.fq'+stderr_ca], shell=True) # trim off SP6 sequences (from VSG PCR step)
		subprocess.call(['rm '+ header + "/"+str(file)+'_trimmed.fq'], shell=True) # removes intermediate trimmed file 


def trinity(header, filebasenames, arguments):
	max_memory_trinity  =arguments.mem
	min_pro_len = arguments.minp
	trim = arguments.start
	for file in filebasenames:
		if trim == 1:
			filepath = str(str(file) +'.fq')
		if trim == 0:
			filepath = str(header + "/"+str(file) +'_trimmed2.fq')
		stderr_tr = ""
		if arguments.stderr == 0 :
			stderr_tr = " > " + header + "/StandardError/Trinity-"+file+".txt"
		subprocess.call(['Trinity --seqType fq --max_memory '+max_memory_trinity+'G --full_cleanup --single ' + str(filepath) + ' --CPU '+arguments.cpu+' --min_contig_length ' + str(min_pro_len*3) + ' --no_path_merging --no_normalize_reads --output ' + header+"/"+file+"_Trinity/" + stderr_tr], shell=True)

def addSeqRecord(recordRD, start, end, count):
	sequence = recordRD.seq[start:end]
	ORF_outfile.write('>'+str(recordRD.id)+'_'+str(count)+'_'+str(start)+'_'+str(end)+"_" + filename + '\n'+str(sequence)+'\n')
	trans_out_file.write('>'+str(recordRD.id)+'_'+str(count)+'_'+str(start)+'_'+str(end)+"_" + filename + '\n'+str(sequence.translate())+'\n')

def addSeqRecord_RC(recordRD, start, end, count):
	sequence = recordRD.seq[start:end].reverse_complement()
	ORF_outfile.write('>'+str(recordRD.id)+'_'+str(count)+'_'+str(start)+'_'+str(end)+'_RC_' + filename +'\n'+str(sequence)+'\n')
	trans_out_file.write('>'+str(recordRD.id)+'_'+str(count)+'_'+str(start)+'_'+str(end)+"_RC_" + filename + '\n'+str(sequence.translate())+'\n')

def findORFs(header, filebasenames, arguments):
	min_pro_len = arguments.minp
	for file in filebasenames:
		record_dict = SeqIO.index(header+ "/" + file+"_Trinity.Trinity.fasta","fasta") # "parses" a fasta file, creating a dictionary-like object of sequences. not everything is kept in memeory. instead it just records where each record is within the file. parses on demand. 
		# the key is the dictionary is the ">" line in the fasta file
		 # minimum protein length, within typical VSG protein length, in a.a.
		global ORF_outfile
		global contig_outfile
		global trans_out_file
		ORF_outfile = open(os.path.join(header, header+"_"+file+"_orf.fa"), 'w') # each file will get their own ORF output file
		contig_outfile = open(os.path.join(header, header+"_"+file+"_contig.fa"), 'w') # same as above, but for contigs the ORFs came from
		trans_out_file = open(os.path.join(header, header+"_"+file+'_orf_trans.fa'), 'w') # translated orf file

		global filename
		filename = file # takes only the name of the file, this will be added to the ORF description, so we know where it came from
		for record in record_dict: # iterates through the sequences
			#print record_dict[record]
			
			if len(record_dict[record].seq) > ((min_pro_len)*3): # if the orf length is > 3 times the specified protein length 
				count=1 
				for strand, nuc in [(+1, record_dict[record].seq), (-1,record_dict[record].seq.reverse_complement())]: # get a list of possible ORF translations, loops twice, each tuple
					for frame in [0,1,2]: #loops 3 times, frame will start at 0, 1, and 2
						wholeAdded = False # will allow for the read to be put in 3 times(once for each strand)
						
						trans = nuc[frame:].translate() # translates the sequence nuc[frame] into protein sequence for the current frame
						# find a way for it to stop at first stop codon
						trans_len = len(trans) # gets length of the protein sequence
						seq_len = len(record_dict[record].seq) # gets the length of the nucleotide sequence from the fasta file 						
						trans_end = trans.find("*", 0)
						trans_start = trans.find("M", 0)
						if trans_start != -1: #start codon exists
							if trans_end != -1: # end codon exists
								if trans_start < trans_end: # there is a start codon before the 1st stop codon
									if trans_end - trans_start > min_pro_len: # check if this protein is long enought to be considered
										#print "A"
										if strand ==1:
											addSeqRecord(record_dict[record], frame+(trans_start*3), frame+(trans_end*3)+3, count)
										else:
											addSeqRecord_RC(record_dict[record], max(frame, seq_len-frame-(trans_end*3)-3), seq_len-frame-(trans_start*3), count)
										count += 1
									elif trans_end > min_pro_len : # protein isn't long enough, so just take the entire front end if thats long enough
										#print "B"
										if strand ==1:
											addSeqRecord(record_dict[record], frame, frame+(trans_end*3)+3, count)
										else:
											addSeqRecord_RC(record_dict[record], max(frame, seq_len-frame-(trans_end*3)-3), seq_len, count)
										count += 1
									# if neither protein is long enough,or even if it was, we will move on to find where the next start codon is, and find the end codon for that one
									trans_start = trans.find("M", trans_end)
									trans_end = trans.find("*", trans_start)		
								else: # there is a start codon AFTER the first appearing stop codon
									if trans_end > min_pro_len: # is the seq_start to end codon is long enough for the a min length protein
										#print "C"
										if strand ==1:
											addSeqRecord(record_dict[record], frame, frame+(trans_end*3)+3, count)
										else:
											addSeqRecord_RC(record_dict[record], max(frame, seq_len-frame-(trans_end*3)-3), seq_len, count)
										count += 1
									trans_end = trans.find("*", trans_start) # will start checking the orf at that first found start codon, sets new trans_end to stop codon after the start	
							else: # start codon exists, but there is no stop codon anywhere! :O
							# this portion of the code will take the entire thing as an ORF, we are being generous
								if wholeAdded == False: # won't add if the entire strand was already added
									wholeAdded = True
									addSeqRecord_RC(record_dict[record], frame, seq_len, count)
									count += 1
							# this part isn't generous and only takes from the start codon till the end of the sequence
		#						if trans_len-trans_start > min_pro_len: # is the remaining protein segment after the start codon long enought to be a protein?
		#							print "D"
		#							print trans_start
		#							if strand ==1:
		#								addSeqRecord(record_dict[record], frame+(trans_start*3), seq_len, count)
		#							else:
		#								addSeqRecord_RC(record_dict[record], frame, seq_len-frame-(trans_start*3), count)
		#							count += 1
								trans_start = trans_len # won't go into loop, bc the entire rest of the sequence just got taken as an orf, so there is no point	in searching further, no more to search
						elif trans_end != -1: # start codon doesn't exist, but end codon does
							if trans_end > min_pro_len: # take entire sequence from start of seq to the found end codon as ORF if longer than min length
								#print "E"
								if strand ==1:
									addSeqRecord(record_dict[record], frame, frame+(trans_end*3)+3, count)
								else:
									addSeqRecord_RC(record_dict[record], max(frame, seq_len-frame-(trans_end*3)-3), seq_len-frame-(trans_start*3), count)
								count += 1	
							trans_start = trans_len # won't search for more, since there aren't any start codons anywhere
						else: # there are no start or stop codons anywhere?!
							if wholeAdded == False:
								wholeAdded == True
								addSeqRecord_RC(record_dict[record], frame, seq_len, count)
								trans_start = trans_len
								count += 1	
						
						trans_max = trans_len - min_pro_len # farthest a codon can be before it's below min protein length

						while trans_start < trans_max:
							if trans_start == -1: # if no more starts are found, doesn't matter if there is a stop or not, game over
								trans_start = trans_len
							elif trans_end == -1: #if start is found but no end is found
								if trans_len-trans_start > min_pro_len: #if ORF from found start to end sequence is long enough
									#print "G"
									if strand ==1:
										addSeqRecord(record_dict[record], frame+(trans_start*3), seq_len, count)
									else:
										addSeqRecord_RC(record_dict[record], frame, seq_len-frame-(trans_start*3), count)
									count += 1
								trans_start = trans_len
							elif trans_end-trans_start > min_pro_len: # if we have both a start and an end and its long enough
								#print "H"
								if strand ==1:
									addSeqRecord(record_dict[record], frame+trans_start*3, frame+(trans_end*3)+3, count)
								else:
									addSeqRecord_RC(record_dict[record], max(frame, seq_len-frame-(trans_end*3)-3), seq_len-frame-(trans_start*3), count)
								count += 1	
								trans_start = trans.find("M", trans_end) # finds index of start aa sequence
								trans_end = trans.find("*", trans_start)	
							else:
								trans_start = trans.find("M", trans_end) 
								trans_end = trans.find("*", trans_start)			

				if count > 1: # have any new orf been added? if so, add this record to file
					SeqIO.write(record_dict[record], contig_outfile, "fasta") 

		ORF_outfile.close()
		contig_outfile.close()
		trans_out_file.close()


def blastCDHIT(header, filebasenames, arguments):
	vsgdDbName = arguments.vsgdb
	seqIdenThresh = arguments.sit
	numCPU = arguments.cpu
	max_memory_cdhit = 1000*int(arguments.mem)
	scoredict = {}
	for file in filebasenames:
		filename = header + "/"+header+"_"+file+"_orf"
		print ' *****analyzing '+filename+".fa"+' *****'
		#blast VSG
		subprocess.call(['blastn -db '+vsgdDbName+' -query '+filename+'.fa -outfmt 5 -out '+filename+'.xml'], shell=True)
		#blast nonVSG
		if arguments.NoNonVSGdb == False:
			subprocess.call(['blastn -db NOTvsgs -query '+filename+'.fa -outfmt 5 -out '+filename+'_nonVSG.xml'], shell=True)
		#get all the blast results which are for ONLY VSGs, get rid of hits which are VSG-similar but not vsgs
		#blast_sort(filename+'.xml', filename+'_nonVSG.xml',filename+".fa", arguments.NoNonVSGdb) # the _VSGs.fa file is produced from this
		
		v = filename+'.xml'
		n = filename+'_nonVSG.xml' 
		s = filename+".fa"
		
		result_handle = open(v)
		
		blast_records = NCBIXML.parse(result_handle) # returns an iterator of the blast results
		record_dict = SeqIO.index(s,"fasta")
		
		outfile = open(s.split('.')[0]+'_VSGs.fa', 'w')
		#scorefile = open(s.split('.')[0]+'_VSGs_scores.txt', 'w')
	
		hit_list = []
		scoredict = {}
	
		if arguments.NoNonVSGdb == False: # blasted against nonVSG database
			nonVSGresult_handle = open(n)
			blast_records_nonVSG = NCBIXML.parse(nonVSGresult_handle)
			 # list of VSGs we have found! 
			exclude_list = []
					
			print 'Now looking for non-VSG transcripts...'
			for blast_record_nonVSG in blast_records_nonVSG:
				for alignment in blast_record_nonVSG.alignments:
					for hsp in alignment.hsps: 
						percent_identity = (100.0 * hsp.identities) / alignment.length # hsp.identities is a tuple(bp matches, total bp in seq) to give percent match of sequence, percent identity is # of bp
						percent_query_identity = (100.0 * hsp.identities) / blast_record_nonVSG.query_letters
						#print blast_record_nonVSG.query+'\t'+alignment.title+'\t'+str(percent_identity)+'\t'+str(percent_query_identity)+'\t'
						if (percent_query_identity > 30 and hsp.identities > 300) or (percent_identity > 90):
							if not blast_record_nonVSG.query in exclude_list:
								exclude_list.append(str(blast_record_nonVSG.query))
								#print 'nonVSG hit!'+'\t'+str(blast_record_nonVSG.query)+' \t '+str(alignment.title)
	
			print 'VSG hits! - maybe?'
			for blast_record in blast_records:
				for alignment in blast_record.alignments:
					for hsp in alignment.hsps:
						if hsp.expect < 1.0e-10: # hsp.expect = e value for the hsp value, the lower the e value, the more statistically significant 
							if not blast_record.query in hit_list: # if this query hasn't already been added to the hit list, add it now
								if not blast_record.query in exclude_list: # if the query isn't a fake VSG hit, add it now!
									hit_list.append(str(blast_record.query))											# percent query aligned										# percent identity
									#scoredict[str(blast_record.query)] = str('\t'+str(alignment.title)+'\t'+str((100.0 * hsp.identities) / blast_record.query_letters)+'\t'+str((100.0 * hsp.identities) / alignment.length)+'\t'+str(alignment.length)+'\n')
									SeqIO.write(record_dict[blast_record.query], outfile, "fasta")
		else: # didn't blast against nonVSG database
			for blast_record in blast_records:
				for alignment in blast_record.alignments:
					for hsp in alignment.hsps:
						if hsp.expect < 1.0e-10: # hsp.expect = e value for the hsp value, the lower the e value, the more statistically significant 
							if not blast_record.query in hit_list: # if this query hasn't already been added to the hit list, add it now
								hit_list.append(str(blast_record.query))											# percent query aligned										# percent identity
								#scoredict[str(blast_record.query)] = ('\t'+str(alignment.title)+'\t'+str((100.0 * hsp.identities) / blast_record.query_letters)+'\t'+str((100.0 * hsp.identities) / alignment.length)+'\t'+str(alignment.length)+'\n')
								SeqIO.write(record_dict[blast_record.query], outfile, "fasta")
		
		
	# cdhit merge, clusters VSGs which are similar into one, so that we dont have replicates of VSGs


	all_VSGs = open(os.path.join(header, header+'_orf_VSGs.fa'), 'w')
	for file in filebasenames:
		subprocess.call(['cat '+header + "/"+header+"_"+file+"_orf_VSGs.fa >> " +header + "/"+header+"_orf_VSGs.fa"], shell=True) # concatenates all the files for MULTo

	stderr_cd = ""
	if arguments.stderr == 0 :
		stderr_cd = " > " + header + "/StandardError/cdhitest.txt"
	subprocess.call(['cd-hit-est -i '+header+"/"+header+'_orf_VSGs.fa '+' -o '+header+"/"+header+'_orf_VSGs_merged.fa -d 0 -c ' + seqIdenThresh + ' -n 8 -G 1 -g 1 -s 0.0 -aL 0.0 -M '+str(max_memory_cdhit)+' -T ' + numCPU + stderr_cd], shell=True)
	
	#blast merged database:
	
	subprocess.call(['blastn -db '+vsgdDbName+' -query '+header+"/"+header+'_orf_VSGs_merged.fa -outfmt 5 -out '+header+"/"+header+'_orf_VSGs_merged.xml'], shell=True)
	
	#process blast
	result_handle = open(header+"/"+header+'_orf_VSGs_merged.xml', 'r')
	blast_records = NCBIXML.parse(result_handle)
	s = header+"/"+header+'_orf_VSGs_merged.fa'
	record_dict = SeqIO.index(s,"fasta")
	hit_list = []
	for blast_record in blast_records:
				for alignment in blast_record.alignments:
					for hsp in alignment.hsps:
						if hsp.expect < 1.0e-10: # hsp.expect = e value for the hsp value, the lower the e value, the more statistically significant 
							if not blast_record.query in hit_list: # if this query hasn't already been added to the hit list, add it now
								hit_list.append(str(blast_record.query))											# percent query aligned										# percent identity
								scoredict[str(blast_record.query)] = ('\t'+str(alignment.title)+'\t'+str((100.0 * hsp.identities) / blast_record.query_letters)+'\t'+str((100.0 * hsp.identities) / alignment.length)+'\t'+str(alignment.length)+'\n')
								#SeqIO.write(record_dict[blast_record.query], outfile, "fasta")
	
	return scoredict



def makeMulto(header, filebasenames, arguments):
	path = arguments.p
	rmulto = arguments.rmulto
	numCPU = arguments.cpu
	numMismatch = arguments.v
	currDir = os.getcwd()	
	if path != '':
		os.chdir(path)
	print os.getcwd()
	os.chdir('MULTo1.0')
	fullmultopath = os.getcwd()
	print fullmultopath

	# make multo dir
	# MULTo identifies, stores and retrieves the minimum length required at each genomic position to be unique across the genome or transcriptome.
	#make multo file heirarchy from given basename
	if rmulto == '': # make multo
		subprocess.call(['mkdir -p '+fullmultopath+'/files/tbb/tb'+header+'/fastaFiles/annotationFiles/'], shell=True)
		subprocess.call(['mkdir -p '+fullmultopath+'/files/tbb/tb'+header+'/fastaFiles/genomeFasta/noRandomChrom'], shell=True)

		# make bed
		record_dict = SeqIO.index(currDir +'/' + header + "/"+header+"_orf_VSGs_merged.fa","fasta")
		concatFA = open('chr1.fa', 'w')
		BEDfile = open('chr1.bed', 'w')  
		#name of chr
		start = 0
		end = 0
		N = 'N'*100
		concatFA.write('>chr1'+'\n')
		for record in record_dict:
			concatFA.write(str(record_dict[record].seq)+str(N)) # adds lots of N's to the end of the sequence
			seqlen = len(record_dict[record].seq)
			end = start + seqlen
			seqname = record_dict[record].id
			BEDfile.write('chr1'+'\t'+str(start)+'\t'+str(end)+'\t'+str(seqname)+'\t 0 \t + \t'+str(start)+'\t'+str(end)+'\n')
			start = end + 100          
		concatFA.close()
		BEDfile.close()
		   
		# move multo 
		#move concat and bedfile into new folders
		subprocess.call(['mv chr1.fa '+fullmultopath+'/files/tbb/tb'+header+'/fastaFiles/genomeFasta/noRandomChrom'], shell=True)
		subprocess.call(['mv chr1.bed '+fullmultopath+'/files/tbb/tb'+header+'/fastaFiles/annotationFiles/'], shell=True)
		#run MULTo 
		subprocess.call(['python '+fullmultopath+'/src/MULTo1.0.py -s tbb -a tb'+header+' -v '+numMismatch+' -O -p '+numCPU], shell=True)
		rmulto = header

	os.chdir(currDir) # sets working directory back to previous folder

	# Step 4

	# file = FASTQ file
	# bowtie aligns short reads to genome
	# we are taking our quality trimmed and adaptercut sequence file (from step before running trinity)
	# and we take all these trimmed sequence fragments and align them to the de novo assembled genome
	trim = arguments.start
	for file in filebasenames:
		stderr_bw = ""
		stderr_rp = ""
		if arguments.stderr == 0:
			stderr_bw = " 2> " + header + "/StandardError/bowtie-"+file+".txt"
			stderr_rp = " > " + header + "/StandardError/rpkmforgenes-"+file+".txt"
		if trim == 1:
			filepath = str(str(file) +'.fq')
		if trim == 0:
			filepath = str(header + "/"+str(file) +'_trimmed2.fq') # replace filepath = str(header + "/"+str(file) +'_trimmed2.fq' with str(filepath)
		subprocess.call(['bowtie -v '+numMismatch+' -m 1 -p '+numCPU+' -S -a '+fullmultopath+'/files/tbb/tb'+rmulto+'/bowtie_indexes/tb'+rmulto+'_genome/tb'+rmulto+'_no_random '+str(filepath)+" " +header + "/" + file+'_align.sam'+stderr_bw], shell=True)
		subprocess.call(['python '+fullmultopath+'/src/rpkmforgenes.py -i '+header + "/"+file+'_align.sam -samu -bedann -a '+fullmultopath+'/files/tbb/tb'+rmulto+'/fastaFiles/annotationFiles/chr1.bed -u '+fullmultopath+'/files/tbb/tb'+rmulto+'/MULfiles/tb'+rmulto+'_20-255/MULTo_files -o '+header + "/"+ file+'_MULTo.txt' + stderr_rp], shell=True)
	
def analyzeMulto(header, filebasenames, arguments, header_list, sampledict, scoredict ):	#create scoredict during blast step!also ID top match database?
	outfile = open(str(header+'/'+header+'_MULTo_analyzed_ALLdata.csv'), 'w')
	score_header = 'hit VSG\tpct_id_query\tpct_id_match\tmatch_VSG_length\n'	
	outfile.write(str('\t'.join(header_list)+'\tVSG\tPercent\tRPKM\t'+score_header+'\n'))
	
	for file in filebasenames:
		filepath = open(str(header + "/"+ file+'_MULTo.txt'), 'rU')
		filepath = filepath.read()
		multofilesplit = filepath.split('\n')
		FPKM = float(0)
		
		file_data = '\t'.join(sampledict[file])
		print file_data
		##calculate total FPKM
		for line in multofilesplit:
			if not line.startswith("#"):
					l = line.split("\t")
					if len(l) > 2:
						FPKM += float(l[2])
		
		###calculate final FPKM - this doesn't filter out those below 0.01% for calculation of percentage.. should it?
		for line in multofilesplit:
			if not line.startswith("#"):	
				l = line.split("\t")
				if len(l) == 3 :	
					percent = (float(l[2])/FPKM)*100
					VSG = str(l[0])
					VSG_data = scoredict[VSG]
					VSG_FPKM = float(l[2])
					outfile.write(file_data+'\t'+VSG+'\t'+str(percent)+'\t'+str(VSG_FPKM)+VSG_data)
		






	













