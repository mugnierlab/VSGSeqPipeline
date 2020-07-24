from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
from Bio import SeqIO 
from Bio import Seq
import subprocess
import argparse
import time
import os

def makeFilesList(filesList): #parses file containing information on files/samples
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
			sampledict['.'.join(linesplit[0].split('.')[:-1])] = line.split('\t')[0:] 
			filebasenames.append('.'.join(linesplit[0].split('.')[:-1]))
	return (filebasenames, header_list, sampledict)

def addSeqRecord(recordRD, start, end, count, ORF_outfile, trans_out_file, annotation): 
	sequence = recordRD.seq[start:end]
	ORF_outfile.write('>'+str(recordRD.id)+'_'+str(count)+'_'+str(start)+'_'+str(end)+"_" + annotation + '\n'+str(sequence)+'\n')
	trans_out_file.write('>'+str(recordRD.id)+'_'+str(count)+'_'+str(start)+'_'+str(end)+"_" + annotation + '\n'+str(sequence.translate())+'\n')

def addSeqRecord_RC(recordRD, start, end, count, ORF_outfile, trans_out_file, annotation):
	sequence = recordRD.seq[start:end].reverse_complement()
	ORF_outfile.write('>'+str(recordRD.id)+'_'+str(count)+'_'+str(start)+'_'+str(end)+'_RC_' + annotation +'\n'+str(sequence)+'\n')
	trans_out_file.write('>'+str(recordRD.id)+'_'+str(count)+'_'+str(start)+'_'+str(end)+"_RC_" + annotation + '\n'+str(sequence.translate())+'\n')

def findORFs(header, file, min_pro_len): #finds open reading frames in each assembled contig
		record_dict = SeqIO.index(header+ "/" + file+"_Trinity.Trinity.fasta","fasta")
		ORF_outfile = open(os.path.join(header, header+"_"+file+"_orf.fa"), 'w') # each file will get their own ORF output file
		contig_outfile = open(os.path.join(header, header+"_"+file+"_contig.fa"), 'w') # same as above, but for contigs the ORFs came from
		trans_out_file = open(os.path.join(header, header+"_"+file+'_orf_trans.fa'), 'w') # translated orf file
		for record in record_dict: # iterates through the sequences		
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
										if strand ==1:
											addSeqRecord(record_dict[record], frame+(trans_start*3), frame+(trans_end*3)+3, count, ORF_outfile, trans_out_file, file)
										else:
											addSeqRecord_RC(record_dict[record], max(frame, seq_len-frame-(trans_end*3)-3), seq_len-frame-(trans_start*3), count, ORF_outfile, trans_out_file, file)
										count += 1
									elif trans_end > min_pro_len : # protein isn't long enough, so just take the entire front end if thats long enough
										if strand ==1:
											addSeqRecord(record_dict[record], frame, frame+(trans_end*3)+3, count, ORF_outfile, trans_out_file, file)
										else:
											addSeqRecord_RC(record_dict[record], max(frame, seq_len-frame-(trans_end*3)-3), seq_len-frame, count, ORF_outfile, trans_out_file, file)
										count += 1
									# if neither protein is long enough,or even if it was, we will move on to find where the next start codon is, and find the end codon for that one
									trans_start = trans.find("M", trans_end)
									trans_end = trans.find("*", trans_start)		
								else: # there is a start codon AFTER the first appearing stop codon
									if trans_end > min_pro_len: # is the seq_start to end codon is long enough for the a min length protein
										if strand ==1:
											addSeqRecord(record_dict[record], frame, frame+(trans_end*3)+3, count, ORF_outfile, trans_out_file, file)
										else:
											addSeqRecord_RC(record_dict[record], max(frame, seq_len-frame-(trans_end*3)-3), seq_len, count, ORF_outfile, trans_out_file, file)
										count += 1
									trans_end = trans.find("*", trans_start) # will start checking the orf at that first found start codon, sets new trans_end to stop codon after the start	
							else: # start codon exists, but there is no stop codon anywhere! :O
							# this portion of the code will take the entire thing as an ORF, we are being generous
								if wholeAdded == False: # won't add if the entire strand was already added
									wholeAdded = True
									if strand ==1:
										addSeqRecord(record_dict[record], frame+(trans_start*3), seq_len, count, ORF_outfile, trans_out_file, file)
									else:
										addSeqRecord_RC(record_dict[record], frame, seq_len-frame-(trans_start*3), count, ORF_outfile, trans_out_file, file)
									count += 1
								trans_start = trans_len # won't go into loop, bc the entire rest of the sequence just got taken as an orf, so there is no point	in searching further, no more to search
						elif trans_end != -1: # start codon doesn't exist, but end codon does
							if trans_end > min_pro_len: # take entire sequence from start of seq to the found end codon as ORF if longer than min length
								if strand ==1:
									addSeqRecord(record_dict[record], frame, frame+(trans_end*3)+3, count, ORF_outfile, trans_out_file, file)
								else:
									addSeqRecord_RC(record_dict[record], max(frame, seq_len-frame-(trans_end*3)-3), seq_len-frame-(trans_end*3), count, ORF_outfile, trans_out_file, file)
								count += 1	
							trans_start = trans_len # won't search for more, since there aren't any start codons anywhere
						else: # there are no start or stop codons anywhere?!
							if wholeAdded == False:
								wholeAdded == True
								if strand ==1:
									addSeqRecord(record_dict[record], frame, seq_len, count, ORF_outfile, trans_out_file, file)
								else:
									addSeqRecord_RC(record_dict[record], frame, seq_len-frame, count, ORF_outfile, trans_out_file, file)
								trans_start = trans_len
								count += 1	
						
						trans_max = trans_len - min_pro_len # farthest a codon can be before it's below min protein length
						while trans_start < trans_max:
							if trans_start == -1: # if no more starts are found, doesn't matter if there is a stop or not, game over
								trans_start = trans_len
							elif trans_end == -1: #if start is found but no end is found
								if trans_len-trans_start > min_pro_len: #if ORF from found start to end sequence is long enough
									if strand ==1:
										addSeqRecord(record_dict[record], frame+(trans_start*3), seq_len, count, ORF_outfile, trans_out_file, file)
									else:
										addSeqRecord_RC(record_dict[record], frame, seq_len-frame-(trans_start*3), count, ORF_outfile, trans_out_file, file)
									count += 1
								trans_start = trans_len
							elif trans_end-trans_start > min_pro_len: # if we have both a start and an end and its long enough
								if strand ==1:
									addSeqRecord(record_dict[record], frame+trans_start*3, frame+(trans_end*3)+3, count, ORF_outfile, trans_out_file, file)
								else:
									addSeqRecord_RC(record_dict[record], max(frame, seq_len-frame-(trans_end*3)-3), seq_len-frame-(trans_start*3), count, ORF_outfile, trans_out_file, file)
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
	
def parseVSGblast(inputfile,filter): #parses blast of ORFs against VSG databases and identifies contigs containing VSGs
		inputfile = filename
		
		v = filename+'.xml'
		n = filename+'_nonVSG.xml' 
		s = filename+".fa"
		
		result_handle = open(v)
		
		blast_records = NCBIXML.parse(result_handle) # returns an iterator of the blast results
		record_dict = SeqIO.index(s,"fasta")
		
		outfile = open(s.split('.')[0]+'_VSGs.fa', 'w')
	
		hit_list = []
	
		if filter  == False: # blasted against nonVSG database
			nonVSGresult_handle = open(n)
			blast_records_nonVSG = NCBIXML.parse(nonVSGresult_handle)
			 # list of VSGs we have found! 
			exclude_list = []
					
			for blast_record_nonVSG in blast_records_nonVSG:
				for alignment in blast_record_nonVSG.alignments:
					for hsp in alignment.hsps: 
						percent_identity = (100.0 * hsp.identities) / alignment.length # hsp.identities is a tuple(bp matches, total bp in seq) to give percent match of sequence, percent identity is # of bp
						percent_query_identity = (100.0 * hsp.identities) / blast_record_nonVSG.query_letters
						if (percent_query_identity > 30 and hsp.identities > 300) or (percent_identity > 90):
							if not blast_record_nonVSG.query in exclude_list:
								exclude_list.append(str(blast_record_nonVSG.query))

			for blast_record in blast_records:
				for alignment in blast_record.alignments:
					for hsp in alignment.hsps:
						if hsp.expect < 1.0e-10: # hsp.expect = e value for the hsp value, the lower the e value, the more statistically significant 
							if not blast_record.query in hit_list: # if this query hasn't already been added to the hit list, add it now
								if not blast_record.query in exclude_list: # if the query isn't a fake VSG hit, add it now!
									hit_list.append(str(blast_record.query))
									SeqIO.write(record_dict[blast_record.query], outfile, "fasta")
									
		else: # didn't blast against nonVSG database
			for blast_record in blast_records:
				for alignment in blast_record.alignments:
					for hsp in alignment.hsps:
						if hsp.expect < 1.0e-10: # hsp.expect = e value for the hsp value, the lower the e value, the more statistically significant 
							if not blast_record.query in hit_list: # if this query hasn't already been added to the hit list, add it now
								hit_list.append(str(blast_record.query))				
								SeqIO.write(record_dict[blast_record.query], outfile, "fasta")
		
def makescoredict(xmlName): #analyzes output from blast against VSG databases, creates dictionary describing each assembled VSG
	
	#process blast
	result_handle = open(xmlName, 'r')
	blast_records = NCBIXML.parse(result_handle)
	hit_list = []
	scoredict = {}
	for blast_record in blast_records:
				for alignment in blast_record.alignments:
					for hsp in alignment.hsps:
						if hsp.expect < 1.0e-10: # hsp.expect = e value for the hsp value, the lower the e value, the more statistically significant 
							if not blast_record.query in hit_list: # if this query hasn't already been added to the hit list, add it now
								hit_list.append(str(blast_record.query))											# percent query aligned										# percent identity                     #hit VSG length                     #assembled VSG length 
								scoredict[str(blast_record.query)] = ('\t'+str(alignment.title)+'\t'+str((100.0 * hsp.identities) / blast_record.query_letters)+'\t'+str((100.0 * hsp.identities) / alignment.length)+'\t'+str(alignment.length)+'\t'+str(blast_record.query_letters)+'\n')
	
	return scoredict

def makeMulto(header, path, referencefasta, numMismatch, numCPU ): #Makes the MULTo database
	currDir = os.getcwd()
	fullmultopath = path
	os.chdir(fullmultopath)
	# make multo dir
	# MULTo identifies, stores and retrieves the minimum length required at each genomic position to be unique across the genome or transcriptome.
	#make multo file heirarchy from given basename
	os.makedirs(fullmultopath+'/files/tbb/tb'+header+'/fastaFiles/annotationFiles')
	os.makedirs(fullmultopath+'/files/tbb/tb'+header+'/fastaFiles/genomeFasta/noRandomChrom')

	# make bed
	record_dict = SeqIO.index(referencefasta,"fasta")
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
	#copy concat and bedfile into new foldersand move to current directory with name based on header
	subprocess.call(['cp chr1.fa '+fullmultopath+'/files/tbb/tb'+header+'/fastaFiles/genomeFasta/noRandomChrom'], shell=True)
	subprocess.call(['cp chr1.bed '+fullmultopath+'/files/tbb/tb'+header+'/fastaFiles/annotationFiles/'], shell=True)
	subprocess.call(['mv chr1.fa '+currDir+'/'+header+'/'+header+'.fa'], shell=True)
	subprocess.call(['mv chr1.bed '+currDir+'/'+header+'/'+header+'.bed'], shell=True)
	
	#run MULTo 
	subprocess.call(['python '+fullmultopath+'/src/MULTo1.0.py -s tbb -a tb'+header+' -v '+numMismatch+' -O -p '+numCPU], shell=True)

	os.chdir(currDir) # sets working directory back to previous folder

def analyzeMulto(header, filebasenames, header_list, sampledict, scoredict ):#parse the MULTo output files to make final results file
		
	outfile = open(str(header+'/'+header+'_RESULTS.txt'), 'w')
	if scoredict == False:
		outfile.write(str('\t'.join(header_list)+'\tVSG\tPercent\tRPKM\n'))
	else:
		outfile.write(str('\t'.join(header_list)+'\tVSG\tPercent\tRPKM\thit_VSG\tpct_id_vs_query\tpct_id_vs_match\tmatch_VSG_length\tassembled_VSG_length\n'))
	
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
		
		###calculate final percentage - this doesn't filter out those below 0.01% for calculation of percentage.. should it?
		for line in multofilesplit:
			if not line.startswith("#"):	
				l = line.split("\t")
				if len(l) == 3 :	
					percent = (float(l[2])/FPKM)*100
					VSG = str(l[0])
					VSG_FPKM = float(l[2])
					if scoredict == False:
						outfile.write(file_data+'\t'+VSG+'\t'+str(percent)+'\t'+str(VSG_FPKM)+'\n')
					else:
						VSG_data = scoredict[VSG]
						outfile.write(file_data+'\t'+VSG+'\t'+str(percent)+'\t'+str(VSG_FPKM)+VSG_data)



parser = argparse.ArgumentParser()

# arguments about input files and naming
parser.add_argument('-i',help='text file with .fastq names of sequence files to be put through the pipeline, same for entering any step of pipeline, fq', action="store", dest="i", default='')
parser.add_argument('-header', help='name for output folder, overrides default Y-m-d-H_M format ', action="store", dest='head', default='')

# trimming setting
parser.add_argument('-trim', help='perform quality trimming, default is to skip this step', action='store_false', dest='trim')
parser.add_argument('-g', help='trim_galore, stringency for trim galore', action ="store", dest = "g", default="3") 
parser.add_argument('-trimlen', help='cutadapt/trim_galore, minimum length of read for trimming', action ="store", dest = "trimlen", default="50") 

# trinity settings
parser.add_argument('-minp', help='Trinity, minimum protein length you are filtering for', action ="store", type=int, dest = "minp", default=300) 
parser.add_argument('-mem', help='Trinity, max memory allocation for trinity, G', action ="store", dest = "mem", default="10") 
parser.add_argument('-cpu', help='Trinity/cd-hit-est/MULTo/bowtie, number of processors to use', action="store", dest='cpu', default='2')

# Blast settings
parser.add_argument('-vsgdb', help='BLAST, name of the vsg database to blast against', action ="store", dest = "vsgdb", default="concatAnTattb427")
parser.add_argument('-NoNonVSGdb', help='BLAST, adding this flag means you dont want to BLAST ORFS against the non-VSG database' , dest='NoNonVSGdb', action='store_true')

# cd-hit-est parameters
parser.add_argument('-sit', help='cd-hit-est, sequence identiy threshold - how much the alignment has to match, percentage. value is 0.0 through 1.0 ', action ="store", dest = "sit", default=".98")
parser.add_argument('-t', help='cd-hit-est, number of threads for cd-hit-est to use, 0 default(all CPU will be used)', action ="store", dest = "t", default="0") #check this

# MULTO settings
parser.add_argument('-p', help='MULTo, path to MULTo1.0 folder. default is current directory, please dont use "~/", python doesnt like this in the path', action="store", dest='p', default='') # default assumes MULTo is in your home dir
parser.add_argument('-v', help='number of mismatches allowed for bowtie, default 2', action="store", dest='v', default='2')
parser.add_argument('-reuseMulto', help='MULTo, name of the multo files to be reused, typically same as header, default will make new MULTo files. Final output won\'t provide information on VSG', action="store", dest='rmulto', default='')

arguments = parser.parse_args()

sample_info =  makeFilesList(arguments.i)
filebasenames = sample_info[0]
header_list = sample_info[1]
sampledict = sample_info[2]
min_pro_len = arguments.minp
trim = arguments.trim
seqIdenThresh = arguments.sit
numCPU = arguments.cpu
max_memory_cdhit = 1000*int(arguments.mem)
path = arguments.p
numMismatch = arguments.v
rmulto = arguments.rmulto
trimlen = arguments.trimlen
stringency = arguments.g
currDir = os.getcwd()


#checking for all required software
if path != '':
	fullmultopath = path + '/MULTo1.0'
else:
	fullmultopath = currDir + '/MULTo1.0'

try:
	subprocess.check_call(["bowtie","--version"])
	subprocess.check_call(["blastn","-version"])
	subprocess.check_call(["trim_galore","--version"]) 
	subprocess.check_call(["cutadapt","--version"])
	subprocess.check_call(["python", str(fullmultopath+'/src/MULTo1.0.py'),"-h"])
	subprocess.check_call(["python", str(fullmultopath+'/src/rpkmforgenes.py'),"-h"])
except:
	raise SystemExit, "Problem opening required software. Is everything in your path?"
try:	
	subprocess.check_output(["cd-hit","--version"])
except subprocess.CalledProcessError as cpe:
		out = cpe.output
		if out[9:15] != 'CD-HIT':
			raise SystemExit, "Problem opening required software. Is cd-hit in your path?"	
except OSError:
	raise SystemExit, "Problem opening required software. Is cd-hit in your path?"
try:	
	subprocess.check_output(["Trinity","--version"])
except subprocess.CalledProcessError as cpe:
		out = cpe.output
		print out[0:7]
		if out[0:7] != 'Trinity':
			raise SystemExit, "Problem opening required software. Is Trinity in your path?"	
except OSError:
	raise SystemExit, "Problem opening required software. Is Trinity in your path?"	

header = time.strftime("%Y-%m-%d-%H_%M")  # Year/Month/Day-Hour:Minute , names the folder for output files

if header != "":
	header = str(arguments.head)
	
if not os.path.exists(header):
		os.makedirs(header) # creates the folder
		
if not os.path.exists(header+"/StandardError"):
	os.makedirs(header+"/StandardError")

scoredict = False

if trim == True:
	for file in filebasenames:
		if os.path.exists(str(file+'.fq')): 
			filepath = str(str(file) +'.fq')
		elif os.path.exists(str(file+'.fastq')):
			filepath = str(str(file) +'.fastq')
		stderr_tg = " 2> " + header + "/StandardError/trim_galore-"+file+".txt"
		stderr_ca = " 2> " + header + "/StandardError/cutadapt-"+file+".txt"
		subprocess.call(['trim_galore --stringency '+str(stringency)+' --length '+trimlen+' --dont_gzip --output_dir ' + header + "/ " +str(filepath)+stderr_tg], shell=True) # trim off sequenceing adapters
		subprocess.call(['cutadapt -m '+trimlen+' -b ATTTAGGTGACACTATAG -b CTATAGTGTCACCTAAAT '+ header + "/"+str(file)+'_trimmed.fq > '+ header + "/"+str(file)+'_trimmed2.fq'+stderr_ca], shell=True) # trim off SP6 sequences (from VSG PCR step)
		subprocess.call(['rm '+ header + "/"+str(file)+'_trimmed.fq'], shell=True) # removes intermediate trimmed file 

if rmulto == '': #if you need to make a multo database from your VSGs
	for file in filebasenames:
		if trim == True:
			filepath = str(header + "/"+str(file) +'_trimmed2.fq')
		else:
			if os.path.exists(str(file+'.fq')): 
				filepath = str(str(file) +'.fq')
			elif os.path.exists(str(file+'.fastq')):
				filepath = str(str(file) +'.fastq')
		print filepath
		stderr_tr = " > " + header + "/StandardError/Trinity-"+str(file)+".txt" 
		subprocess.call(['Trinity --seqType fq --max_memory '+arguments.mem+'G --full_cleanup --single ' + str(filepath) + ' --CPU '+numCPU+' --min_contig_length ' + str(min_pro_len*3) + ' --no_path_merging --no_normalize_reads --output ' + header+"/"+file+"_Trinity/" + stderr_tr], shell=True)
		
		findORFs(header, file, min_pro_len)
			
	
	for file in filebasenames:
		filename = header + "/"+header+"_"+file+"_orf"
		vsgdDbName = arguments.vsgdb
		subprocess.call(['blastn -db '+vsgdDbName+' -query '+filename+'.fa -outfmt 5 -out '+filename+'.xml'], shell=True)
		#blast nonVSG
		if arguments.NoNonVSGdb == False:
			subprocess.call(['blastn -db NOTvsgs -query '+filename+'.fa -outfmt 5 -out '+filename+'_nonVSG.xml'], shell=True)
		parseVSGblast(filename, arguments.NoNonVSGdb) 
	
	##concatenate all VSG ORFs, merge with cd-hit, blast and make dictionary of results
	for file in filebasenames:
		subprocess.call(['cat '+header + "/"+header+"_"+file+"_orf_VSGs.fa >> " +header + "/"+header+"_orf_VSGs.fa"], shell=True) # concatenates all the files for MULTo
	
	stderr_cd = " > " + header + "/StandardError/cdhitest.txt"
	subprocess.call(['cd-hit-est -i '+header+"/"+header+'_orf_VSGs.fa '+' -o '+header+"/"+header+'_orf_VSGs_merged.fa -d 0 -c ' + seqIdenThresh + ' -n 8 -G 1 -g 1 -s 0.0 -aL 0.0 -M '+str(max_memory_cdhit)+' -T ' + numCPU + stderr_cd], shell=True)
	
	fastaName = header+"/"+header+'_orf_VSGs_merged.fa'
	xmlName = str(header+"/"+header+'_orf_VSGs_merged.xml')
	subprocess.call(['blastn -db '+arguments.vsgdb+' -query '+fastaName+' -outfmt 5 -out '+xmlName], shell=True)
	scoredict = makescoredict(xmlName)

	referencefasta = currDir +'/' + fastaName
	
	makeMulto(header, fullmultopath, referencefasta, numMismatch, numCPU)

for file in filebasenames:
	stderr_bw = " 2> " + header + "/StandardError/bowtie-"+file+".txt"
	stderr_rp = " > " + header + "/StandardError/rpkmforgenes-"+file+".txt"
	if trim == False:
		filepath = str(str(file) +'.fq')
	if trim == True:
		filepath = str(header + "/"+str(file) +'_trimmed2.fq') 
	if rmulto != "":
		multoname = rmulto
	else:
		multoname = header
	subprocess.call(['bowtie -v '+numMismatch+' -m 1 -p '+numCPU+' -S -a '+fullmultopath+'/files/tbb/tb'+multoname+'/bowtie_indexes/tb'+multoname+'_genome/tb'+multoname+'_no_random '+str(filepath)+" " +header + "/" + file+'_align.sam'+stderr_bw], shell=True)
	subprocess.call(['python '+fullmultopath+'/src/rpkmforgenes.py -i '+header + "/"+file+'_align.sam -samu -bedann -a '+fullmultopath+'/files/tbb/tb'+multoname+'/fastaFiles/annotationFiles/chr1.bed -u '+fullmultopath+'/files/tbb/tb'+multoname+'/MULfiles/tb'+multoname+'_20-255/MULTo_files -o '+header + "/"+ file+'_MULTo.txt' + stderr_rp], shell=True)

analyzeMulto(header, filebasenames, header_list, sampledict, scoredict)
	

#tidy up all the output files

os.makedirs(header+"/bowtie_alignments")
subprocess.call(["mv "+header+"/*.sam "+header+"/bowtie_alignments/"], shell=True)
os.makedirs(header+"/MULTo_output")
subprocess.call(["mv "+header+"/*_MULTo.txt "+header+"/MULTo_output/"], shell=True)
if trim == True:
	os.makedirs(header+"/trimmed_files")
	subprocess.call(["mv "+header+"/*trimmed2.fq "+header+"/trimmed_files/"], shell=True)
	subprocess.call(["mv "+header+"/*trimming* "+header+"/trimmed_files/"], shell=True)
if rmulto == '':
	os.makedirs(header+"/Trinity_assemblies")
	subprocess.call(["mv "+header+"/*Trinity.fasta "+header+"/Trinity_assemblies/"], shell=True)
	os.makedirs(header+"/VSG_reference_fasta/")
	subprocess.call(["mv "+header+"/*merged* "+header+"/VSG_reference_fasta/"], shell=True)
	subprocess.call(["mv "+header+"/"+header+".fa "+header+"/VSG_reference_fasta/"], shell=True)
	subprocess.call(["mv "+header+"/"+header+".bed "+header+"/VSG_reference_fasta/"], shell=True)
	os.makedirs(header+"/intermediate_VSG_identification_files/")
	subprocess.call(["mv "+header+"/*_orf* "+header+"/intermediate_VSG_identification_files/"], shell=True)
	subprocess.call(["mv "+header+"/*_contig.fa "+header+"/intermediate_VSG_identification_files/"], shell=True)












