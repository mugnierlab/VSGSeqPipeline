# VSG_Pipeline

Organization Structure

To start:

Current Directory{
	Pipeline.py
	VSGFunctions.py
	[SequencingFile1].fastq
	[SequencingFile2].fastq
}

When Finished:

Current Directory{
	Pipeline.py
	VSGFunctions.py
	[SequencingFile1].fastq
	[SequencingFile2].fastq
	[Y-M-D-H_M]-[OptionalDescriptiveHeaderNames]{ # "Header"
		[SequencingFile1]{
				[SequencingFile1]_trimmed2.fq
			}
		[SequencingFile2]{
				[SequencingFile2]_trimmed2.fq
			}

		[Header]_[SequencingFile1]_contig.fa
		[Header]_[SequencingFile1]_orf_clean.fa
		[Header]_[SequencingFile1]_orf_trans_clean.fa
		[Header]_[SequencingFile1]_orf_trans.fa
		[Header]_[SequencingFile1]_orf_VSGs_scores.txt
		[Header]_[SequencingFile1]_orf_VSGs.fa
		[Header]_[SequencingFile1]_orf.fa
		[Header]_[SequencingFile2]_contig.fa
		[Header]_[SequencingFile2]_orf_clean.fa
		[Header]_[SequencingFile2]_orf_trans_clean.fa
		[Header]_[SequencingFile2]_orf_trans.fa
		[Header]_[SequencingFile2]_orf_VSGs_scores.txt
		[Header]_[SequencingFile2]_orf_VSGs.fa
		[Header]_[SequencingFile2]_orf.fa
		[Header]_orf_VSGs_merged.fa
		[Header]_orf_VSGs.fa 

		[Header]_[SequencingFile1]_orf.xml
		[Header]_[SequencingFile2]_orf.xml

		[SequencingFile1]_MULTo.txt
		[SequencingFile1]_Trinity.Trinity.fasta
		[SequencingFile2]_MULTo.txt
		[SequencingFile2]_Trinity.Trinity.fasta
		[Header]_orf_VSGs_merged.fa.clstr
		
		[SequencingFile1]_align.sam
		[SequencingFile2]_align.sam
	}
	
}