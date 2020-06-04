#!/usr/bin/env python


import argparse
import csv
from Bio import SeqIO

parser = argparse.ArgumentParser(description='writes a table of CG content in a genome (although technically it will work with any fasta)')
parser.add_argument("-i","--input",
					type=str,
					default='',
					help="Sequence file in fasta format")
parser.add_argument("-wsi","--windowsize",
					type=int,
					default=10000,
					help="window size in basepairs")
parser.add_argument("-wst","--windowstep",
					type=int,
					default=10000,
					help="window step length in basepairs")
parser.add_argument("-o","--output",
					type=str,
					default='',
					help="name of cytoband output file")
args=parser.parse_args()

seqfile = args.input
output = args.output
windowsize = args.windowsize
windowstep = args.windowstep

#########################################################################################

def Calc_CG(test_seq) :
	C = test_seq.count('C')
	G = test_seq.count('G')
	A = test_seq.count('A')
	T = test_seq.count('T')
	CG = C + G
	Total = C + G + A + T
	CG_perc = CG / Total
	return(CG_perc)
	
	
def Calc_CG_Window(sequences, window_size, step_size) :
	values = []
	for seq in sequences:
		if len(seq.seq) >= windowsize :
			start = 0
			end = window_size
			while ( end <= len(seq.seq) ) :
				test_seq = seq.seq[start:end]
				CG_perc = Calc_CG(test_seq)
				entry = [seq.id,start,end,CG_perc]
				values.append(entry)
				start = start + step_size
				end = end + step_size
		elif len(seq.seq) < windowsize :
			start = 0
			end = len(seq.seq)
			test_seq = seq.seq[start:end]
			CG_perc = Calc_CG(test_seq)
			entry = [seq.id,start,end,CG_perc]
			values.append(entry)
	return(values)



def Calc_gene_density(sequence_list, gff_list, windowsize, windowstep) :
	density_list = []
	for seq in sequence_list:
		if len(seq.seq) >= windowsize :
			start = 0
			end = windowsize
			while end < len(seq.seq) :
				reduced_gff = [gene for gene in gff_list if gene.seqid == seq.id and gene.start >= start and gene.start <= end]
				num_genes = len(reduced_gff)
				entry = [seq.id, start, end, num_genes]
				density_list.append(entry)
				start = start + windowstep
				end = end + windowstep
		elif len(seq.seq) < windowsize :
			reduced_gff = [gene for gene in gff_list if gene.seqid == seq.id]
			num_genes = len(reduced_gff)
			start = 0
			end = len(seq.seq)
			entry = [seq.id, start, end, num_genes]
			density_list.append(entry)
	return(density_list)
		
#########################################################################################		
		
sequences = list(SeqIO.parse(seqfile,"fasta"))


CG_entry = Calc_CG_Window(sequences, windowsize, windowstep)


with open(output,'w') as outFile:
	outfile_writer = csv.writer(outFile, delimiter = '\t')
	outfile_writer.writerow(['scaffold', 'start', 'end', 'CG'])
	for entry in CG_entry:
		outfile_writer.writerow(entry)


outFile.close()
	