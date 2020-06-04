#!/usr/bin/env python


import argparse
import csv
from Bio import SeqIO

parser = argparse.ArgumentParser(description='writes a table of CG content in a genome (although technically it will work with any fasta)')
parser.add_argument("-s","--seqfile",
					type=str,
					default='',
					help="Sequence file in fasta format")
parser.add_argument("-g","--gff",
					type=str,
					default='',
					help="gff3 file containing genes as a type")
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

seqfile = args.seqfile
gff = args.gff
windowsize = args.windowsize
windowstep = args.windowstep
output = args.output

#########################################################################################

class GFF:
	def __init__(self, gff_entry):
		self.seqid = gff_entry[0]
		self.source = gff_entry[1]
		self.type = gff_entry[2]
		self.start = int(gff_entry[3])
		self.end = int(gff_entry[4])
		self.score = gff_entry[5]
		self.strand = gff_entry[6]
		self.phase = gff_entry[7]
		self.attributes = gff_entry[8]



def GFF_parse(gff) :
	gff_list = []
	with open(gff) as OF:
		reader = csv.reader(OF, delimiter='\t')
		for row in reader :
			if (row[0][0] != '#') and (len(row) > 1) :
				gff_entry = GFF(row)
				gff_list.append(gff_entry)
	return(gff_list)


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

gff_list = GFF_parse(gff)


gff_list = [ entry for entry in gff_list if entry.type == 'gene']

		
sequences = list(SeqIO.parse(seqfile,"fasta"))

gene_density = Calc_gene_density(sequences, gff_list, windowsize, windowstep)

	
with open(output,'w') as outFile:
	outfile_writer = csv.writer(outFile, delimiter = '\t')
	outfile_writer.writerow(['scaffold', 'start', 'end', 'genes'])
	for entry in gene_density:
		outfile_writer.writerow(entry)


outFile.close()









