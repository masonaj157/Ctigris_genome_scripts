#!/usr/bin/env python


import argparse
import csv
from Bio import SeqIO

parser = argparse.ArgumentParser(description='writes a table of CG content in a genome (although technically it will work with any fasta)')
parser.add_argument("-s","--seqfile",
					type=str,
					default='',
					help="Sequence file in fasta format")
parser.add_argument("-b","--bedgraph",
					type=str,
					default='',
					help="bedfile with coverage for each position")
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
					help="name of output file")
args=parser.parse_args()

seqfile = args.seqfile
bedgraph = args.bedgraph
windowsize = args.windowsize
windowstep = args.windowstep
output = args.output

#########################################################################################

class BED_GRAPH:
	def __init__(self, bed_entry):
		self.seqid = bed_entry[0]
		self.start = int(bed_entry[1])
		self.end = int(bed_entry[2])
		self.coverage = float(float(bed_entry[3]) / 100)


def BED_parse(bed) :
	bed_list = []
	with open(bed) as OF:
		reader = csv.reader(OF, delimiter='\t')
		for row in reader :
			bed_entry = BED_GRAPH(row)
			bed_list.append(bed_entry)
	return(bed_list)



def window_filter(bed_list, start, end):
	#print('start and end are', start, end)
	test_region = []
	for x in bed_list:
		#print(x.start,x.end)
		if x.start >= start and x.end <= end :
			#print('x is within window')
			test_region.append(x)
		elif x.start <= start and x.end >= end:
			#print('x is larger than window')
			x.start = start
			x.end =end
			test_region.append(x)
		elif x.start <= start and x.end > start and x.end <= end :
			#print('x overlaps start')
			x.start = start
			test_region.append(x)
		elif x.start >= start and x.start < end and x.end >= end :
			#print('x overlaps end')
			x.end = end
			test_region.append(x)
		elif x.start < start and x.end <= start :
			#print('x is lower than window')
			continue
		elif x.start >= end and x.end > end :
			#print('x is lower than window')
			continue
	return(test_region)



def Calc_window_coverage(bed_list, sequences, windowsize, windowstep):
	cov_list = []
	for seq in sequences:
		print(seq.id)
		reduced_bed = [x for x in bed_list if x.seqid == seq.id]		
		if len(seq.seq) >= windowsize :
			start = 0
			end = (windowsize - 1)
			while end <= (len(seq.seq) - 1):
				test_region = []
				test_region = window_filter(reduced_bed, start, end)
				if len(test_region) == 0:
					entry = [seq.id, start, end, 'NA']
					cov_list.append(entry)
				else:			
					values = []
					num_calls = []
					for entry in test_region:
						section_len = (entry.end - entry.start)
						num_calls.append(section_len)
						exp_methylated_C = float(section_len * entry.coverage)
						values.append(exp_methylated_C)
					sum_C = sum(values)
					sum_calls = sum(num_calls)
					ave_cov = float(sum_C / sum_calls)
					entry = [seq.id, start, end, ave_cov]
					print(entry)
					print(sum_C, sum_calls)
					cov_list.append(entry)
					start = start + windowstep
					end = end + windowstep
		elif len(seq.seq) < windowsize:
			test_region = [x for x in reduced_bed if x.seqid == seq.id]
			if len(test_region) == 0:
				entry = entry = [seq.id, 0, (len(seq.seq)-1), 'NA']
				cov_list.append(entry)
			else:
				values = []
				num_calls = []
				for entry in test_region:
					section_len = (entry.end - entry.start)
					num_calls.append(section_len)
					exp_methylated_C = float(section_len * entry.coverage)
					values.append(exp_methylated_C)
				sum_C = sum(values)
				sum_calls = sum(num_calls)
				ave_cov = float(sum_C / sum_calls)
				entry = [seq.id, 0, (len(seq.seq)-1), ave_cov]
				print(entry)
				print(sum_C, sum_calls)
				cov_list.append(entry)
	return(cov_list)
						



#########################################################################################

bed_list = BED_parse(bedgraph)

sequences = list(SeqIO.parse(seqfile,"fasta"))

window_coverage = Calc_window_coverage(bed_list, sequences, windowsize, windowstep)
#window_coverage = Calc_window_coverage(test_bed_list, test_seq_list, windowsize, windowstep)
	
with open(output,'w') as outFile:
	outfile_writer = csv.writer(outFile, delimiter = '\t')
	outfile_writer.writerow(['scaffold', 'start', 'end', 'coverage'])
	for entry in window_coverage:
		outfile_writer.writerow(entry)


outFile.close()