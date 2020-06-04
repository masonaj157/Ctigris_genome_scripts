#!/usr/bin/env python


import argparse
import csv
from Bio import SeqIO

parser = argparse.ArgumentParser(description='writes a table of CG content in a genome (although technically it will work with any fasta)')
parser.add_argument("-b","--bed",
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

bed = args.bed
windowsize = args.windowsize
windowstep = args.windowstep
output = args.output

#########################################################################################

class BED:
	def __init__(self, bed_entry):
		self.seqid = bed_entry[0]
		self.position = int(bed_entry[1])
		self.coverage = int(bed_entry[2])
		


def BED_parse(bed) :
	bed_list = []
	with open(bed) as OF:
		reader = csv.reader(OF, delimiter='\t')
		for row in reader :
			bed_entry = BED(row)
			bed_list.append(bed_entry)
	return(bed_list)


def Calc_window_coverage(bed_list, windowsize, windowstep):
	cov_list = []
	scaffolds = list(set([x.seqid for x in bed_list]))
	for scaffold in scaffolds:
		print(scaffold)
		reduced_bed = [x for x in bed_list if x.seqid == scaffold]
		if len(reduced_bed) >= windowsize :
			start = 0
			end = windowsize
			while end < len(reduced_bed):
				test_region = reduced_bed[start:end]
				test_region = [x.coverage for x in test_region]
				sum_cov = sum(test_region)
				ave_cov = (sum_cov / len(test_region))
				entry = [scaffold,start,(end-1),ave_cov]
				cov_list.append(entry)
				start = start + windowstep
				end = end + windowstep
		elif len(reduced_bed) < windowsize :
			test_region = reduced_bed
			test_region = [x.coverage for x in test_region]
			sum_cov = sum(test_region)
			ave_cov = (sum_cov / len(test_region))
			entry = [scaffold,0,(len(test_region)-1),ave_cov]
			cov_list.append(entry)
	return(cov_list)



#########################################################################################

bed_list = BED_parse(bed)


window_coverage = Calc_window_coverage(bed_list, windowsize, windowstep)

	
with open(output,'w') as outFile:
	outfile_writer = csv.writer(outFile, delimiter = '\t')
	outfile_writer.writerow(['scaffold', 'start', 'end', 'coverage'])
	for entry in window_coverage:
		outfile_writer.writerow(entry)


outFile.close()









