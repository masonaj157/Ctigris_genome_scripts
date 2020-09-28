#!/usr/bin/env python


import argparse
import csv
from Bio import SeqIO
import time

parser = argparse.ArgumentParser(description='writes a table of repeat content in a genome (although technically it will work with any fasta) given a gff file containing repeats.')
parser.add_argument("-s","--seqfile",
					type=str,
					default='',
					help="Sequence file in fasta format")
parser.add_argument("-g","--gff",
					type=str,
					default='',
					help="gff3 file containing repeats as a type")
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




def Remove_encompassed_repeats(reduced_gff):
	if len(reduced_gff) < 2 :
		return(reduced_gff)
	else:
		#start_time = time.time()
		reduced_gff.sort(key=lambda x: (x.end),reverse=True)
		reduced_gff.sort(key=lambda x: (x.start))
		reduced_gff = reduced_gff + ['end']
		entry_to_keep = []
		focal_entry = reduced_gff[0]
		if len(reduced_gff) > 1:
			test_entry = reduced_gff[1]
			#print(focal_entry.start, focal_entry.end)
			while test_entry != 'end':
				#print(test_entry.start, test_entry.end)
				if test_entry.end <= focal_entry.end :
					index = reduced_gff.index(test_entry)
					test_entry = reduced_gff[(index + 1)]
				else:
					entry_to_keep.append(focal_entry)
					focal_entry = test_entry 
					#print(reduced_gff.index(focal_entry))
					index = reduced_gff.index(test_entry)
					test_entry = reduced_gff[(index + 1)]
		else:
			entry_to_keep.append(focal_entry)
		#end_time = time.time()
		#print((end_time - start_time))
		return(entry_to_keep)
			


def Remove_overlapping_repeats(reduced_gff):
	if len(reduced_gff) < 2 :
		return(reduced_gff)
	else:
		#start_time = time.time()
		#if len(reduced_gff) == 1:
		#	print(reduced_gff)
		reduced_gff.sort(key=lambda x: (x.end),reverse=True)
		reduced_gff.sort(key=lambda x: (x.start))
		reduced_gff = reduced_gff + ['end']
		entry_to_keep = []
		focal_entry = reduced_gff[0]
		if len(reduced_gff) > 1:
			i = 1
			#entry_to_keep.append(focal_entry)
			while reduced_gff[i] != 'end':
				if reduced_gff[i].start < focal_entry.end:
					gff_entry = GFF([reduced_gff[i].seqid, reduced_gff[i].source, reduced_gff[i].type, focal_entry.start, reduced_gff[i].end, '','','',''])
					focal_entry = gff_entry
					i = i + 1
				else:
					entry_to_keep.append(focal_entry)
					focal_entry = reduced_gff[i]
					i = i + 1
		entry_to_keep.append(focal_entry)
		#end_time = time.time()
		#print((end_time - start_time))	
		return(entry_to_keep)



def Calc_repeat_content(sequence_list, gff_list, windowsize, windowstep) :
	start_time = time.time()
	content_list = []
	for seq in sequence_list:
		print(seq.id)
		reduced_gff = [repeat for repeat in gff_list if repeat.seqid == seq.id]
		endtime1 = time.time()
		print('time to reduce list = ' + str(endtime1 - start_time))
		#print(len(reduced_gff))
		reduced_gff = Remove_encompassed_repeats(reduced_gff)
		endtime2 = time.time()
		print('time to remove encompassing = ' + str(endtime2 - endtime1))
		#print(len(reduced_gff))
		#print(reduced_gff[-1])
		reduced_gff = Remove_overlapping_repeats(reduced_gff)
		endtime3 = time.time()
		print('time to remove overlap = ' + str(endtime3 - endtime2))
		if len(seq.seq) < windowsize:
			sum_repeat_length = 0
			for repeat in reduced_gff :
				repeat_length = abs(repeat.end - repeat.start)
				sum_repeat_length = sum_repeat_length + repeat_length
			repeat_content = sum_repeat_length / len(seq.seq)
			start = 1
			end = len(seq.seq)
			entry = [seq.id, start, end, repeat_content]
			print(repeat_content, start, end, seq.id)
			content_list.append(entry)
		else:
			start = 1
			end = windowsize
			while end <= len(seq.seq) :
				sum_repeat_length = 0
				for repeat in reduced_gff:
					if repeat.start >= start and repeat.end <= end :
						repeat_length = abs((repeat.end + 1) - repeat.start)
						sum_repeat_length = sum_repeat_length + repeat_length
					elif repeat.start < start and repeat.end > start and repeat.end <= end :
						repeat_length = abs((repeat.end + 1) - start)
						sum_repeat_length = sum_repeat_length + repeat_length
					elif repeat.start >= start and repeat.start < end and repeat.end > end :
						repeat_length = abs((end + 1) - repeat.start)
						sum_repeat_length = sum_repeat_length + repeat_length
				repeat_content = sum_repeat_length / windowsize
				entry = [seq.id, start, end, repeat_content]
				print(repeat_content,start,end, seq.id)
				content_list.append(entry)
				start = start + windowstep
				end = end + windowstep
		print()
	end_time = time.time()
	print((end_time - start_time))
	return(content_list)			
				
			
#########################################################################################

gff_list = GFF_parse(gff)


gff_list = [ entry for entry in gff_list if entry.source == 'repeat_gff:repeatmasker' and entry.type == 'match']

		
sequences = list(SeqIO.parse(seqfile,"fasta"))

repeat_list = Calc_repeat_content(sequences, gff_list, windowsize, windowstep)
#test = Calc_repeat_content(reduced_sequences, gff_list, windowsize, windowstep)
	
with open(output,'w') as outFile:
	outfile_writer = csv.writer(outFile, delimiter = '\t')
	outfile_writer.writerow(['scaffold', 'start', 'end', 'repeat_content'])
	for entry in repeat_list:
		outfile_writer.writerow(entry)


outFile.close()
				
			
#########################################################################################
