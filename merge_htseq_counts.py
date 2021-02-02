#!/usr/bin/env python
# coding=utf-8
# A python script that merges single htsqe -count files into one.
# Usage: python merge_htseq_counts.py -i <directory> -o <output_file>

import argparse
import glob
import os

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input directory taht contains htseq *.counts files", required=True)
parser.add_argument("-o", "--output", help="Output file name", required=True)
args = parser.parse_args()

all_counts = {}
header = ['gene_id']

file_paths = glob.glob(args.input + "/*.counts")

for record in sorted(file_paths):

	full_file_name = os.path.split(record)[1]
	sample_id = full_file_name.split('.')[0]
	header.append(sample_id)

	with open(record, 'r') as input_file:

		for line in input_file:

			if line.startswith('__'):
				next
			else: 
				line = line.strip().split('\t')
				gene_id, count = line[0], line[1]

				if gene_id not in all_counts:
					all_counts[gene_id] = []

				all_counts[gene_id].append(count)

with open(args.output, 'w') as output_file:

	output_file.write('\t'.join(header) + '\n')

	for gene_id in sorted(all_counts.keys()):
		all_counts[gene_id].insert(0, gene_id)
		output_file.write('\t'.join(all_counts[gene_id]) + '\n')
