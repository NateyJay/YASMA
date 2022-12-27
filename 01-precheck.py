#!/usr/bin/env python3 

import sys
import os
from pprint import pprint

from subprocess import PIPE, Popen
from pathlib import Path

from os.path import isfile

from collections import Counter
from time import time

# import pyBigWig



import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-b', '--bam_file', 
	nargs="?",
	required=True,
	help='bamfile of aligned small RNAs (tested with shortstack)')

# parser.add_argument('-r', '--replicate_groups', 
# 	nargs='+', 
# 	required=True,
# 	help='list of readgroup names, following the order of their invocation in the bam header')


parser.add_argument('-o', '--output_directory', 
	nargs="?",
	default=f'SmoothLoc_{round(time())}',
	help='working folder for locus analysis')


# parser.add_argument('-s', '--sizes', 
	# nargs='+', 
	# help='list of sRNA sizes to be analyzed separately', 
	# default = [20,21,22,23,24])


parser.add_argument('-f', '--force', 
	action='store_true',
	default=False,
	help='force remake of supporting files')



args = parser.parse_args()
input_bam  = args.bam_file
# replicate_groups      = args.replicate_groups
out_dir    = args.output_directory
force = args.force
# rgs        = args.readgroups
# sizes      = args.sizes

Path(out_dir).mkdir(parents=True, exist_ok=True)


# input_bam = "/Users/jax/Desktop/+Colleagues/2022 - Consuelo - TxB/03-D-medfilt_loci/01out-BC.bocin.bam"
# input_bam = '/Volumes/HeavyPoint/+SeqLibraries/Control_libraries/Artha.3/merged_alignments.bam'
# input_bam = '/Volumes/HeavyPoint/+SeqLibraries/Control_libraries/Artha.2/merged_alignments.bam'
# out_dir   = "01out-wigs"



def samtools_stat(bam):

	stat_file = f"{out_dir}/01out-stats.txt"

	if not isfile(stat_file) and not force:


		call = ['samtools','stat', bam]

		p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')
		out,err=p.communicate()
		# print(out)
		print(err)

		with open(stat_file, 'w') as outf:
			outf.write(out)

	return(stat_file)





def process_stats(stat_file):
	length_c = Counter()

	read_count =0
	with open(stat_file, 'r') as f:
		for line in f:
			line = line.strip().split("\t")


			if line[0] == "RL":
				# print(line)
				read_count += int(line[2])
				length_c[int(line[1])] = int(line[2])

	print(f"{read_count:,} aligned reads")
	print()
	print("lengths:")
	for r in range(15, 30+1):
		d = length_c[r]
		p = round(d / read_count, 3)
		print(f"  {r}\t{p}\t{d:,}")

stat_file = samtools_stat(input_bam)

process_stats(stat_file)

















