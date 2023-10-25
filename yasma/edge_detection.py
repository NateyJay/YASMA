
import sys
import os

import click
from click_option_group import optgroup

from pathlib import Path
import shutil
from os.path import isfile, isdir
from collections import Counter#, deque
from pprint import pprint
from random import sample

import numpy as np
# from statistics import quantiles
import math

from .generics import *
from .cli import cli

from statistics import mean, median, StatisticsError

from time import time

import re







@cli.command(group='Utilities', help_priority=3)



@optgroup.group('\n  Basic options',
                help='')

@optgroup.option("-a", "--alignment_file", 
	required=False, 
	type=click.UNPROCESSED, callback=validate_path,
	help='Alignment file input (bam or cram).')

@optgroup.option('-r', '--annotation_readgroups', 
	required=False,
	multiple=True,
	help="List of read groups (RGs, libraries) to be considered for the annotation. 'ALL' uses all readgroups for annotation, but often pertainent RGs will need to be specified individually.")

@optgroup.option("-o", "--output_directory", 
	# default=f"Annotation_{round(time())}", 
	required=True,
	type=click.Path(),
	help="Directory name for annotation output.")


@optgroup.group('\n  Coverage options',
                help='')

@optgroup.option('--kernel_window',
	default=40,
	help='This is the bandwidth smoothing window for the kernel coverage method. Default 40 nt.')

@optgroup.option('--coverage_method',
	default='kernel',
	type=click.Choice(['kernel', 'depth']),
	help="Choice between two methods for finding genomic coverage. 'Kernel' uses a density smoothing based on the center of the read. 'depth' uses samtools depth and can be sensitive to read-length.")


@optgroup.option("--rpm_threshold",
	default=0.5,
	help="Depth threshold in reads per million for discovery of a locus peak. Default 0.5 RPM.")

@optgroup.option('--peak_threshold',
	default=0.05,
	help='Minimum depth, as a proportion of the maximum coverage depth in a locus, used for trimming low-depth edges from a locus. Default 0.05 proportion of peak.')


def edge_detection(**params):
	'''Annotator using peak identification and similarity merging.'''

	rc = requirementClass()
	rc.add_samtools()
	rc.check()

	ic = inputClass(params)
	ic.check(['alignment_file', 'annotation_readgroups'])

	output_directory        = str(ic.output_directory)
	alignment_file          = ic.inputs["alignment_file"]
	annotation_readgroups   = ic.inputs['annotation_readgroups']
	project_name            = ic.inputs['project_name']

	kernel_window           = params['kernel_window']
	coverage_method         = params['coverage_method']
	peak_threshold          = params['peak_threshold']
	rpm_threshold           = params['rpm_threshold']

	
	dir_name = 'edge'

	import scipy


	Path(output_directory+ f"/{dir_name}/").mkdir(parents=True, exist_ok=True)

	log_file = f"{output_directory}/{dir_name}/log.txt"
	sys.stdout = Logger(log_file)


	cl_i = 0

	# assert 0 < sensitivity <= 1.0, "sensitivity must be 0 < S <= 1.0"
	assert 0 < peak_threshold < 1, "peak_trim must be between 0 and 1."

	# half_window = int(window/2)

	chromosomes, bam_rgs = get_chromosomes(alignment_file)
	# print(alignment_file)
	annotation_readgroups = check_rgs(annotation_readgroups, bam_rgs)

	# chromosomes = [c for c in chromosomes if c[0] == 'NC_037320.1']

	chrom_depth_c = get_global_depth(output_directory, alignment_file, aggregate_by=['rg','chrom'])

	keys = list(chrom_depth_c.keys())
	for key in keys:
		if key[0] in annotation_readgroups:
			chrom_depth_c[key[1]] += chrom_depth_c[key]

		del chrom_depth_c[key]


	## preparing output files

	# gff_file = f"{output_directory}/{dir_name}/loci.gff3"

	# with open(gff_file, 'w') as outf:
	# 	print("##gff-version 3", file=outf)

	# 	for chrom, chrom_length in chromosomes:
	# 		print(f"##sequence-region   {chrom} 1 {chrom_length}", file=outf)


	# results_file = f"{output_directory}/{dir_name}/loci.txt"
	# with open(results_file, 'w') as outf:
	# 	print("\t".join(assessClass().header), file=outf)


	# reads_file = f"{output_directory}/{dir_name}/reads.txt"
	# with open(reads_file, 'w') as outf:
	# 	print(TOP_READS_HEADER, file=outf)



	# region_file = f"{output_directory}/{dir_name}/regions.gff3"
	# with open(region_file, 'w') as outf:
	# 	print("##gff-version 3", file=outf)

	# 	for chrom, chrom_length in chromosomes:
	# 		print(f"##sequence-region   {chrom} 1 {chrom_length}", file=outf)

	# merge_file = f"{output_directory}/{dir_name}/merges.txt"
	# with open(merge_file, 'w') as outf:
	# 	outf.write('')



	overall_d = {}

	all_loci = []
	total_read_count = 0
	cl_i = 0


	# total_reads = sum(chrom_depth_c.values())
	read_equivalent = 1 / sum(chrom_depth_c.values()) * 1000000

	# depth_wig = wiggleClass(f"{output_directory}/peak/Depths.wig", rpm_threshold, total_reads)
	# peak_wig = wiggleClass(f"{output_directory}/peak/Depths_peak.wig", rpm_threshold, total_reads)


	## iterating through chromosomes and reads

	# print("Annotating loci for each chromosome...")
	# print()

	# for chrom_count, chrom_and_length in enumerate(chromosomes):


	print(f" {len(chromosomes)} chromosomes/scaffolds/contigs")
	print(f" {sum([l for c,l in chromosomes]):,} bp total")
	print(f" {sum(chrom_depth_c.values()):,} reads")
	print()

	loci = []
	# chrom, chrom_length = chrom_and_length

	# print()
	# print(f"{chrom_count+1} / {len(chromosomes)}")
	# print(f"chrom: {chrom}")
	# print(f"       {chrom_length:,} bp")
	# print(f"       {chrom_depth_c[chrom]:,} reads")


	def get_peaks(bam, rgs):

		depth_c = Counter()

		perc = percentageClass(1, chrom_length)


		c1 = ['samtools', 'view', '-h', '-F', '4']
		for rg in rgs:
			c1 += ['-r', rg]
		c1 += [bam]#, chrom]

		p1 = Popen(c1, stdout=PIPE, stderr=PIPE, encoding=ENCODING)

		c2 = ['samtools', 'depth', '-']
		p2 = Popen(c2, stdin=p1.stdout, stdout=PIPE, stderr=PIPE, encoding=ENCODING)


		
		print(f"   reading position depths... 0%   ", end='\r', flush=True)
		for i,line in enumerate(p2.stdout):

			chrom, pos, depth = line.strip().split('\t')


			depth_c[(chrom, int(pos))] = int(depth)


			perc_out = perc.get_percent(i)
			if perc_out:
				print(f"   reading position depths... {perc_out}%", end='\r', flush=True)

		p2.wait()

		print()

		return(depth_c)


	out_file = Path(output_directory, "edge", "edges.txt")
	with open(out_file, 'w') as outf:
		print('chrom','position','density', 'last_avg', 'next_avg','delta','p_delta', 'f_delta', sep='\t', file=outf)



	def get_kernel_peaks(chrom, bam, rgs):

		half_window = math.floor(kernel_window/2)


		perc = percentageClass(1, chrom_depth_c[chrom])


		depth_c = Counter()


		reads = samtools_view(bam, rgs=rgs, locus = chrom)

		print(f"   reading position depths ..... 0%", end='\r', flush=True)

		i = 0
		for i, read in enumerate(reads):
			_, length, _, pos, chrom, _, _, _ = read
		# call = ['samtools', 'view', '-F', '4']
		# for rg in rgs:
		# 	call += ['-r', rg]
		# call += [bam, chrom]

		# p = Popen(call, stdout=PIPE, stderr=PIPE, encoding=ENCODING)

		# for i,line in enumerate(p.stdout):
		# 	line = line.strip().split("\t")


		# 	pos    = int(line[3])
		# 	length = int(line[5].rstrip("M"))


			pos += math.floor(length / 2)

			depth_c[pos-half_window] += 1

			perc_out = perc.get_percent(i)
			if perc_out:
				print(f"   reading position depths ..... {perc_out}%\t{i:,} reads", end='\r', flush=True)

		# p.wait()

		# print(f"   reading position depths ..... 100% ({i} reads)", end='\r', flush=True)

		# total_reads = i + 1


		# pprint(chrom_depth_c)
		# print(chrom)
		# print(length(depth_c[chrom]))
		# print(length(list(depth_c.items())))
		# print()
		# sys.exit()

		print()


		## Performing kernel density smoothing based on mean or median or sum?

		window_dq = deque()


		perc = percentageClass(1, sum([l for c,l in chromosomes if c == chrom]))
		ii =0

		# for chrom, chrom_length in chromosomes:
		kernel_c = Counter()


		for i in range(chrom_length):
			ii += 1

			window_dq.append(depth_c[i])

			summary_value = round(sum(window_dq),2)

			if summary_value > 0:
				kernel_c[i] = summary_value

			



			perc_out = perc.get_percent(ii)
			if perc_out:
				print(f"   calculating kernel density .. {perc_out}% \t{round(ii/1000000, 1)}M nt ", end='\r', flush=True)


			if len(window_dq) > kernel_window:
				window_dq.popleft()



		print()

		return(kernel_c)


	delta_bw = 100

	for chrom_count, chrom_and_length in enumerate(chromosomes):


		chrom, chrom_length = chrom_and_length
		print(f"{chrom} {chrom_length:,} nt)")
		peak_c = get_kernel_peaks(
			chrom,
			bam=alignment_file, 
			rgs=annotation_readgroups)


		last_set = deque()
		next_set = deque()

		for r in range(delta_bw):
			last_set.append(peak_c[0])
			next_set.append(peak_c[r])

		with open(out_file, 'a') as outf:

			perc = percentageClass(1, sum([l for c,l in chromosomes if c == chrom]))
			for pos in range(chrom_length):

				perc_out = perc.get_percent(pos)
				if perc_out:
					print(f"   calculating deltas ....... {perc_out}% \t{round(pos/1000000, 1)}M nt ", end='\r', flush=True)


				last_set.popleft()
				last_set.append(next_set.popleft())
				next_set.append(peak_c[pos])

				last_avg = mean(last_set)
				next_avg = mean(next_set)

				delta = next_avg - last_avg

				try:
					fold_delta = max(last_set + next_set) / min(last_set)
				except ZeroDivisionError:
					fold_delta = 0

				try:
					p_delta = delta / max(last_set + next_set)
				except ZeroDivisionError:
					p_delta = 0.0

				# print(last_set)
				# print(next_set)


				# print(chrom, pos, peak_c[pos], last_avg, next_avg, round(delta,4), round(p_delta,4), sep='\t')
				print(chrom, pos, peak_c[pos], last_avg, next_avg, round(delta,4), round(p_delta,4), round(fold_delta), sep='\t', file=outf)

				# input()






