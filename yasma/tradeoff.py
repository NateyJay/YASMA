

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
from itertools import chain

import numpy as np
# from statistics import quantiles
import math

from .generics import *
from .cli import cli

from statistics import mean, median, StatisticsError, stdev

from time import time

import re

from tqdm import tqdm


from datetime import datetime
from pympler import asizeof


class elapsedClass():
	def __init__(self):

		self.start_time = datetime.now()

	def __str__(self):
		current_time = datetime.now()
		elapsed = current_time - self.start_time

		return f" elapsed: {elapsed}"



class sizeClass():
	def __init__(self, 
		sizes=[],
		min_size=15,
		max_size=30):


		self.min_size=min_size
		self.max_size=max_size

		self.depth = 0

		self.size_c = Counter()

		if len(sizes) > 0:
			self.update(sizes)

	def get_keys(self, size):
		keys = set()

		keys.add((size,))

		keys.add((size-1, size+0,))
		keys.add((size-0, size+1,))

		keys.add((size-2, size-1, size+0,))
		keys.add((size-1, size+0, size+1,))
		keys.add((size-0, size+1, size+2,))

		keys = [k for k in keys if min(k) >= self.min_size or max(k) <= self.max_size]
		return(keys)



	def update(self, sizes):

		# if type(sizes) == int:
		# 	sizes = [sizes]

		for size in sizes:


			if 15 <= size <= 30:

				self.depth += 1
				self.size_c.update(self.get_keys(size))

				# for mer in [1,2,3]:
				# 	self.size_d[mer].update(self.size_key_d[mer][size])


	def get(self):

		mc = self.size_c.most_common()

		self.size_1_depth, self.size_2_depth, self.size_3_depth = 0,0,0
		self.size_1_key, self.size_2_key, self.size_3_key = None, None, None

		for key, depth in mc:
			if self.size_1_depth == 0 and len(key) == 1:
				self.size_1_key, self.size_1_depth = key, depth

			if self.size_2_depth == 0 and len(key) == 2:
				self.size_2_key, self.size_2_depth = key, depth

			if self.size_3_depth == 0 and len(key) == 3:
				self.size_3_key, self.size_3_depth = key, depth

			if self.size_1_depth * self.size_2_depth * self.size_3_depth > 0:
				break

		# pprint(self.size_c.most_common(10))
		# print(self.depth, "->", round(self.depth/2), "min")
		# print(self.size_1_key, self.size_1_depth, sep="\t")
		# print(self.size_2_key, self.size_2_depth, sep="\t")
		# print(self.size_3_key, self.size_3_depth, sep="\t")

		if self.size_1_depth > self.depth * 0.5:
			sizecall = self.size_1_key

		elif self.size_2_depth > self.depth * 0.5 and self.depth > 30:
			sizecall = self.size_2_key

		elif self.size_3_depth > self.depth * 0.5 and self.depth > 60:
			sizecall = self.size_3_key

		else:
			sizecall = tuple("N")

		self.sizecall = sizecall
		# print(self.depth)
		# print(self.sizecall)
		# sys.exit()

		return(sizecall)


	def __str__(self):
		out = self.get()
		# print(out)
		out = "_".join(map(str, out))

		return(out)


	def __eq__(self, other):

		self.get()
		other.get()

		# print(self.sizecall)
		# print(other.sizecall)


		scall = set(self.sizecall)
		ocall = set(other.sizecall)


		def expand_call(call):
			if len(call) == 3:
				call.add("N")
				return(call)

			if "N" in call:
				return(call)

			call.add(min(call)-1)
			call.add(max(call)+1)

			return(call)


		scall = expand_call(scall)
		ocall = expand_call(ocall)


		common = scall.intersection(ocall)

		if len(common) > 1:
			return True
		elif "N" in common:
			return True
		else:
			return False

	def __add__(self, other):
		self.size_c += other.size_c
		self.depth += other.depth
		return(self)




class assessClass():
	'''produces a line assessment of a locus, similar to ShortStack3'''

	def __init__(self):

		self.header = ['Locus','Name','Length','Reads','RPM']
		self.header += ['UniqueReads','FracTop','Strand','MajorRNA','MajorRNAReads','Complexity']
		self.header += ['Gap', 'size_1n','size_1n_depth', 'size_2n','size_2n_depth', 'size_3n','size_3n_depth', 'sizecall']



	def format(self, locus, seq_c, strand_c, sizecall, aligned_depth, last_stop):

		name, chrom, start, stop = locus


		### Basic information

		# input_len = len(reads)

		# reads = [r for r in reads if not r[3] + r[1] < start or not r[3] > stop ]
		depth = sum(seq_c.values())
		rpm = depth / aligned_depth * 1000000





		### ShortStack standard metrics

		# seq_c       = Counter()
		# strand_c    = Counter()
		# size_d      = {1 : Counter(), 2 : Counter(), 3 : Counter()}

		# sizecall = sizeClass()

		# for read in reads:
		# 	sam_strand, sam_length, sam_size, sam_pos, sam_chrom, sam_rg, sam_read, sam_name = read



		# 	if 15 <= sam_length <= 30:

		# 		# print(sam_length)
		# 		# print(self.get_size_keys(sam_length, 1))
		# 		# print(self.get_size_keys(sam_length, 2))
		# 		# print(self.get_size_keys(sam_length, 3))
		# 		# print(self.get_size_keys(sam_length, 4))
		# 		# sys.exit()



		# 		seq_c[sam_read] += 1
		# 		strand_c[sam_strand] += 1

		# 		sizecall.update([sam_length])

		# 		# size_d[1][sam_length] += 1
		# 		# size_d[2].update(size_key_d[2][sam_length])
		# 		# size_d[3].update(size_key_d[3][sam_length])

		# if sum(strand_c.values()) == 0:
		# 	# this solution is a stop-gap. Apparently, in some species low-lambdas can lead to trace-regions being annotated, which only have reads outside of accepted ranges. 
		# 	return(None,None)



		unique_reads = len(seq_c.keys())
		frac_top = strand_c["+"] / sum(strand_c.values())

	
		if frac_top > 0.8:
			strand = "+"
		elif frac_top < 0.2:
			strand = "-"
		else:
			strand = "."

		major_rna = seq_c.most_common()[0][0]
		major_rna_depth = seq_c.most_common()[0][1]

		complexity = unique_reads / depth



		### More derived metrics


		gap = start - last_stop


		



		frac_top   = round(frac_top,3)
		complexity = round(complexity,3)
		rpm        = round(rpm,3)

		sizecall.get()

		result_line = [f"{chrom}:{start}-{stop}", name, stop-start, depth, rpm]
		result_line += [unique_reads, frac_top, strand, major_rna, major_rna_depth, complexity]
		result_line += [
			gap, 
			sizecall.size_1_key, sizecall.size_1_depth,
			sizecall.size_2_key, sizecall.size_2_depth,
			sizecall.size_3_key, sizecall.size_3_depth,
			sizecall
		]


		if 'N' in sizecall.sizecall:
			feature_type = "OtherRNA"
		else:
			feature_type = f"RNA_{sizecall}"

		if start < 1:
			start = 1
		gff_line = [
			chrom, 'yasma_locus',feature_type, start, stop, '.', strand, '.',
			f'ID={name};sizecall={sizecall};depth={depth};rpm={rpm};fracTop={frac_top};majorRNA={major_rna}'
		]


		return(result_line, gff_line)



def get_bin_threshold(cdf_c, to_save=False, to_print=False):

	def cumsum(l):
		total = 0
		for x in l:
			total += x
			yield total

	### try 1
	depths     = [d      for d in range(max(cdf_c.keys())+1) if cdf_c[d] > 0]
	bin_counts = [cdf_c[d] for d in range(max(cdf_c.keys())+1) if cdf_c[d] > 0]
	p_genome = [1-(c / sum(bin_counts)) for c in cumsum(bin_counts)]
	total_reads = [depths[i] * bin_counts[i] for i in range(len(depths))]
	p_reads = [(sum(total_reads) - c) / sum(total_reads) for c in cumsum(total_reads)]
	averages = [((1-p_genome[i]) + p_reads[i])/2 for i in range(len(depths))]


	for highest_average,a in enumerate(averages):
		if a == max(averages):
			break

	if to_save:
		with open("comparison_table.txt", 'w') as outf:
			print("p_gen\tp_read\tavg\tbcnt\tdpth\ttotr", file=outf)
			for i in range(0, len(depths)):

				print(
					round(p_genome[i], 4), 
					round(p_reads[i], 4),
					round(averages[i], 4), 
					round(bin_counts[i], 4), 
					round(depths[i], 4), 
					round(total_reads[i], 4),
					sep='\t', file=outf)


	if to_print:
		print()
		print()
		print("\tp_gen\tp_read\tavg\tbcnt\tdpth\ttotr")

		for i in range(0, highest_average + 25):
			if i == highest_average:
				print(color.BOLD)
				print("",
					round(p_genome[i], 4), 
					round(p_reads[i], 4),
					round(averages[i], 4), 
					round(bin_counts[i], 4), 
					round(depths[i], 4), 
					round(total_reads[i], 4),
					sep='\t')
			if i == highest_average:
				print(color.END)

		# print("...")
		# for i in range(len(depths) -50, len(depths)):
		# 	print("",
		# 		round(p_genome[i], 4), 
		# 		round(p_reads[i], 4),
		# 		round(averages[i], 4), 
		# 		round(bin_counts[i], 4), 
		# 		round(depths[i], 4), 
		# 		round(total_reads[i], 4),
		# 		sep='\t')
	return p_genome[highest_average], depths[highest_average]






def get_kernel_coverage(bam, rgs, params, chrom_depth_c, chromosomes, out_dir):
	""" Produces coverages for alignment based off a kernel density estimation. 

	Sum of all reads within the user defined kernel bandwith are used to generate a coverage. This is meant to normalize the effect of read length on coverage.

	Returns a dict() of Counter() objects, where the dict keys are chromosomes, counter keys are chromosome positions, and counter values are depth at positions.
	"""


	## Counting read depth by position.


	cumulative_chrom = {}
	i = 0
	for c,l in chromosomes:
		cumulative_chrom[c] = i
		i += l

	cov_window = params['coverage_window']
	half_cov_window = math.floor(cov_window/2)

	# kde_window = params['kernel_window']
	# half_kde_window = int(kde_window/2)

	# pos_c    = dict()
	depth_c  = dict()
	pos_d    = dict()

	genome_length = sum([c[1] for c in chromosomes])


	print(" processing alignment...")

	for chrom, chrom_length in chromosomes:
		# pos_c[chrom]       = Counter()
		depth_c[chrom]     = Counter()
		# kernel_c[chrom]    = Counter()
		pos_d[chrom]       = dict()



	# ec = elapsedClass()
	iterables = []
	for c,l in chromosomes:
		iterables.append(samtools_view(bam, rgs=rgs, locus=c))

	reads = chain.from_iterable(iterables)

	print(f"    encoding reads ............... 0%", end='\r', flush=True)
	perc = percentageClass(1, sum(chrom_depth_c.values()))

	aligned_read_count = 0
	for i, read in enumerate(reads):
		aligned_read_count+=1
		strand, length, _, pos, chrom, _, _, _ = read

		pos += math.floor(length / 2)

		# pos_c[chrom][pos] += 1

		try:
			pos_d[chrom][pos][0] += 1
			pos_d[chrom][pos][1].append(strand)
			pos_d[chrom][pos][2].append(length)
		except KeyError:
			pos_d[chrom][pos] = [1, [strand],[length]]



		# try:
		# 	strand_d[chrom][pos].append(strand)
		# except KeyError:
		# 	strand_d[chrom][pos] = [strand]

		# try:
		# 	size_d[chrom][pos].append(length)
		# except KeyError:
		# 	size_d[chrom][pos] = [length]

		perc_out = perc.update()
		if perc_out:
			print(f"    encoding reads ............... {perc_out}%\t {i+1:,} reads", end='\r', flush=True)

		# bin_c[int((chrom_index[chrom] + pos) / params['kernel_window'])] += 1
	print()
	# print(ec)




	def inline_method():
		ec = elapsedClass()
		kernel_c = dict()
		for chrom, chrom_length in chromosomes:
			kernel_c[chrom]    = Counter()

		kernel   = deque([0]* cov_window)
		coverage = deque([0]* cov_window)
		rolling_sum = 0
		rolling_max = 0


		# print()


		gen_c = Counter()
		read_c = Counter()


		class trackClass():
			def __init__(self, bw_file, chromosomes):
				self.bw = pyBigWig.open(str(bw_file), 'w')
				self.bw.addHeader(chromosomes)

				self.last_start      = 0
				self.interval_length = 1
				self.last_chrom = chromosomes[0][0]

			def write(self):
				stop = self.last_start + self.interval_length
				self.bw.addEntries(
								[self.last_chrom], 
								[self.last_start], 
								ends= [stop], 
								values= [float(self.last_val)]
								)


			def add(self, chrom, pos, val):

				if chrom != self.last_chrom:
					self.write()
					self.last_start      = 0
					self.interval_length = 1

				elif pos > self.last_start:

					if val != self.last_val:
						self.write()
						self.last_start = pos
						self.interval_length = 1

					else:
						self.interval_length += 1


				self.last_val   = val
				self.last_chrom = chrom

			def close(self):
				self.write()
				self.bw.close()



		cov_track = trackClass(Path(out_dir, "coverage.bw"), chromosomes)
		ker_track = trackClass(Path(out_dir, "kernel.bw"), chromosomes)

		perc = percentageClass(1, sum([l for c,l in chromosomes]))
		for chrom, chrom_length in chromosomes:

			last_val = -1
			start_p  = -1


			coverage_pos = -1 * int(cov_window/2)
			kernel_pos   = -1 * cov_window

			for pos in range(chrom_length+cov_window):

				## calculating coverage

				try:
					depth = pos_d[chrom][pos][0]
				except KeyError:
					depth = 0
				
				coverage.append(depth)
				rolling_sum += depth
				rolling_sum -= coverage.popleft()

				cov_track.add(chrom, coverage_pos, rolling_sum)


				## max-smoothing coverage

				kernel.append(rolling_sum)
				outgoing_sum = kernel.popleft()


				if outgoing_sum == rolling_max:
					rolling_max = max(kernel)
				elif rolling_sum > rolling_max:
					rolling_max = rolling_sum

				ker_track.add(chrom, kernel_pos, rolling_max)



				if 0 <= kernel_pos <= chrom_length:

					gen_c[rolling_max] += 1
					try:
						read_c[rolling_max] += pos_d[chrom][kernel_pos][0]
					except KeyError:
						pass



				perc_out = perc.update()
				if perc_out:
					chrom_i = [c for c,l in chromosomes].index(chrom) + 1
					chrom_n = len(chromosomes)
					print(f"    computing coverage ........... {perc_out}%\t {chrom_i}/{chrom_n} {chrom} : {pos:,}    ", end='\r', flush=True)

				coverage_pos += 1
				kernel_pos   += 1


		cov_track.close()
		ker_track.close()

		print()
		print(ec)
		print()

		return(gen_c, read_c)


	def counter_method():

		kernel_c = dict()
		for chrom, chrom_length in chromosomes:
			kernel_c[chrom]    = Counter()


		##counter method

		print()
		print(f"    calculating coverage ......... 0%", end='\r', flush=True)
		perc = percentageClass(1, sum([len(c.values()) for c in pos_d.values()]))
		i=0
		for chrom in pos_d.keys():
			for pos, item in pos_d[chrom].items():
				depth, _, _ = item
				i += 1
				for p in range(pos-half_cov_window, pos+half_cov_window+1):
					depth_c[chrom][p] += depth

				perc_out = perc.get_percent(i)
				if perc_out:
					print(f"    calculating coverage ......... {perc_out}%\t {i:,} positions", end='\r', flush=True)

		counters_to_bigwig(depth_c, chromosomes, Path(out_dir, "coverage.bw"), verbose=False)


		# bw = pyBigWig.open(str(Path(out_dir, "coverage.bw")))

		# for chrom, chrom_length in chromosomes:
		# 	print(chrom)
		# 	for p in range(chrom_length):
		# 		try:
		# 			depth = bw.stats(chrom, p, p+cov_window, type="max")
		# 		except RuntimeError:
		# 			pass
		# 		# print(depth)
		# 		# depth = depth[0]


		gen_c = Counter()
		read_c = Counter()

		print()
		perc = percentageClass(1, genome_length)
		print(f"    calculating kernel max ....... 0%", end='\r', flush=True)

		for chrom, chrom_length in chromosomes:

			kernel = deque([0]*cov_window)
			rolling_max = 0

			for pos in range(chrom_length):

				kernel_pos = pos-half_cov_window


				depth = depth_c[chrom][pos]
				del depth_c[chrom][pos]
				kernel.append(depth)



				outgoing_avg = kernel.popleft()


				if outgoing_avg == rolling_max:
					rolling_max = max(kernel)
				elif depth > rolling_max:
					rolling_max = depth

				kernel_c[chrom][kernel_pos] = rolling_max

				gen_c[rolling_max] += 1
				try:
					read_c[rolling_max] += pos_d[chrom][kernel_pos][0]
				except KeyError:
					pass

				perc_out = perc.update()
				if perc_out:
					print(f"    calculating kernel max ....... {perc_out}%\t {chrom} : {pos:,}           ", end='\r', flush=True)

		# pprint(depth_c)

		print(" writing to bigwig file...")
		counters_to_bigwig(kernel_c, chromosomes, Path(out_dir, "kernel.bw"), verbose=False)

		return(gen_c, read_c)


	# print()
	# print(f"    tallying thresholds .......... 0%", end='\r', flush=True)
	# perc = percentageClass(1, genome_length)
	# measured_genome = 0
	# # gen_c = Counter()
	# # read_c = Counter()
	# found_depths = set()
	# for c,l in chromosomes:
	# 	for p in range(l):
	# 		measured_genome += 1
	# 		depth = kernel_c[c][p]
	# 		# gen_c[depth] += 1
	# 		# read_c[depth] += pos_c[c][p]
	# 		found_depths.add(depth)


	# 		perc_out = perc.update()
	# 		if perc_out:
	# 			print(f"    tallying thresholds .......... {perc_out}%           ", end='\r', flush=True)

	# gen_c, read_c = inline_method()

	def tally(gen_c, read_c):
		found_depths = list(gen_c.keys())
		found_depths.sort()

		averages = []
		table    = []

		total_genomic_space = genome_length
		total_read_space    = sum(read_c.values())

		for depth_threshold in found_depths:

			total_genomic_space -= gen_c[depth_threshold]
			total_read_space    -= read_c[depth_threshold]

			# print(depth_threshold, total_genomic_space, sep='\t')

			gen_score = total_genomic_space/genome_length
			read_score = total_read_space/aligned_read_count
			avg_score = ((1-gen_score) + read_score) /2

			# if not out and avg_score < last_avg:
			# 	peak = "1"
			# 	out = {
			# 		'depth_threshold' : depth_threshold,
			# 		'gen_score' : gen_score,
			# 		'read_score' : read_score,
			# 		'avg_score' : avg_score
			# 	}

			# else:
			# 	peak = '0'
			# last_avg = avg_score

			averages.append(avg_score)

			table.append([depth_threshold, total_genomic_space, round(gen_score,4),
				total_read_space, round(total_read_space/aligned_read_count,4),
				round(avg_score, 4)])

		peak_index = averages.index(max(averages))

		with open(Path(out_dir, 'thresholds.txt'), 'w') as outf:
			print('depth\tannotated_space\tp_genome\tannotated_reads\tp_reads\taverage_score\tpeak', file=outf)
			for i,t in enumerate(table):

				if i == peak_index:
					peak = 1
					out = {
						'depth_threshold' : t[0],
						'gen_score' : t[2],
						'read_score' : t[4],
						'avg_score' : t[5]
						}
				else:
					peak = 0

				print("\t".join(map(str, t)), peak, sep='\t', file=outf)
				if total_genomic_space > genome_length:
					sys.exit("problem!!")

		# pprint(out)
		return(out)
	
	# print()
	# print(f'   pos_d: {asizeof.asizeof(pos_d):,}')
	# print(f' depth_c: {asizeof.asizeof(depth_c):,}')
	# print(f'kernel_c: {asizeof.asizeof(kernel_c):,}')

	gen_c, read_c = inline_method()
	out = tally(gen_c, read_c)
	# print(out)


	# gen_c, read_c = counter_method()
	# out = tally(gen_c, read_c)
	# print(out)
	# sys.exit()



	return(pos_d, out)












@cli.command(group='Annotation', help_priority=2)



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

@optgroup.option("-n", "--name", 
	# default=f"Annotation_{round(time())}", 
	required=False,
	type=str,
	help="Optional name alignment. Useful if comparing annotations.")



@optgroup.group('\n  Coverage options',
				help='')

@optgroup.option('--coverage_window',
	default=250,
	help='This is the bandwidth for accumulating read alignments into coverage. This differs from common methods like samtools depth to avoid bias for longer reads. Default 40 nt.')



@optgroup.option('--subsample',
	help="Allows the user to subsample alignments for the annotation to a defined depth. Accepts an integer number of reads, which can be modified with a 10^3 prefix (ex. 10M).")

@optgroup.option('--subsample_seed',
	type=int,
	default=0,
	help="Seed value used for subsampling (default: 0)")

@optgroup.option('--subsample_n',
	type=int,
	default=0,
	help="The index of which split group from the subsets you want to use for the annotation. For example, a 105M deep alignment will be split into 5 distinct sets when subset by 20M (residual 5M are ignored). This option which pick which to use (base-0)")


@optgroup.option('--subsample_keep_max',
	type=int,
	default=0,
	help="The maximum number of subset alignments that will be written to the disk. Numbers higher than 1 are really only useful for performance comparisons. This value will automatically be raised to a minimum of the subsample_n+1.")





@optgroup.group('\n  Peak finding options',
				help='')

# @optgroup.option('--kernel_window',
# 	default=500,
# 	help='This is the bandwidth for coverage smoothing. Within this window, coverages are smoothed by mean and susequently used to find region boundaries. Default 50 nt.')



# @optgroup.option('--bin_size',
# 	default=50,
# 	help='')


@optgroup.option('--genome_perc',
	default=False,
	help='')
@optgroup.option('--read_perc',
	default=False,
	help='')

@optgroup.option('--trim_regions', is_flag=True, default=False, help='Flag to include trimming of regions')





@optgroup.group('\n  Merging options',
				help='')


@optgroup.option("--clump_dist",
	default=500,
	help="Distance in nucleotides for which sRNA peaks should be considered for 'clumping'. Clumped regions must have sufficient similarity in sRNA-size profile and strand-preference. Default 500 nt.")

@optgroup.option("--clump_strand_similarity",
	default=0.7,
	help="Similarity threshold of strand fraction for clumping two peaks. Difference in fraction must be smaller than threshold. Default 0.5.")

@optgroup.option("--min_locus_length",
	default=30,
	help="Minimum size for a locus to be included. Default 30 nt.")




@optgroup.group('\n Other options',
				help='')

@optgroup.option('--debug', is_flag=True, default=False, help='Debug flag')
@optgroup.option('--test_mode', is_flag=True, default=False, help='test_mode flag')
@optgroup.option('--override', is_flag=True, default=False, help='Overrides config file changes without prompting.')


def tradeoff(**params):
	'''Annotator using large coverage window and prec/sens tradeoff.'''

	rc = requirementClass()
	rc.add_samtools()
	rc.check()

	ic = inputClass(params)
	ic.check(['alignment_file', 'annotation_readgroups'])

	output_directory        = str(ic.output_directory)
	alignment_file          = ic.inputs["alignment_file"]
	annotation_readgroups   = ic.inputs['annotation_readgroups']
	project_name            = ic.inputs['project_name']

	# coverage_window         = params['coverage_window']
	# kernel_window           = params['kernel_window']
	# peak_threshold          = params['peak_threshold']
	# rpm_threshold           = params['rpm_threshold']
	# pad                     = params['pad']
	clump_dist              = params['clump_dist']
	clump_strand_similarity = params['clump_strand_similarity']
	min_locus_length        = params['min_locus_length']
	debug                   = params['debug']
	annotation_name         = params['name']
	target_depth            = params['subsample']
	seed                    = params['subsample_seed']


	params['output_directory'] = output_directory
	params['alignment_file'] = alignment_file
	params['project_name'] = project_name
	pprint(params)

	if annotation_name:
		dir_name = f"tradeoff_{annotation_name}"
	else:
		dir_name = 'tradeoff'




	if debug: 
		show_warnings = True
	else:
		show_warnings = False


	Path(output_directory, dir_name).mkdir(parents=True, exist_ok=True)

	log_file = f"{output_directory}/{dir_name}/log.txt"
	sys.stdout = Logger(log_file)


	start_time = datetime.now()
	date_time = start_time.strftime("%Y/%m/%d, %H:%M:%S")
	print()
	print("Start time:",date_time)	

	params['tool'] = 'tradeoff'
	parameters_file = Path(output_directory, dir_name, "params.json")
	with open(parameters_file, 'w') as outf:
		oparams = params
		oparams['alignment_file'] = str(params['alignment_file'])
		outf.write(json.dumps(params, indent=2))
		del oparams
		# sys.exit()

	cl_i = 0

	# assert 0 < sensitivity <= 1.0, "sensitivity must be 0 < S <= 1.0"
	# assert 0 < peak_threshold < 1, "peak_trim must be between 0 and 1."

	# half_window = int(window/2)

	chromosomes, bam_rgs = get_chromosomes(alignment_file)

	if params['test_mode']:
		# chromosomes = chromosomes[3:5]
		# chromosomes = chromosomes[6:10]
		# chromosomes = chromosomes[2:5]
		chromosomes = chromosomes[:2]
		# chromosomes = chromosomes[:1]

	genome_length = sum([l for c,l in chromosomes])

	# print(alignment_file)
	annotation_readgroups = check_rgs(annotation_readgroups, bam_rgs)
	params['annotation_readgroups'] = annotation_readgroups

	print()
	print(" Annotation readgroups:", annotation_readgroups)

	# chromosomes = [c for c in chromosomes if c[0] == 'NC_037320.1']

	chrom_depth_c = get_global_depth(output_directory, alignment_file, aggregate_by=['rg','chrom'])

	keys = list(chrom_depth_c.keys())
	for key in keys:
		if key[0] in annotation_readgroups:
			chrom_depth_c[key[1]] += chrom_depth_c[key]

		del chrom_depth_c[key]

	for key in list(chrom_depth_c.keys()):
		if key not in [c for c,l in chromosomes]:
			del chrom_depth_c[key]

	aligned_read_count = sum(chrom_depth_c.values())




	if target_depth:
		subsample = parse_subsample(target_depth, alignment_file, "bam", aligned_read_count, seed=seed,
			n=params['subsample_n'])

		alignment_file = perform_subsample(subsample, subsample_keep_max=params['subsample_keep_max'])

		chrom_depth_c = get_global_depth(output_directory, alignment_file, aggregate_by=['rg','chrom'])

		keys = list(chrom_depth_c.keys())
		for key in keys:
			if key[0] in annotation_readgroups:
				chrom_depth_c[key[1]] += chrom_depth_c[key]

			del chrom_depth_c[key]

		for key in list(chrom_depth_c.keys()):
			if key not in [c for c,l in chromosomes]:
				del chrom_depth_c[key]

		aligned_read_count = subsample.target


	def get_frac_top(c):
		try:
			ftop = c['+'] / sum(c.values())
		except ZeroDivisionError:
			ftop = -1

		return(ftop)

	def get_basic_dimensions(l):
		name, chrom, start, stop = l
		coords = f"{chrom}:{start}-{stop}"


		reads = [r for r in samtools_view(alignment_file, 
			locus=coords, 
			rgs=annotation_readgroups, 
			boundary_rule = 'tight')]


		strand_c = Counter()
		lens = []

		for read in reads:
			sam_strand, sam_length, _, _, _, _, _, _ = read

			strand_c[sam_strand] += 1
			lens.append(sam_length)

		try:
			frac_top = strand_c["+"] / sum(strand_c.values())
		except ZeroDivisionError:
			frac_top = -1

		# try:
		# 	most_common = length_c.most_common(1)[0][0]
		# 	frac_size = length_c[most_common] / sum(strand_c.values())
		# except ZeroDivisionError:
		# 	most_common = False
		# 	frac_size = -1
		
		return(frac_top, lens)



	## preparing output files

	def init_gff(file_name):
		with open(file_name, 'w') as outf:
			print("##gff-version 3", file=outf)

			for chrom, chrom_length in chromosomes:
				print(f"##sequence-region   {chrom} 1 {chrom_length}", file=outf)

	gff_file = Path(output_directory, dir_name, 'loci.gff3')
	init_gff(gff_file)

	revised_gff = Path(output_directory, dir_name, "revised_regions.gff3")
	init_gff(revised_gff)

	region_file = Path(output_directory, dir_name, 'regions.gff3')
	init_gff(region_file)


	results_file = Path(output_directory, dir_name, 'loci.txt')
	with open(results_file, 'w') as outf:
		print("\t".join(assessClass().header), file=outf)


	reads_file = Path(output_directory, dir_name, 'reads.txt')
	with open(reads_file, 'w') as outf:
		print(TOP_READS_HEADER, file=outf)


	merge_file = Path(output_directory, dir_name, 'merges.txt')
	with open(merge_file, 'w') as outf:
		outf.write('')


	stats_file = f"{output_directory}/{dir_name}/stats_by_chrom.txt"
	with open(stats_file, 'w') as outf:
		print("chromosome\tregions\tregions\tchromosome_length\tproportion_genome_annotated\tmean_length\tmedian_length\treads\tproportion_libraries_annotated\tmean_abundance\tmedian_abundance", file=outf)


	overall_file = f"{output_directory}/{dir_name}/stats.txt"
	with open(overall_file, 'w') as outf:
		outf.write('')

	overall_d = {}

	all_regions = []
	total_read_count = 0
	cl_i = 0


	# total_reads = sum(chrom_depth_c.values())
	read_equivalent = 1 / sum(chrom_depth_c.values()) * 1000000


	## iterating through chromosomes and reads

	print()
	print(f" {len(chromosomes)} chromosomes/scaffolds/contigs")
	print(f" {sum([l for c,l in chromosomes]):,} bp total")
	print(f" {aligned_read_count:,} reads")
	print()





	# cdf = get_bin_coverage_2(alignment_file=alignment_file, chromosomes=chromosomes, 
	# 	bin_size = params['bin_size'])




	## Getting positional coverage accross the alignment



	# if len(all_regions) == 0:

	start = time()

	pos_d, threshold_stats = get_kernel_coverage(
		bam=alignment_file, 
		rgs=annotation_readgroups, 
		params=params, 
		chrom_depth_c=chrom_depth_c,
		chromosomes=chromosomes,
		out_dir=Path(output_directory, dir_name))

	# sys.exit()

	depth_threshold = threshold_stats['depth_threshold']
	gen_score       = threshold_stats['gen_score']
	read_score      = threshold_stats['read_score']


	print()
	print()
	print(" annotation parameters...")
	print(f"    depth threshold: ......... {depth_threshold} reads")
	print(f"    exp. genome proportion: .. {gen_score}")
	print(f"    exp. read proportion: .... {read_score}")


	def get_regions_dep(kernel_c, depth_threshold, chromosomes):

		regions = []
		reg_i = 0
		for chrom, chrom_length in chromosomes:
			# print(chrom)

			positions = []
			for p,d in kernel_c[chrom].most_common():
				if d < depth_threshold:
					break
				positions.append(p)

			positions.sort()

			last_p = None
			for p in positions:

				if not last_p:
					start = p

				elif p > last_p + 1:
					reg_i += 1
					if start < last_p:
						regions.append([f"region_{reg_i}", chrom, start, last_p])
					start = p


				last_p = p



			if last_p:

				if last_p > chrom_length - params['coverage_window']:
					last_p = chrom_length
				reg_i += 1
				if start < last_p:
					regions.append([f"region_{reg_i}", chrom, start, last_p])




		return(regions)

	def get_regions_dep2(depth_threshold, chromosomes):

		regions = []

		reg_i = 0
		bw = pyBigWig.open(str(Path(output_directory, dir_name, "kernel.bw")))
		# print(bw.isBigWig())
		# print(bw.chroms())
		# print(bw.header())

		for chrom, chrom_length in chromosomes:
			# print(chrom)

			reg_start = False
			reg_stop  = False

			for inv_start, inv_stop, depth in bw.intervals(chrom, 0, chrom_length):

				if depth >= depth_threshold:

					if inv_start != reg_stop:

						if reg_start:
							reg_i += 1
							region = [f"region_{reg_i}", chrom, reg_start, reg_stop]
							# print(region)
							regions.append(region)

						reg_start = inv_start
					reg_stop  = inv_stop
				# last_stop = inv_stop



		bw.close()


		return(regions)


	def get_regions(depth_threshold, chromosomes):

		regions = []

		reg_i = 0
		bw = pyBigWig.open(str(Path(output_directory, dir_name, "kernel.bw")))
		# print(bw.isBigWig())
		# print(bw.chroms())
		# print(bw.header())

		for chrom, chrom_length in chromosomes:
			# print(chrom)

			in_region = False

			for inv_start, inv_stop, depth in bw.intervals(chrom, 0, chrom_length):
				if depth >= depth_threshold:
					if not in_region:
						reg_start = inv_start
					reg_stop  = inv_stop

					in_region = True

				else:
					if in_region:
						reg_i += 1
						region = [f"region_{reg_i}", chrom, reg_start, reg_stop]
						regions.append(region)

					in_region = False




		bw.close()


		return(regions)



	print()
	print(' finding regions...')
	all_regions = get_regions(depth_threshold, chromosomes)

	# pprint(regions)


	# all_regions = [
	# 	['test_1', 'NC_037310.1', 2877634, 2880258],
	# 	['test_2', 'NC_037310.1', 2815538, 2815794],
	# 	['locus_68', 'NC_037310.1', 2748416, 2749269]
	# ]

	print(f"    {len(all_regions):,} regions")
	print()

	def get_region_space(regions):
		tot = 0
		for name, chrom, start, stop in all_regions:
			length = stop - start

			tot += length
			# print(start, stop, length, tot, sep='\t')
		return tot

	total_region_space = get_region_space(all_regions)
	print(f"    {total_region_space:,} ({round(total_region_space/genome_length *100,1)}%) genome in regions")
	print(f"        expected: {round(100*threshold_stats['gen_score'],1)}%")



	total_annotated_reads = 0
	for l in all_regions:
		# print(l)

		name, chrom, start, stop = l
		for r in range(start, stop+1):
			try:
				total_annotated_reads += pos_d[chrom][r][0]
			except KeyError:
				pass



	print(f"    {total_annotated_reads:,} ({round(total_annotated_reads/aligned_read_count *100,1)}%) total reads in regions")
	print(f"        expected: {round(100*threshold_stats['read_score'],1)}%")




	def write_regions_to_file():
		print()
		print(" writing regions to gff file...")
		with open(region_file, 'a') as outf:

			for l in all_regions:

				name, chrom, start, stop = l
				if start < 1:
					start = 1

				line = [chrom, 'yasma_peak','region', start, stop, '.', '.', '.',
			f'ID={name}']


				line = "\t".join(map(str,line))
				print(line, file=outf)

				# reads = samtools_view(alignment_file, locus=f"{chrom}:{start}-{stop}", rgs=annotation_readgroups)
				# for read in reads:
				# 	total_annotated_reads += 1

	write_regions_to_file()


	now = datetime.now()

	# del peak_c
	# del pos_c
	# del kernel_c


	def expand_region(claim_d, region):

		locus_name, chrom, start, stop = region
		coords = f"{chrom}:{start}-{stop}"

		# reads = samtools_view(alignment_file, rgs=annotation_readgroups, locus=coords)


		strands = Counter()
		sizes   = sizeClass()
		names   = set()


		for p in range(start, stop+1):
			try:
				strands.update(pos_d[chrom][p][1])
				sizes.update(pos_d[chrom][p][2])
			except:
				pass


		try:
			frac_top = strands['+'] / sum(strands.values())
		except ZeroDivisionError:
			print(region)

			sys.exit()

		cov_window = params['coverage_window']

		region_size = stop - start


		# print()

		def expand_window(window_l, window_r, direction, increment=50, max_increments=False):

			window_size = abs(window_l - window_r)

			if direction == -1:
				i = 0
			else:
				i = 1

			increment_count =0
			while True:

				increment_count += 1
				if max_increments and increment_count > max_increments:
					break

				window_l += increment * direction
				window_r += increment * direction

				w_strands = Counter()
				w_sizes   = sizeClass()

				for w in range(window_l, window_r):
					try:
						w_strands.update(pos_d[chrom][w][1])
						w_sizes.update(pos_d[chrom][w][2])
					except:
						pass



				gen = [range(window_r, window_l,-1), range(window_l, window_r)][i]
				claimed = False
				for r in gen:
					if r in claim_d[chrom]:
						if  name != claim_d[chrom][r]:
							if i == 0:
								window_l = r + 1
							else:
								window_r = r - 1

							claimed=True
							break



				try:
					w_fractop = w_strands['+'] / sum(w_strands.values())
				except ZeroDivisionError:
					w_fractop = 1




				if claimed:
					# print('claim break', end='\t')
					break

				if sum(w_strands.values()) == 0:
					# pprint(w_strands)
					# print('empty break', end='\t')
					# sys.exit()
					break
					
				if sum(w_strands.values())/window_size < sum(strands.values())/region_size * 0.05:
					# print('depth break', end='\t')
					break

				if not w_sizes == sizes:
					# print('size break', end='\t')
					break

				if not abs(frac_top - w_fractop) < 0.5:
					# print('strand break', end='\t')
					break


			if direction == -1:
				return(window_l)
			else:				
				return(window_r)


		# for direction in enumerate([-1,1]):


		next_window_l = expand_window(start, start+cov_window, -1, increment=50)
		next_window_r = expand_window(stop-cov_window, stop, 1, increment=50)

		if next_window_l < next_window_r:
			window_l, window_r = next_window_l, next_window_r

		# print(window_l, window_r, end='\t\t')

		next_window_l = expand_window(window_l+cov_window, window_l+cov_window-50, -1, increment=15)
		next_window_r = expand_window(window_r-cov_window, window_r-cov_window+50, 1,  increment=15)
		# print(window_l, window_r)

		if next_window_l < next_window_r:
			window_l, window_r = next_window_l, next_window_r

		for w in range(window_l, window_r+1):
			claim_d[chrom][w] = locus_name


		
		return(window_l, window_r, claim_d)


	claim_d = {}
	for c, l in chromosomes:
		claim_d[c] = {}

	for region in all_regions:
		name, chrom, start, stop = region
		for r in range(start, stop+1):
			claim_d[chrom][r] = name



	
	print()

	print(' revising regions...', end='\r')

	revised_genomic_space = 0
	perc = percentageClass(increment=5, total=len(all_regions))

	with open(revised_gff, 'a') as outf:

		for region in all_regions:
			name, chrom, start, stop = region

			new_start, new_stop, claim_d = expand_region(claim_d, region)

			if new_start > new_stop:
				print(region)
				print(new_start, new_stop)
				print("expand failed!! stop before start")
				sys.exit()

			revised_genomic_space += new_stop - new_start

			if new_start < 1:
				new_start = 1


			region[2] = new_start
			region[3] = new_stop


			gff_line = [
				chrom, 'test_region','region', new_start, new_stop, '.', '.', '.',
				f'ID={name}']

			print('\t'.join(map(str,gff_line)), file=outf)

			perc_out = perc.update()
			if perc_out:
				print(f' revising regions ... {perc_out}%  ', end='\r', flush=True)

		# if name == 'region_2':
		# 	sys.exit()


	total_revised_reads = 0
	for l in all_regions:
		# print(l)

		name, chrom, start, stop = l
		for r in range(start, stop+1):
			try:
				total_revised_reads += pos_d[chrom][r][0]
			except KeyError:
				pass

	print()

	print(f"    {revised_genomic_space:,} ({round(revised_genomic_space/genome_length *100,1)}%) genome in revised regions")
	print(f"    {total_revised_reads:,} ({round(total_revised_reads/aligned_read_count *100,1)}%) total reads in revised regions")
	diff = round(total_revised_reads/aligned_read_count *100,1) - round(total_annotated_reads/aligned_read_count *100,1)
	print(f"      +{round(diff,2)}% over unrevised regions")

	max_chrom_word_length = max([len(c) for c,l in chromosomes])
	def print_progress_string(i, n, chrom, input_loci, output_loci, assess=False):

		chrom = chrom + (max_chrom_word_length - len(chrom)) * " "
		if not assess:
			assess = ''
		else:
			assess = f"\t{assess}%"
		print(f"{i+1}/{n}\t{chrom}\t{unclumped_regions_count}\t{len(regions)}{assess}", end='\r', flush=True)

	print()
	print('prog', "chrom"+(max_chrom_word_length-5)*" ",'regions\tloci\tassess', sep='\t')

	total_region_space = 0
	regions_name_i = 0
	total_annotated_reads = 0

	for chrom_count, chrom_and_length in enumerate(chromosomes):


		chrom, chrom_length = chrom_and_length

		status_line = f'{chrom_count+1}/{len(chromosomes)}\t'



		regions = [r for r in all_regions if r[1] == chrom]


		## Performing clumping by reading whole chromosome alignments

		clump_set = set()


		def get_considered_regions(i):

			if i >= len(regions):
				return(False)

			out = [i]
			n = 1
			current_end = regions[i][3]
			while i+n < len(regions):

				# try:
				next_start = regions[i+n][2]
				# except:

				# 	print(regions[i+n])
				# 	sys.exit()

				if next_start < current_end + clump_dist:
					out.append(i+n)

				else:
					break

				n += 1

			return(out)





		# regions_copy = regions[:]


		# start = time()


		locus_i = 0
		considered_regions = get_considered_regions(locus_i)
		unclumped_regions_count = len(regions)


		print_progress_string(chrom_count, len(chromosomes), chrom, unclumped_regions_count, len(regions))
		# print(f"{chrom_count+1}/{len(chromosomes)}\t{chrom}\t100%\t{unclumped_regions_count} -> {len(regions)}", end='\r', flush=True)
		# print(f"   clumping similar neighbors... {unclumped_regions_count} -> {len(regions)} regions    ", end='\r', flush=True)

		last_claim = 'start'

		# test_d = {}
		in_locus = False


		sizecall_d = {}
		strand_d   = {}
		seq_d      = {}
		n = False


		def final_check_and_increment(i):
			locus_length = regions[i][3] - regions[i][2]
			locus_depth  = sum(strand_d[regions[i][0]].values())
			if locus_length < min_locus_length:
				del regions[i]
			# elif locus_depth / locus_length < params['min_locus_resolution']:
			# 	del regions[i]
			else:
				i += 1
			return(i)






		for read in samtools_view(alignment_file, locus=chrom, rgs=annotation_readgroups):

			## Breaks for the final region
			if not considered_regions:
				break


			## Processing sam output
			sam_strand, sam_length, _, sam_lbound, _, _, sam_seq, read_name = read
			sam_rbound = sam_lbound + sam_length


			## Identifies the current claim

			try:
				rclaim = claim_d[chrom][sam_rbound]
			except KeyError:
				rclaim = None

			try:
				lclaim = claim_d[chrom][sam_lbound]
			except KeyError:
				lclaim = None


			if not lclaim or not rclaim:
				claim = f"after_{last_claim}"

			elif lclaim != rclaim:
				if lclaim:
					last_claim = lclaim
				claim = f"after_{last_claim}"

			else:
				in_locus = True
				claim = lclaim
				last_claim = lclaim



			## logs the claim to the data containers
			try:
				sizecall_d[claim]
			except KeyError:
				sizecall_d[claim] = sizeClass()
				strand_d[claim] = Counter()
				seq_d[claim] = Counter()


			## Adding values to current claim
			if in_locus:


				sizecall_d[claim].update([sam_length])
				strand_d[claim].update([sam_strand])
				seq_d[claim].update([sam_seq])


			verbose = False



			if sam_lbound > regions[considered_regions[-1]][3]:
				## Checking for any broken regions

				for c in considered_regions[::-1]:
					c_name = regions[c][0]

					if c_name not in strand_d:
						to_loc   = f"after_{c_name}"
						from_loc = f"after_{regions[c-1][0]}"

						try:
							sizecall_d[to_loc] += sizecall_d[from_loc]
							strand_d[to_loc] += strand_d[from_loc]
							seq_d[to_loc] += seq_d[from_loc]
						except KeyError:
							pass

						del regions[c]

						if show_warnings:
							print(f"Warning: bad locus {c_name} removed")


				considered_regions = get_considered_regions(locus_i)




			if not considered_regions:
				break

			# try:
			# 	sam_lbound > regions[considered_regions[-1]][3]
			# except:

			# 	print()
			# 	print("all regions:")
			# 	print(regions)
			# 	print()
			# 	print("considered_regions:")
			# 	print(considered_regions)
			# 	print()
			# 	print("response:")
			# 	print(regions[considered_regions[-1]])
			# 	print("Warning: unknown problem in locus merging...")
			# 	sys.exit()
			# 	break

			if sam_lbound > regions[considered_regions[-1]][3]:
				## if all good, moving forward



				if len(considered_regions) == 1:

					with open(merge_file, 'a') as outf:
						print('', file=outf)
						pprint([regions[c] for c in considered_regions], outf)
						print('      -> locus has no other regions in range', file=outf)

					locus_i = final_check_and_increment(locus_i)
					# print("\nonly_one_considered <- increment")

					considered_regions = get_considered_regions(locus_i)



				else:

					n = 0
					while True:
						with open(merge_file, 'a') as outf:
							print('', file=outf)
							pprint([regions[c] for c in considered_regions], outf)
							# input()

						none_merged = True
						
						current_locus = regions[locus_i]

						# print(considered_regions)
						for n in considered_regions[1:]:
							next_locus = regions[n]

							# print()
							# print()
							# print(current_locus)
							# print(next_locus)
							# print(sam_lbound)
							# print()

							try:
								current_sc = sizecall_d[current_locus[0]]
								current_ft = get_frac_top(strand_d[current_locus[0]])
							except KeyError:
								current_sc = sizeClass()
								current_ft = -1
								# print(f"\nWarning: region {current_locus} may not have counts")


							try:
								next_sc = sizecall_d[next_locus[0]]
								next_ft = get_frac_top(strand_d[next_locus[0]])
							except KeyError:
								next_sc = sizeClass()
								next_ft = -1
								# print(f"\nWarning: region {next_locus} may not have counts")





							size_test = current_sc == next_sc
							frac_test = abs(current_ft - next_ft) < clump_strand_similarity
							# print(size_test, current_sc, next_sc, sep='\t')
							# print(frac_test, current_ft, next_ft, sep='\t')

							acc_to_merge = size_test and frac_test





							## Code to test if the methods produce equivalent results

							# sam_sc1, sam_ft1, sizes1, strands1 = test_by_samtools(current_locus)
							# sam_sc2, sam_ft2, sizes2, strands2 = test_by_samtools(next_locus)

							# sam_to_merge = sam_sc1 == sam_sc2 and abs(sam_ft1 - sam_ft2) < clump_strand_similarity

							# if acc_to_merge and not sam_to_merge:
							# 	print()
							# 	print()
							# 	print(current_locus)
							# 	print(next_locus)
							# 	print(sam_lbound, "<- sam_lbound")
							# 	print()

							# 	print("acc current:")
							# 	print(current_sc)
							# 	print(current_sc.depth)
							# 	print(current_ft)
							# 	print()
							# 	print("sam current:")
							# 	print(sam_sc1)
							# 	print(sam_sc1.depth)
							# 	print(sam_ft1)
							# 	print()
							# 	print("acc next:")
							# 	print(next_sc)
							# 	print(next_sc.depth)
							# 	print(next_ft)
							# 	print()
							# 	print("sam next:")
							# 	print(sam_sc2)
							# 	print(sam_sc2.depth)
							# 	print(sam_ft2)
							# 	print()

							# 	print(size_test, "<- size_test")
							# 	print(frac_test, "<- frac_test")
							# 	print()

							# 	input()

							# if current_locus[0] == 'Cluster_401' or next_locus[0] == 'Cluster_401' or last_one:
							# 	print()
							# 	print(locus_i)
							# 	last_one = True
							# 	print()
							# 	print()
							# 	print(current_locus)
							# 	print(next_locus)
							# 	print(size_test, frac_test)

							# 	print(strand_d[current_locus[0]])
							# 	print(strand_d[next_locus[0]])
							# 	input()
							# else:

							# 	last_one=False

							with open(merge_file, 'a') as outf:
								print(current_locus, "->", next_locus, file=outf)
								print("", file=outf)
								print("", size_test, current_sc, next_sc, sep='\t', file=outf)
								print("", frac_test, round(current_ft,4), round(next_ft,4), sep='\t', file=outf)
								print("", file=outf)

							if size_test and frac_test:
								## Then merge!

								none_merged = False

								with open(merge_file, 'a') as outf:
									print(f"merging {current_locus[0]} to {next_locus[0]}", file=outf)
									print("", file=outf)

									clump_set.add((current_locus[0], next_locus[0]))


									## Packing counts into first locus

									for r in range(locus_i, n):

										to_loc   = current_locus[0]
										from_loc = regions[r+1][0]


										try:
											sizecall_d[to_loc] += sizecall_d[from_loc]
											strand_d[to_loc] += strand_d[from_loc]
											seq_d[to_loc] += seq_d[from_loc]
											print(f"  packing -> {to_loc} with {sizecall_d[from_loc].depth} reads from {from_loc}", file=outf)
										except KeyError:
											print("Warning: unknown KeyError in locus packing...")
											# sys.exit("error!")
											# print(f"\nWarning: region {next_locus} may not have counts")

											# del regions[locus_i]
											# considered_regions = get_considered_regions(locus_i) 
											pass

										to_loc   = current_locus[0]
										from_loc = f"after_{regions[r][0]}"

										try:
											sizecall_d[to_loc] += sizecall_d[from_loc]
											strand_d[to_loc] += strand_d[from_loc]
											seq_d[to_loc] += seq_d[from_loc]
											print(f"  packing -> {to_loc} with {sizecall_d[from_loc].depth} reads from {from_loc}", file=outf)
										except KeyError:
											# sys.exit("error!")
											# print(f"\nWarning: region {next_locus} may not have counts")
											# del regions[locus_i]
											# considered_regions = get_considered_regions(locus_i) 
											pass

									to_loc   = f"after_{current_locus[0]}"
									from_loc = f"after_{regions[n][0]}"

									try:
										sizecall_d[to_loc] += sizecall_d[from_loc]
										strand_d[to_loc] += strand_d[from_loc]
										seq_d[to_loc] += seq_d[from_loc]
										print(f"  packing -> {to_loc} with {sizecall_d[from_loc].depth} reads from {from_loc}", file=outf)
									except KeyError:
										# print(f"\nWarning: region {next_locus} may not have counts")

										# sys.exit("error!")
										# del regions[locus_i]
										# considered_regions = get_considered_regions(locus_i) 
										pass



								## Redefining current locus boundaries
								regions[locus_i][3] = regions[n][3]


								## Eliminating regions up to the one that was merged
								for r in range(locus_i, n)[::-1]:
									del regions[r+1]


								## Filling in claim_d
								for r in range(regions[locus_i][2],regions[locus_i][3]+1):
									claim_d[chrom][r] = regions[locus_i][0]


								break


						## updating considered regions after merging
						considered_regions = get_considered_regions(locus_i) 

						if none_merged:
							with open(merge_file, 'a') as outf:
								print("      -> no valid merges for considered_regions", file=outf)
							locus_i = final_check_and_increment(locus_i)
							# print("\nnone_merged <- increment")
							considered_regions = get_considered_regions(locus_i) 
							break

						if len(considered_regions) == 1:
							with open(merge_file, 'a') as outf:
								print('      -> new locus has no other regions in range', file=outf)
							locus_i = final_check_and_increment(locus_i)
							# print("\nlen(considered_regions) <- increment")
							considered_regions = get_considered_regions(locus_i) 
							break

						if sam_lbound < regions[considered_regions[-1]][3]:
							with open(merge_file, 'a') as outf:
								print("      -> regions passed read location", file=outf)
							last_claim = regions[considered_regions[-1]][0]
							break




				# print(f"   clumping similar neighbors... {unclumped_regions_count} -> {len(regions)} regions    ", end='\r', flush=True)

				# print(f"{chrom_count+1}/{len(chromosomes)}\t{chrom}\t100%\t{unclumped_regions_count} -> {len(regions)}      ", end='\r', flush=True)

				print_progress_string(chrom_count, len(chromosomes), chrom, unclumped_regions_count, len(regions))


		# print()
		to_delete = []
		for i,locus in enumerate(regions):
			name = locus[0]

			if name not in strand_d:
				print("Serious problem - regions with no reads broke through:", name)

				to_delete.append(i)

		for i in to_delete[::-1]:
			del regions[i]




		### Trimming regions

		if params['trim_regions']:
			for i,locus in enumerate(regions):
				locus[2:] = trim_locus(locus, peak_c[locus[1]], params)





		## Assessing locus dimensions and making annotations

		# read_depths = []

		perc = percentageClass(increment=5, total=len(regions))


		last_stop = 0
		for i,locus in enumerate(regions):

			old_name = locus[0]
			regions_name_i += 1
			locus[0] = f"locus_{regions_name_i}"


			print_percentage = perc.get_percent(i)
			if print_percentage:

				# print(f"{chrom_count+1}/{len(chromosomes)}\t{chrom}\t100%\t{unclumped_regions_count} -> {len(regions)}\t{print_percentage}%", end='\r', flush=True)

				print_progress_string(chrom_count, len(chromosomes), chrom, unclumped_regions_count, len(regions), print_percentage)
				# print(f"   assessing regions .............. {print_percentage}%", end="\r", flush=True)

			name, chrom, start, stop = locus
			coords = f"{chrom}:{start}-{stop}"

			# print(locus)

			strand_d[name]   = strand_d.pop(old_name)
			sizecall_d[name] = sizecall_d.pop(old_name)
			seq_d[name]      = seq_d.pop(old_name)

			strand_c = strand_d[name]
			sizecall = sizecall_d[name]
			read_c   = seq_d[name]

			# reads = [r for r in samtools_view(alignment_file, locus=coords, rgs=annotation_readgroups)]

			# read_depths.append(len(reads))

			# read_c = Counter()
			# for read in reads:
			# 	sam_strand, _, _, _, _, _, sam_read, _ = read

			# 	if sam_strand == '-':
			# 		sam_read = complement(sam_read[::-1])

			# 	read_c[sam_read] += 1



			
			# results_line, gff_line = assessClass().format(name, chrom, start, stop, reads, sum(chrom_depth_c.values()), last_stop)

			results_line, gff_line = assessClass().format(locus, read_c, strand_c, sizecall, sum(chrom_depth_c.values()), last_stop)


			last_stop = stop

			with open(results_file, 'a') as outf:
				print("\t".join(map(str, results_line)), file=outf)

			with open(gff_file, 'a') as outf:
				print("\t".join(map(str, gff_line)), file=outf)

			top_reads_save(read_c, reads_file, read_equivalent, name)

			# print(results_line)
			# print(gff_line)

			# print(locus)
			# pprint(read_c)
			# pprint(strand_c)
			# print(sizecall)
			# sys.exit()


		if regions:
			print_progress_string(chrom_count, len(chromosomes), chrom, unclumped_regions_count, len(regions), print_percentage)

		# all_regions += regions
		# sys.exit()

		# print()
		# print()


		locus_lengths = [l[3]-l[2] for l in regions]
		try:
			mean_length   = round(mean(locus_lengths),1)
		except StatisticsError:
			mean_length = None
		try:
			median_length = median(locus_lengths)
		except StatisticsError:
			median_length = None

		if len(regions) > 0:
			proportion_chromosome_annotated = round(sum(locus_lengths) / chrom_length, 4)
		else:
			proportion_chromosome_annotated = None

		# print(f"     length:")
		# print(f"       mean ---> {mean_length} bp")
		# print(f"       median -> {median_length} bp")
		# print(f"       proportion of chromosome annotated -> {proportion_chromosome_annotated}")

		read_depths = [sum(strand_d[l[0]].values()) for l in regions]
		try:
			mean_depth   = round(mean(read_depths),1)
		except StatisticsError:
			mean_depth = None
		try:
			median_depth = median(read_depths)
		except StatisticsError:
			median_depth = None
		if len(regions) > 0:
			proportion_libraries_annotated = round(sum(read_depths) / chrom_depth_c[chrom], 4)
		else:
			proportion_libraries_annotated = None

		# if mean_depth:
		# 	print(f"     abundance:")
		# 	print(f"       mean ---> {mean_depth} reads ({round(mean_depth * read_equivalent, 2)} rpm)")
		# 	print(f"       median -> {median_depth} reads ({round(median_depth * read_equivalent, 2)} rpm)")
		# 	print(f"       proportion of reads annotated -> {proportion_libraries_annotated}")
			
		# sys.exit()

		# • Formalize the coverage file outputs
		# • Formalize all other intermediate file outputs
		# • Include a log file (any tables?)
		# • Produce an alignment process which extracts tables of log.

		with open(stats_file, 'a') as outf:
			out = [chrom]
			out += [unclumped_regions_count, len(regions)]
			out += [chrom_length, proportion_chromosome_annotated, mean_length, median_length]
			out += [chrom_depth_c[chrom], proportion_libraries_annotated, mean_depth, median_depth]
			print("\t".join(map(str, out)), file=outf)

		# print()


		try:
			overall_d['region_count']  += unclumped_regions_count
			overall_d['regions_count']    += len(regions)
			overall_d['genome_length'] += chrom_length
			overall_d['locus_lengths'] += locus_lengths
			overall_d['total_depth']   = aligned_read_count
			overall_d['read_depths']   += read_depths

		except KeyError:
			overall_d['region_count']   = unclumped_regions_count
			overall_d['regions_count']     = len(regions)
			overall_d['genome_length']  = chrom_length
			overall_d['locus_lengths']  = locus_lengths
			overall_d['total_depth']    = aligned_read_count
			overall_d['read_depths']    = read_depths


		print()



	print()
	# print(f"{total_region_space:,} bp ({round(total_region_space/genome_length*100,1)}%) total region space ")
	# print(f"{total_annotated_reads:,} ({round(total_annotated_reads/aligned_read_count *100,1)}%) total reads in regions")
	print()
	print()


	with open(overall_file, 'a') as outf:


		print('project\tannotation_name\tregion_count\tregions_count\tgenome_length\tproportion_genome_annotated\tmean_length\tmedian_length\ttotal_depth\tproportion_library_annotated\tmean_depth\tmedian_depth', file=outf)

		line = [
			project_name,
			annotation_name,
			overall_d['region_count'], 
			overall_d['regions_count'], 
			overall_d['genome_length']
		]

		if overall_d['regions_count'] == 0:
			line += ['NA', "NA", 'NA']
		else:
			line += [
				round(sum(overall_d['locus_lengths'])/overall_d['genome_length'], 4),
				round(mean(overall_d['locus_lengths']),1),
				median(overall_d['locus_lengths'])
			]

		line += [
			overall_d['total_depth']
		]

		if overall_d['regions_count'] == 0:
			line += ['NA', "NA", 'NA']
		else:
			line += [
				round(sum(overall_d['read_depths'])/overall_d['total_depth'], 4),
				round(mean(overall_d['read_depths']),1),
				median(overall_d['read_depths'])
			]

		print("\t".join(map(str, line)), file=outf)


		# print('region_count ..................', overall_d['region_count'], file=outf)
		# print('regions_count ....................', overall_d['regions_count'], file=outf)
		# print()

		# print('genome_length: ................', overall_d['genome_length'], 'bp', file=outf)
		# print('proportion_genome_annotated ...', 
		# 	round(sum(overall_d['locus_lengths'])/overall_d['genome_length'], 4), file=outf)
		# print('mean_length ...................', round(mean(overall_d['locus_lengths']),1), 'bp', file=outf)
		# print('median_length .................', median(overall_d['locus_lengths']), 'bp', file=outf)
		# print()


		# print('total_depth ...................', overall_d['total_depth'], 'reads', file=outf)
		# print('proportion_library_annotated ..', 
		# 	round(sum(overall_d['read_depths'])/overall_d['total_depth'], 4), file=outf)
		# print('mean_depth ....................', round(mean(overall_d['read_depths']),1), 'reads', file=outf)
		# print('median_depth ..................', median(overall_d['read_depths']), 'reads', file=outf)
	# print("converting coverage files to bigwig...")

	# depth_wig.convert(output_directory=output_directory)
	# peak_wig.convert(output_directory=output_directory)

	end_time = datetime.now()


	date_time = end_time.strftime("%Y/%m/%d, %H:%M:%S")
	elapsed = end_time - start_time
	print(f"Run completed: {date_time}  ({elapsed} elapsed)")	















