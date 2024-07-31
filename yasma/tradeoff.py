

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

# import numpy as np
# from statistics import quantiles
import math

from .generics import *
from .cli import cli

from statistics import mean, median, StatisticsError, stdev

from time import time

import re

# from tqdm import tqdm


from datetime import datetime


import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
# from pympler import asizeof


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
		minmax=(15,30)):


		self.min_size=minmax[0]
		self.max_size=minmax[1]

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

		####################################
		## revision Jul 12 2024
		### these used to just require majority for all... (> 0.5), but that is a really weak standard. For example, to have contiguous sizes make up just a bare majority?
		### I think it is more reasonable to say that it is -> freq(n.sizes) > n.sizes / (n.sizes+1)
		### depth safe guards are still important... otherwise we will get some weird loci
		### size_1 did not have a depth threshold... i have added it to 15 to make sure we're not calling loci with virtually no reads to be selective. (8/15 reads must be one size)
		####################################
		## revision Jul 18 2024
		## On second thought, i have opted for just a bare majority to consider a locus size specific.
		## This increasing threshold looks to hold many loci ~just outside~ of consideration. This makes some sense, where there is probably a single predominant size, and adding in peripheral off-sized reads is unlikely to add 1/6 (1-size) or 1/4 (2-sizes) of total locus abundance. 
		## I also upped the minimums abundances a bit. Seems like high-duplication loci (loci skewed towards one or a few sequences) are a problem and maybe this can help it.	
		####################################

		if self.size_1_depth > self.depth * 0.5 and self.depth > 30:
			sizecall = self.size_1_key

		elif self.size_2_depth > self.depth * 0.5 and self.depth > 45:
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
		self.header += ['Gap', 'skew', 'size_1n','size_1n_depth', 'size_2n','size_2n_depth', 'size_3n','size_3n_depth', 'sizecall']



	def format(self, locus, seq_c, strand_c, sizecall, aligned_depth, last_stop):

		name, chrom, start, stop = locus


		### Basic information

		depth = sum(seq_c.values())
		rpm = depth / aligned_depth * 1000000


		### ShortStack standard metrics

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


		# complexity = unique_reads / depth


		### More derived metrics


		complexity = unique_reads / (stop - start)

		skew = major_rna_depth / depth

		gap = start - last_stop


		



		frac_top   = round(frac_top,3)
		complexity = round(complexity,3)
		rpm        = round(rpm,3)
		skew       = round(skew, 3)

		sizecall.get()

		result_line = [f"{chrom}:{start}-{stop}", name, stop-start, depth, rpm]
		result_line += [unique_reads, frac_top, strand, major_rna, major_rna_depth, complexity]
		result_line += [
			gap, skew, 
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
			f'ID={name};sizecall={sizecall};depth={depth};rpm={rpm};fracTop={frac_top};complexity={complexity};skew={skew};majorRNA={major_rna}'
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

	ker_window = params['kernel_window']
	half_ker_window = math.floor(ker_window/2)

	# kde_window = params['kernel_window']
	# half_kde_window = int(kde_window/2)

	# pos_c    = dict()
	pos_depth_d = dict()
	pos_strand_d = dict()
	pos_size_d = dict()

	genome_length = sum([c[1] for c in chromosomes])


	print(" processing alignment...")

	for chrom, chrom_length in chromosomes:
		pos_depth_d[chrom]  = [0] * chrom_length
		pos_strand_d[chrom] = [None] * chrom_length
		pos_size_d[chrom]   = [None] * chrom_length





	# ec = elapsedClass()
	iterables = []
	for c,l in chromosomes:
		iterables.append(samtools_view(bam, rgs=rgs, contig=c))#, read_minmax=(params['min_read_length'], params['max_read_length'])))

	reads = chain.from_iterable(iterables)

	print(f"    encoding reads ............... 0%", end='\r', flush=True)
	perc = percentageClass(1, sum(chrom_depth_c.values()))
	perc.update()

	aligned_read_count = 0
	for i, read in enumerate(reads):
		aligned_read_count+=1
		strand, length, _, pos, chrom, _, _, _ = read

		pos += math.floor(length / 2)


		pos_depth_d[chrom][pos] += 1

		try:
			pos_strand_d[chrom][pos][int(strand=="+")] += 1
			pos_size_d[chrom][pos].append(length)

		except TypeError:
			pos_strand_d[chrom][pos] = [0,0]
			pos_strand_d[chrom][pos][int(strand=="+")] += 1
			pos_size_d[chrom][pos] = [length]

		# try:
		# 	pos_d[chrom][pos][0] += 1
		# 	pos_d[chrom][pos][1].append(strand)
		# 	pos_d[chrom][pos][2].append(length)
		# except KeyError:
		# 	pos_d[chrom][pos] = [1, [strand],[length]]




		perc_out = perc.update()
		if perc_out:
			sys.stdout.write(f"    encoding reads ............... {perc_out}%\t {i+1:,} reads   \n", terminal_only=True)
			sys.stdout.flush()
			sys.stdout.overwrite_lines(1)

		# bin_c[int((chrom_index[chrom] + pos) / params['kernel_window'])] += 1

	print(f"    encoding reads ............... {perc.last_percent}%\t {i+1:,} reads   ", end='\n', flush=True)
	# print(ec)


	class trackClass():
		def __init__(self, bw_file, chromosomes):
			self.bw = pyBigWig.open(str(bw_file), 'w')
			self.bw.addHeader(chromosomes)

			self.last_start      = 0
			self.interval_length = 1
			self.last_chrom = chromosomes[0][0]
			self.last_val   = 0

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


	def numpy_method():
		ec = elapsedClass()

		gen_c = Counter() # a counter of kernel depths at positions in the genome 
		## sum(gen_c) = genome_length
		read_c = Counter() # a counter of the max kernel depth across the span of a read.



		cov_track = trackClass(Path(out_dir, "coverage.bw"), chromosomes)
		ker_track = trackClass(Path(out_dir, "kernel.bw"), chromosomes)


		for chrom_i, chromosome_entry in enumerate(chromosomes):
			chrom, chrom_length = chromosome_entry

			sys.stdout.write(f"    computing coverage ........... {chrom_i+1}/{len(chromosomes)} {chrom}                \n", terminal_only=True)
			sys.stdout.flush()
			sys.stdout.overwrite_lines(1)

			coverage = np.sum(sliding_window_view(np.asarray(pos_depth_d[chrom]), cov_window), -1)
			coverage = np.concatenate((np.array([coverage[0]] * half_cov_window), coverage, np.array([coverage[-1]] * (half_cov_window-1))), axis=0)

			kernel  = np.max(sliding_window_view(coverage, ker_window), -1)
			kernel = np.concatenate((np.array([kernel[0]] * half_ker_window), kernel, np.array([kernel[-1]] * (half_ker_window-1))), axis=0)

			for i, c in enumerate(coverage):
				cov_track.add(chrom, i+1, c)

			for i, k in enumerate(kernel):
				ker_track.add(chrom, i+1, k)

			gen_c.update(kernel)

			for i, d in enumerate(pos_depth_d[chrom]):
				if d > 0:
					read_c[kernel[i]] += d

		print(f"    computing coverage ........... {chrom_i+1}/{len(chromosomes)} {chrom}                ", end='\n', flush=True)

		print()
		cov_track.close()
		ker_track.close()

		print()
		print(ec)
		print()

		return(gen_c, read_c)

	# gen_c[rolling_max] += 1
	
	# read_c[rolling_max] += pos_depth_d[chrom][kernel_pos]

	gen_c, read_c = numpy_method()




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



		gen_c = Counter() # a counter of kernel depths at positions in the genome 
		## sum(gen_c) = genome_length
		read_c = Counter() # a counter of the max kernel depth across the span of a read.
		## sum(read_c) = total_aligned readcount



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
					depth = pos_depth_d[chrom][pos]
				except IndexError:
					print(pos, "<- index error 1")
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
					
					read_c[rolling_max] += pos_depth_d[chrom][kernel_pos]
					



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



	def tally(gen_c, read_c):
		found_depths = list(gen_c.keys())
		found_depths.sort()

		averages = []
		table    = []

		total_genomic_space  = genome_length
		total_possible_space = genome_length - gen_c[0]
		total_read_space     = sum(read_c.values())
		genp_thresholds      = []
		adj_genp_thresholds  = []
		readp_thresholds     = []

		for depth_threshold in found_depths:

			total_genomic_space -= gen_c[depth_threshold]
			total_read_space    -= read_c[depth_threshold]



			# print(depth_threshold, total_genomic_space, sep='\t')

			gen_score     = total_genomic_space / genome_length
			adj_gen_score = total_genomic_space / total_possible_space
			read_score    = total_read_space / aligned_read_count

			genp_thresholds.append((gen_score, depth_threshold))
			adj_genp_thresholds.append((adj_gen_score, depth_threshold))
			readp_thresholds.append((read_score, depth_threshold))

			## unweighted avg
			avg_score = ((1-gen_score) + read_score) / 2


			## weight avg
			weight_score = ((1-gen_score) * params['genome_weight'] + read_score * params['read_weight']) / sum([params['genome_weight'], params['read_weight']])



			geom_score = math.sqrt(((1-gen_score) * read_score))


			averages.append(round(weight_score, params['tradeoff_round']))

			table.append([depth_threshold, total_genomic_space, round(gen_score,4), round(adj_gen_score, 4),
				total_read_space, round(total_read_space/aligned_read_count,4),
				round(avg_score, 4), round(geom_score, 4), round(weight_score, 4)])

		peak_index = averages.index(max(averages))

		with open(Path(out_dir, 'thresholds.txt'), 'w') as outf:
			print('depth\tannotated_space\tp_genome\tadj_p_genome\tannotated_reads\tp_reads\taverage_score\tgeom_score\tweighted_avg\tpeak', file=outf)
			for i,t in enumerate(table):

				if i == peak_index:
					peak = 1
					out = {
						'depth_threshold' : t[0],
						'gen_score' : t[2],
						'adj_gen_score' : t[3],
						'read_score' : t[5],
						'avg_score' : t[6],
						'weighted_avg' : t[8]
						}
				else:
					peak = 0

				print("\t".join(map(str, t)), peak, sep='\t', file=outf)
				if total_genomic_space > genome_length:
					sys.exit("problem!!")

		# pprint(out)
		return(out, readp_thresholds, genp_thresholds)
	
	# print()
	# print(f'   pos_d: {asizeof.asizeof(pos_d):,}')
	# print(f' depth_c: {asizeof.asizeof(depth_c):,}')
	# print(f'kernel_c: {asizeof.asizeof(kernel_c):,}')



	# gen_c, read_c = inline_method()
	# tally_new(gen_c, read_c)
	out, readp_thresholds, genp_thresholds = tally(gen_c, read_c)
	# print(out)


	# gen_c, read_c = counter_method()
	# out = tally(gen_c, read_c)
	# print(out)
	# sys.exit()



	return(pos_depth_d, pos_strand_d, pos_size_d, out, readp_thresholds, genp_thresholds)






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
	default=['all'],
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
	help='This is the bandwidth for accumulating read alignments into coverage, which is used instead of the normal read length. By default, this is very large (250 nt), basically meaning that depth summed across 250 nt windows are used for region annotation.')

@optgroup.option('--kernel_window',
	default=250,
	help="This is a max filter for the coverage, which extends coverages by a default 250 nt. This is built-in padding for regions, which will then be revised to find boundaries.")



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
	default=1,
	help="The maximum number of subset alignments that will be written to the disk. Numbers higher than 1 are really only useful for performance comparisons. This value will automatically be raised to a minimum of the subsample_n+1.")





@optgroup.group('\n  Peak finding options',
				help='')

@optgroup.option("--genome_weight",
	default=1,
	help=f"along with --read_weight, these determine the weighted averages for considering the tradeoff proportion of reads and genome annotated. By default, this is weighted 2 for pReads and 1 for pGenome, meaning that the annotator tries do place more reads at the expense of more genome annotated. Default 1.")

@optgroup.option("--read_weight",
	default=1,
	help=f"Default 2. See above.")

# @optgroup.option("--tradeoff_weight",
# 	default = 0.65,
# 	help=f'Weighting factor applied to tradeoff averages. Higher specificity > 0.5 > higher sensitivity. Basically, the higher the value the more reads will be incorporated into annotations and resultingly more of the genome will be considered part of a locus. Default 0.65 (incorporating reads is 2x more important than being selective with the genome).')

@optgroup.option('--target_genome_perc',
	type=float,
	default=False,
	help='')
@optgroup.option('--target_read_perc',
	type=float,
	default=False,
	help='')

# @optgroup.option('--trim_regions', is_flag=True, default=False, help='Flag to include trimming of regions')


@optgroup.option('--tradeoff_round', default=3, help='Significance rounding for tradeoff average. Defaults to 3 digits (e.g. 0.977 or 97.7%)')





@optgroup.group('\n  Merging options',
				help='')

@optgroup.option("--merge_dist",
	default=500,
	help="Distance in nucleotides for which sRNA peaks should be considered for 'clumping'. Clumped regions must have sufficient similarity in sRNA-size profile and strand-preference. Default 500 nt.")

@optgroup.option("--merge_strand_similarity",
	default=0.7,
	help="Similarity threshold of strand fraction for clumping two peaks. Difference in fraction must be smaller than threshold. Default 0.5.")

@optgroup.option("--min_locus_length",
	default=30,
	help="Minimum size for a locus to be included. Default 30 nt.")



@optgroup.group('\n  Read options',
				help='')

@optgroup.option("--min_read_length",
	default=15,
	help="An override filter to ignore aligned reads which are smaller than a min length in locus calculations.")

@optgroup.option("--max_read_length",
	default=30,
	help="The same as above, but with a max length.")



@optgroup.group('\n  Locus options',
				help='')

@optgroup.option("--max_skew",
	default=0.90,
	type=float,
	help="A filter to help remove loci which are skewed toward only one sequence in abundance. By default (0.95), if more than 1 in 20 reads for a locus are a single sequence, they are excluded from the annotation.")


@optgroup.option("--min_complexity",
	default=0.01,
	type=float,
	help="")




@optgroup.group('\n Other options',
				help='')

@optgroup.option('--force', is_flag=True, default=False, help='force resubsample')
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

	clump_dist              = params['merge_dist']
	clump_strand_similarity = params['merge_strand_similarity']
	min_locus_length        = params['min_locus_length']
	debug                   = params['debug']
	annotation_name         = params['name']
	target_depth            = params['subsample']
	seed                    = params['subsample_seed']

	read_minmax = (params['min_read_length'], params['max_read_length'])

	params['output_directory'] = output_directory
	params['alignment_file'] = alignment_file
	params['project_name'] = project_name

	if annotation_name:
		dir_name = f"tradeoff_{annotation_name}"
	else:
		dir_name = 'tradeoff'



	if params['target_genome_perc'] and params['target_read_perc']:
		sys.exit("ERROR: cannot specify target read AND genome percentages (one is dependent on the other)")



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

	chromosome_max_lengths = {}
	for c,l in chromosomes:
		chromosome_max_lengths[c] = l

	if params['test_mode']:
		# chromosomes = chromosomes[3:5]
		# chromosomes = chromosomes[20:30]
		# chromosomes = chromosomes[2:5]
		# chromosomes = chromosomes[:2]
		# chromosomes = chromosomes[:1]
		chromosomes = chromosomes[4:5]


	genome_length = sum([l for c,l in chromosomes])

	# print(alignment_file)
	annotation_readgroups = check_rgs(annotation_readgroups, bam_rgs)
	params['annotation_readgroups'] = annotation_readgroups

	print()
	print(" Annotation readgroups:", annotation_readgroups)

	# chromosomes = [c for c in chromosomes if c[0] == 'NC_037320.1']

	chrom_depth_c = get_global_depth(alignment_file, aggregate_by=['rg','chrom'])

	keys = list(chrom_depth_c.keys())
	for key in keys:
		if key[0] in annotation_readgroups:
			chrom_depth_c[key[1]] += chrom_depth_c[key]

		del chrom_depth_c[key]

	for key in list(chrom_depth_c.keys()):
		if key not in [c for c,l in chromosomes]:
			del chrom_depth_c[key]

	aligned_read_count = sum(chrom_depth_c.values())




	if params['subsample']:

		alignment_file = subsample(aligned_read_count, alignment_file, params)
		chrom_depth_c = get_global_depth(alignment_file, aggregate_by=['rg','chrom'])
		aligned_read_count = sum(chrom_depth_c.values())



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
			contig=chrom, start=start, stop=stop, 
			rgs=annotation_readgroups, 
			boundary_rule = 'tight',
			read_minmax=read_minmax)]


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


	stats_file = f"{output_directory}/{dir_name}/stats_by_ref.txt"
	with open(stats_file, 'w') as outf:
		print("project\tchromosome\tregions\tregions\tchromosome_length\tproportion_genome_annotated\tmean_length\tmedian_length\treads\tproportion_libraries_annotated\tmean_abundance\tmedian_abundance", file=outf)


	overall_file = f"{output_directory}/{dir_name}/stats.txt"
	with open(overall_file, 'w') as outf:
		outf.write('')

	overall_d = {}

	all_regions = []
	total_read_count = 0
	cl_i = 0


	read_equivalent = 1 / sum(chrom_depth_c.values()) * 1000000


	## iterating through chromosomes and reads

	print()
	print(f" {len(chromosomes)} chromosomes/scaffolds/contigs")
	print(f" {sum([l for c,l in chromosomes]):,} bp total")
	print(f" {aligned_read_count:,} reads")
	print()




	## Getting positional coverage accross the alignment


	start = time()


	pos_depth_d, pos_strand_d, pos_size_d, threshold_stats, readp_thresholds, genp_thresholds = get_kernel_coverage(
		bam=alignment_file, 
		rgs=annotation_readgroups, 
		params=params, 
		chrom_depth_c=chrom_depth_c,
		chromosomes=chromosomes,
		out_dir=Path(output_directory, dir_name))

	# sys.exit()



	if params['target_genome_perc']:
		for p, t in genp_thresholds:
			if p < params['target_genome_perc']:
				break

		depth_threshold = t
		gen_score       = p

		for p, t in readp_thresholds:
			if t == depth_threshold:
				read_score = p
				break


		print(" annotation parameters...")
		print(f"    depth threshold: ......... {depth_threshold} reads")
		print(f" -> set genome proportion: ... {gen_score}")
		print(f"    exp. read proportion: .... {read_score}")


	elif params['target_read_perc']:
		for p, t in readp_thresholds:
			if p < params['target_read_perc']:
				break

		depth_threshold = t
		read_score      = p

		for p, t in genp_thresholds:
			if t == depth_threshold:
				gen_score = p
				break


		print(" annotation parameters...")
		print(f"    depth threshold: ......... {depth_threshold} reads")
		print(f"    exp. genome proportion: .. {gen_score}")
		print(f" -> set read proportion: ..... {read_score}")

	else:
		print(f"Finding threshold through weighted tradeoff. Weight: [{params['read_weight']}] reads to [{params['genome_weight']}] genome")

		depth_threshold = threshold_stats['depth_threshold']
		gen_score       = threshold_stats['gen_score']
		adj_gen_score   = threshold_stats['adj_gen_score']
		read_score      = threshold_stats['read_score']

		print(" annotation parameters...")
		print(f"    depth threshold: ......... {depth_threshold} reads")
		print(f"    exp. genome proportion: .. {gen_score}")
		print(f"    exp. read proportion: .... {read_score}")




	def get_regions(depth_threshold, chromosomes):


		regions = []

		def check_and_cash_region(in_region, reg_i, chrom, reg_start, reg_stop, chrom_length):
			if in_region:

				if reg_stop >= chrom_length:
					reg_stop = chrom_length - 1

				reg_i += 1
				region = [f"region_{reg_i}", chrom, reg_start, reg_stop]
				regions.append(region)

			return(reg_i)

		reg_i = 0
		reg_start = -1
		reg_stop  = -1

		bw = pyBigWig.open(str(Path(output_directory, dir_name, "kernel.bw")))
		# print(bw.isBigWig())
		# print(bw.chroms())
		# print(bw.header())

		for chrom, chrom_length in chromosomes:

			in_region = False

			for inv_start, inv_stop, depth in bw.intervals(chrom, 0, chrom_length):
				if depth >= depth_threshold:
					if not in_region:
						reg_start = inv_start
					reg_stop  = inv_stop

					in_region = True

				else:


					reg_i = check_and_cash_region(in_region, reg_i, chrom, reg_start, reg_stop, chrom_length)
					in_region = False



			reg_i = check_and_cash_region(in_region, reg_i, chrom, reg_start, reg_stop, chrom_length)


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
	print(f"    {total_region_space:,} genomic nt ({round(total_region_space/genome_length *100,1)}%) in regions")
	print(f"        expected: {round(100*gen_score,1)}%")



	total_annotated_reads = 0
	for l in all_regions:
		# print(l)

		name, chrom, start, stop = l
		for r in range(start, stop+1):
			total_annotated_reads += pos_depth_d[chrom][r]



	print(f"    {total_annotated_reads:,} reads ({round(total_annotated_reads/aligned_read_count *100,1)}%) in regions")
	print(f"        expected: {round(100*read_score,1)}%")




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
		sizes   = sizeClass(minmax=read_minmax)
		names   = set()




		for p in range(start, stop+1):
			try:
				strands['+'] += pos_strand_d[chrom][p][0]
				strands['-'] += pos_strand_d[chrom][p][1]
				sizes.update(pos_size_d[chrom][p])
			except TypeError:
				pass
			except IndexError:
				print(p, "<- index error 2")
				pass


		try:
			frac_top = strands['+'] / sum(strands.values())
		except ZeroDivisionError:
			print(region)

			sys.exit()


		r_depth = sum(strands.values())

		# cov_window = params['coverage_window']

		region_size = stop - start

		# print()


		def window_gen(start, size, direction, increment, inset=False):

			if inset:
				start = start - size * direction

			if start < 0:
				start = 0

			if start >= chromosome_max_lengths[chrom]:
				start = chromosome_max_lengths[chrom]

			while True:
				end =  start + size * direction

				if start < 0:
					yield sorted([0,end])
					return

				if end < 0:
					yield sorted([0, start])
					return

				if start >= chromosome_max_lengths[chrom]:
					yield sorted([end, chromosome_max_lengths[chrom]-1])
					return

				if end >= chromosome_max_lengths[chrom]:
					yield sorted([start, chromosome_max_lengths[chrom]-1])
					return

				yield sorted([start,end])

				start += increment * direction



		def test_extend(window):

			window_start, window_end = window

			window = list(range(window_start, window_end)) 
			window_size = len(window)

			w_strands = Counter()
			w_sizes   = sizeClass(minmax=read_minmax)
			w_depths  = 0

			for w in window:
				# print(w)

				# print(w)
				try:
					depth = pos_depth_d[chrom][w]
				except IndexError:
					print(w)
					print(window_start, window_end)
					print(chromosome_max_lengths[chrom])
					sys.exit()

				w_depths += depth

				if depth > 0:
					w_strands["+"] += pos_strand_d[chrom][w][0]
					w_strands["-"] += pos_strand_d[chrom][w][1]
					w_sizes.update(pos_size_d[chrom][w])

				# print(' ', w, depth)

				if w in claim_d[chrom] and claim_d[chrom][w] != name:
					# print('claim break')
					return(False, ['claim'])

			try:
				w_fractop = w_strands['+'] / sum(w_strands.values())
			except ZeroDivisionError:
				w_fractop = 1


			expand_fails = []

			if w_depths == 0:
				# pprint(w_strands)
				# print('empty break')
				# sys.exit()

				expand_fails.append('empty')

			else:

				if window_size > 100:

					if not abs(frac_top - w_fractop) < 0.5:
						# print('strand break')
						expand_fails.append('strand')
					

					if not w_sizes == sizes:
						# print('size break')
						expand_fails.append('size')


				# print(r_depth, region_size, w_depths, window_size)
				# print(round(r_depth / region_size * window_size, 1), "read threshold")
				if w_depths/window_size < r_depth/region_size * 0.05:
					# print('depth break')
					expand_fails.append('depth')

			# print("pass")

			if len(expand_fails) > 0:
				return(False, expand_fails)

			return(True, [])




		def find_outer_boundaries(gen):

			last_window = next(gen)
			for window in gen:
				test, fail_list = test_extend(window)

				# print(" ", window, test, fail_list)

				if not test:
					return(last_window)

				last_window = window
			return last_window



		def find_inner_boundaries(gen):

			for window in gen:
				test, fail_list = test_extend(window)

				# print(" ", window, test, fail_list)

				if test:
					return(window)
			return window


		# if chrom == 'NC_003076.8':
		# 	if start > 3456817:
		# 		sys.exit()




		# print(region)
		# print("boundaries outward")
		# boundaries outward - coarse
		window_size = 250
		# print("  right outward")
		new_stop  = find_outer_boundaries(window_gen(start=stop,  size=window_size, direction=1,  increment=30, inset=True))[1]
		# print("  left outward")
		new_start = find_outer_boundaries(window_gen(start=start, size=window_size, direction=-1, increment=30, inset=True))[0]

		# print("boundaries inward")
		## boundaries inward - fine
		window_size = 50
		# print("  right inward")
		new_stop  = find_inner_boundaries(window_gen(start=new_stop,  size=window_size, direction=-1, increment=5))[1]
		# print("  left inward")
		new_start = find_inner_boundaries(window_gen(start=new_start, size=window_size, direction=1,  increment=5))[0]



		## cleaning up claims if locus shrank
		for p in range(start, new_start+1):
			claim_d[p] = None
		for p in range(new_stop, stop+1):
			claim_d[p] = None

		## adding claims if locus expanded
		for p in range(new_start, start+1):
			claim_d[p] = name
		for p in range(stop, new_stop+1):
			claim_d[p] = name

		return(new_start, new_stop, claim_d)


	claim_d = {}
	for c, l in chromosomes:
		claim_d[c] = {}

	for region in all_regions:
		name, chrom, start, stop = region
		for r in range(start, stop+1):
			claim_d[chrom][r] = name


	class reviseClass():
		def __init__(self, region):

			self.locus_name, self.chrom, self.start, self.stop = region

			self.coords  = f"{self.chrom}:{self.start}-{self.stop}"
			# self.claim_d = claim_d

			self.strands = Counter()
			self.sizes   = sizeClass(minmax=read_minmax)
			self.names   = set()


			for p in range(self.start, self.stop+1):
				try:
					self.strands['+'] += pos_strand_d[self.chrom][p][0]
					self.strands['-'] += pos_strand_d[self.chrom][p][1]
					self.sizes.update(pos_size_d[self.chrom][p])
				except TypeError:
					pass
				except IndexError:
					print(p, "<- index error 2")
					pass


			try:
				self.frac_top = self.strands['+'] / sum(self.strands.values())
			except ZeroDivisionError:
				print(self.region)

				sys.exit()


			self.r_depth = sum(self.strands.values())







		def expand(self):
			# print(region)
			# print("boundaries outward")
			# boundaries outward - coarse
			window_size = 250
			# print("  right outward")
			new_stop  = self.find_outer_boundaries(self.window_gen(start=self.stop,  size=window_size, direction=1,  increment=30, inset=True))[1]
			# print("  left outward")
			new_start = self.find_outer_boundaries(self.window_gen(start=self.start, size=window_size, direction=-1, increment=30, inset=True))[0]

			self.cleanup(new_start, new_stop)

			return(new_start, new_stop)

		def trim(self):
			# print("boundaries inward")
			## boundaries inward - fine
			window_size = 50
			# print("  right inward")
			new_stop  = self.find_inner_boundaries(self.window_gen(start=self.stop,  size=window_size, direction=-1, increment=5))[1]
			# print("  left inward")
			new_start = self.find_inner_boundaries(self.window_gen(start=self.start, size=window_size, direction=1,  increment=5))[0]

			if new_stop < self.start + window_size or new_start > self.stop - window_size:
				# print(f"  warning: trim error {self.locus_name}")
				return(self.start, self.stop)

			self.cleanup(new_start, new_stop)

			return(new_start, new_stop)


		def cleanup(self, new_start, new_stop):

			## cleaning up claims if locus shrank
			for p in range(self.start, new_start+1):
				claim_d[self.chrom][p] = None
			for p in range(new_stop, self.stop+1):
				claim_d[self.chrom][p] = None

			## adding claims if locus expanded
			for p in range(new_start, self.start+1):
				claim_d[self.chrom][p] = self.locus_name
			for p in range(self.stop, new_stop+1):
				claim_d[self.chrom][p] = self.locus_name






		def window_gen(self, start, size, direction, increment, inset=False):

			if inset:
				start = start - size * direction

			if start < 0:
				start = 0

			if start >= chromosome_max_lengths[chrom]:
				start = chromosome_max_lengths[chrom]

			while True:
				end =  start + size * direction

				if start < 0:
					yield sorted([0,end])
					return

				if end < 0:
					yield sorted([0, start])
					return

				if start >= chromosome_max_lengths[chrom]:
					yield sorted([end, chromosome_max_lengths[chrom]-1])
					return

				if end >= chromosome_max_lengths[chrom]:
					yield sorted([start, chromosome_max_lengths[chrom]-1])
					return

				yield sorted([start,end])

				start += increment * direction


		def test_extend(self, window, expand=False):

			window_start, window_end = window

			if window_start < self.start or window_end > self.stop:
				return(False, ['outofbounds'])


			window = list(range(window_start, window_end)) 
			window_size = len(window)

			w_strands = Counter()
			w_sizes   = sizeClass(minmax=read_minmax)
			w_depths  = 0

			for w in window:
				# print(w)

				# print(w)
				try:
					depth = pos_depth_d[chrom][w]
				except IndexError:
					print(w)
					print(window_start, window_end)
					print(chromosome_max_lengths[chrom])
					sys.exit()

				w_depths += depth

				if depth > 0:
					w_strands["+"] += pos_strand_d[chrom][w][0]
					w_strands["-"] += pos_strand_d[chrom][w][1]
					w_sizes.update(pos_size_d[chrom][w])

				# print(' ', w, depth)

				if expand: 
					## trimming doesnt worry about hitting other loci
					if w in claim_d[chrom] and claim_d[chrom][w] != name:
						# print('claim break')
						return(False, ['claim'])

			try:
				w_fractop = w_strands['+'] / sum(w_strands.values())
			except ZeroDivisionError:
				w_fractop = 1


			expand_fails = []

			if w_depths == 0:
				# pprint(w_strands)
				# print('empty break')
				# sys.exit()

				expand_fails.append('empty')

			else:

				if expand:
					## trimming doesnt consider fractop or size distributions (smaller windows make this challenging)

					if not abs(self.frac_top - w_fractop) < 0.5:
						# print('strand break')
						expand_fails.append('strand')
					

					if not w_sizes == self.sizes:
						# print('size break')
						expand_fails.append('size')


				# print(r_depth, region_size, w_depths, window_size)
				# print(round(r_depth / region_size * window_size, 1), "read threshold")
				if w_depths/window_size < self.r_depth/(self.stop - self.start) * 0.05:
					# print('depth break')
					expand_fails.append('depth')

			# print("pass")

			if len(expand_fails) > 0:
				return(False, expand_fails)

			return(True, [])




		def find_outer_boundaries(self, gen):

			last_window = next(gen)
			for window in gen:
				test, fail_list = self.test_extend(window, expand=True)

				# print(" ", window, test, fail_list)

				if not test:
					return(last_window)

				last_window = window
			return last_window



		def find_inner_boundaries(self, gen):

			for window in gen:
				test, fail_list = self.test_extend(window)


				# print(" ", window, test, fail_list)

				if 'outofbounds' in fail_list:
					print(self.start, self.stop)
					sys.exit()

				if test:
					return(window)
			return window


	## revising regions
	
	print()
	sys.stdout.write(f' revising regions ... 0%  \r', terminal_only=True)
	sys.stdout.flush()

	revised_genomic_space = 0
	perc = percentageClass(increment=5, total=len(all_regions))

	with open(revised_gff, 'a') as outf:

		for region in all_regions:
			name, chrom, start, stop = region

			# print(region)


			rc = reviseClass(region)
			new_start, new_stop = rc.expand()

			if new_stop > chromosome_max_lengths[chrom]:
				new_stop = chromosome_max_lengths[chrom] - 1

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

			# print(region)
			# print()


			gff_line = [
				chrom, 'test_region','region', new_start, new_stop, '.', '.', '.',
				f'ID={name}']

			print('\t'.join(map(str,gff_line)), file=outf)

			perc_out = perc.update()
			if perc_out:
				sys.stdout.write(f' revising regions ... {perc_out}%  \r', terminal_only=True)
				sys.stdout.flush()

		# if name == 'region_2':
		# 	sys.exit()

	print(f' revising regions ... {perc.last_percent}%   ', flush=True)


	total_revised_reads = 0
	for l in all_regions:
		# print(l)

		name, chrom, start, stop = l

		for r in range(start, stop+1):
			try:
				total_revised_reads += pos_depth_d[chrom][r]
			except IndexError:
				print(chrom, r, "<- index error 4")
				pass


	def string_plus_white(s, length = 7):
		s = str(s)
		return s + " " * (length - len(s))



	print(f"    {revised_genomic_space:,} genomic nt ({round(revised_genomic_space/genome_length *100,1)}%) in revised regions")
	print(f"    {total_revised_reads:,} reads ({round(total_revised_reads/aligned_read_count *100,1)}%) in revised regions")
	diff = round(total_revised_reads/aligned_read_count *100,1) - round(total_annotated_reads/aligned_read_count *100,1)
	print(f"      +{round(diff,2)}% over unrevised regions")

	max_chrom_word_length = max([len(c) for c,l in chromosomes])
	def print_progress_string(i, n, chrom, input_loci, output_loci, assess=False, terminal_only=False):

		chrom = chrom + (max_chrom_word_length - len(chrom)) * " "

		if not assess:
			assess = ''
		else:
			assess = f"\t{assess}%"

		sys.stdout.write(f"{i+1}/{n}\t{chrom}\t{string_plus_white(unclumped_regions_count)}\t{string_plus_white(len(regions))}{assess}  \r", 
			terminal_only=terminal_only)
		sys.stdout.flush()

	print()
	print('prog', "chrom"+(max_chrom_word_length-5)*" ",'regions\tloci\tassess', sep='\t')

	total_region_space = 0
	regions_name_i = 0
	total_annotated_reads = 0
	complexity_filter = 0
	skew_filter       = 0


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




		locus_i = 0
		considered_regions = get_considered_regions(locus_i)
		unclumped_regions_count = len(regions)


		print_progress_string(chrom_count, len(chromosomes), chrom, unclumped_regions_count, len(regions), terminal_only=True)

		last_claim = 'start'

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








		for read in samtools_view(alignment_file, contig=chrom, rgs=annotation_readgroups):#, read_minmax=read_minmax):

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
				sizecall_d[claim] = sizeClass(minmax=read_minmax)
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
							print(f"Warning: bad locus {c_name} removed                      ")


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
								current_sc = sizeClass(minmax=read_minmax)
								current_ft = -1
								# print(f"\nWarning: region {current_locus} may not have counts")


							try:
								next_sc = sizecall_d[next_locus[0]]
								next_ft = get_frac_top(strand_d[next_locus[0]])
							except KeyError:
								next_sc = sizeClass(minmax=read_minmax)
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

				print_progress_string(chrom_count, len(chromosomes), chrom, unclumped_regions_count, len(regions), terminal_only=True)


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

		# if params['trim_regions']:
		# 	for i,locus in enumerate(regions):
		# 		locus[2:] = trim_locus(locus, peak_c[locus[1]], params)





		## Assessing locus dimensions and making annotations

		perc = percentageClass(increment=5, total=len(regions))


		last_stop = 0
		for i,locus in enumerate(regions):

			# print()
			# print(locus)

			rc = reviseClass(locus)
			locus[2], locus[3] = rc.trim()
			# print(locus)

			old_name = locus[0]
			regions_name_i += 1
			locus[0] = f"locus_{regions_name_i}"


			print_percentage = perc.get_percent(i)
			if print_percentage:
				print_progress_string(chrom_count, len(chromosomes), chrom, unclumped_regions_count, len(regions), print_percentage, terminal_only=True)

			name, chrom, start, stop = locus
			coords = f"{chrom}:{start}-{stop}"

			# print(locus)

			strand_d[name]   = strand_d.pop(old_name)
			sizecall_d[name] = sizecall_d.pop(old_name)
			seq_d[name]      = seq_d.pop(old_name)

			strand_c = strand_d[name]
			sizecall = sizecall_d[name]
			read_c   = seq_d[name]

			complexity = len(read_c.keys()) / (stop - start)
			skew       = read_c.most_common(1)[0][1] / sum(read_c.values())


			if complexity < params['min_complexity']:
				complexity_filter += 1

			if skew > params['max_skew']:
				skew_filter += 1

			# if complexity < params['min_complexity'] or skew > params['max_skew']:
			# 	regions_name_i -= 1
			# 	continue


			results_line, gff_line = assessClass().format(locus, read_c, strand_c, sizecall, sum(chrom_depth_c.values()), last_stop)


			last_stop = stop

			with open(results_file, 'a') as outf:
				print("\t".join(map(str, results_line)), file=outf)

			with open(gff_file, 'a') as outf:
				print("\t".join(map(str, gff_line)), file=outf)

			top_reads_save(read_c, reads_file, read_equivalent, name)


		if regions:
			print_progress_string(chrom_count, len(chromosomes), chrom, unclumped_regions_count, len(regions), print_percentage, terminal_only=False)
		else:
			print_progress_string(chrom_count, len(chromosomes), chrom, unclumped_regions_count, len(regions), terminal_only=False)


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
			try:
				proportion_libraries_annotated = round(sum(read_depths) / chrom_depth_c[chrom], 4)
			except ZeroDivisionError:
				proportion_libraries_annotated = None
		else:
			proportion_libraries_annotated = None



		with open(stats_file, 'a') as outf:
			out = [project_name, chrom]
			out += [unclumped_regions_count, len(regions)]
			out += [chrom_length, proportion_chromosome_annotated, mean_length, median_length]
			out += [chrom_depth_c[chrom], proportion_libraries_annotated, mean_depth, median_depth]
			print("\t".join(map(str, out)), file=outf)

		# print()


		try:
			overall_d['region_count']  += unclumped_regions_count
			overall_d['regions_count'] += len(regions)
			overall_d['genome_length'] += chrom_length
			overall_d['locus_lengths'] += locus_lengths
			overall_d['total_depth']   = aligned_read_count
			overall_d['read_depths']   += read_depths

		except KeyError:
			overall_d['region_count']   = unclumped_regions_count
			overall_d['regions_count']  = len(regions)
			overall_d['genome_length']  = chrom_length
			overall_d['locus_lengths']  = locus_lengths
			overall_d['total_depth']    = aligned_read_count
			overall_d['read_depths']    = read_depths


		print()


	print()
	print(f"{complexity_filter} loci have low complexity (< {params['min_complexity']})")
	print(f"{skew_filter} loci are extremely skewed (> {params['max_skew']})")
	print()
	print()


	with open(overall_file, 'a') as outf:


		print('project\tannotation_name\tregion_count\tlocus_count\tgenome_length\tproportion_genome_annotated\tmean_length\tmedian_length\ttotal_depth\tproportion_library_annotated\tmean_depth\tmedian_depth', file=outf)

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


	end_time = datetime.now()


	date_time = end_time.strftime("%Y/%m/%d, %H:%M:%S")
	elapsed = end_time - start_time
	print(f"Run completed: {date_time}  ({elapsed} elapsed)")	















