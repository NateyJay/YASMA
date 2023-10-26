

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
		# 	# this solution is a stop-gap. Apparently, in some species low-lambdas can lead to trace-loci being annotated, which only have reads outside of accepted ranges. 
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

		gff_line = [
			chrom, 'yasma_peak',feature_type, start, stop, '.', strand, '.',
			f'ID={name};sizecall={sizecall};depth={depth};rpm={rpm};fracTop={frac_top};majorRNA={major_rna}'
		]


		return(result_line, gff_line)






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

@optgroup.option('--kernel_window',
	default=40,
	help='This is the bandwidth smoothing window for the kernel coverage method. Default 40 nt.')

@optgroup.option('--coverage_method',
	default='kernel',
	type=click.Choice(['kernel', 'depth']),
	help="Choice between two methods for finding genomic coverage. 'Kernel' uses a density smoothing based on the center of the read. 'depth' uses samtools depth and can be sensitive to read-length.")


@optgroup.option('--subsample',
	help="Allows the user to subsample alignments for the annotation to a defined depth. Accepts an integer number of reads, which can be modified with a 10^3 prefix (ex. 10M).")




@optgroup.group('\n  Peak finding options',
                help='')

@optgroup.option("--creep_dist",
	default=50,
	help="Distance used in 'creep' method for boundary finding. Default 50 nt.")

@optgroup.option('--peak_threshold',
	default=0.05,
	help='Minimum depth, as a proportion of the maximum coverage depth in a locus, used for trimming low-depth edges from a locus. Default 0.05 proportion of peak.')

@optgroup.option("--rpm_threshold",
	default=0.5,
	help="Depth threshold in reads per million for discovery of a locus peak. Default 0.5 RPM.")

@optgroup.option('--pad',
	default=10,
	help='Number of bases arbitrarily added to either end of a defined locus. Default 10 nt.')

@optgroup.option('--boundary_method',
	default='creep',
	type=click.Choice(['creep', 'bool']),
	help="Choice between two methods for finding peak boundaries. 'creep' expands peaks based on a creeping window and makes wider peaks with the creep_dist option. 'bool' makes more concise peaks based on the pad distance.")





@optgroup.group('\n  Merging options',
                help='')


@optgroup.option("--clump_dist",
	default=500,
	help="Distance in nucleotides for which sRNA peaks should be considered for 'clumping'. Clumped loci must have sufficient similarity in sRNA-size profile and strand-preference. Default 500 nt.")

@optgroup.option("--clump_strand_similarity",
	default=0.7,
	help="Similarity threshold of strand fraction for clumping two peaks. Difference in fraction must be smaller than threshold. Default 0.5.")

@optgroup.option("--min_locus_length",
	default=30,
	help="Minimum size for a locus to be included. Default 30 nt.")



@optgroup.group('\n Other options',
                help='')

@optgroup.option('--debug', is_flag=True, default=False, help='Debug flag')
@optgroup.option('--override', is_flag=True, default=False, help='Overrides config file changes without prompting.')


def peak(**params):
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
	creep_dist              = params['creep_dist']
	peak_threshold          = params['peak_threshold']
	rpm_threshold           = params['rpm_threshold']
	pad                     = params['pad']
	boundary_method         = params['boundary_method']
	clump_dist              = params['clump_dist']
	clump_strand_similarity = params['clump_strand_similarity']
	min_locus_length        = params['min_locus_length']
	debug                   = params['debug']
	annotation_name         = params['name']
	target_depth            = params['subsample']

	if annotation_name:
		dir_name = f"peak_{annotation_name}"
	else:
		dir_name = 'peak'




	if debug: 
		show_warnings = True
	else:
		show_warnings = False


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


	if target_depth:
		subsample = parse_subsample(target_depth, alignment_file, "bam", sum(chrom_depth_c.values()))

		perform_subsample(subsample)

		alignment_file = subsample_file



	def test_by_samtools(c):
		coords = f"{c[1]}:{c[2]}-{c[3]}"

		sc = sizeClass()
		strand_c = Counter()
		sizes, strands = [], []

		for read in samtools_view(alignment_file, locus=coords, rgs=annotation_readgroups, boundary_rule='tight'):
			sam_strand, sam_length, _, sam_lbound, _, _, _, read_name = read
			sam_rbound = sam_lbound + sam_length

			sc.update([sam_length])
			strand_c[sam_strand] += 1

			sizes.append(sam_length)
			strands.append(sam_strand)


		# print()
		# print("sam read count:", sum(strand_c.values()))

		ft = get_frac_top(strand_c)

		return(sc, ft, sizes, strands)


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

	gff_file = f"{output_directory}/{dir_name}/loci.gff3"

	with open(gff_file, 'w') as outf:
		print("##gff-version 3", file=outf)

		for chrom, chrom_length in chromosomes:
			print(f"##sequence-region   {chrom} 1 {chrom_length}", file=outf)


	results_file = f"{output_directory}/{dir_name}/loci.txt"
	with open(results_file, 'w') as outf:
		print("\t".join(assessClass().header), file=outf)


	reads_file = f"{output_directory}/{dir_name}/reads.txt"
	with open(reads_file, 'w') as outf:
		print(TOP_READS_HEADER, file=outf)



	region_file = f"{output_directory}/{dir_name}/regions.gff3"
	with open(region_file, 'w') as outf:
		print("##gff-version 3", file=outf)

		for chrom, chrom_length in chromosomes:
			print(f"##sequence-region   {chrom} 1 {chrom_length}", file=outf)

	merge_file = f"{output_directory}/{dir_name}/merges.txt"
	with open(merge_file, 'w') as outf:
		outf.write('')


	stats_file = f"{output_directory}/{dir_name}/stats_by_chrom.txt"
	with open(stats_file, 'w') as outf:
		print("chromosome\tregions\tloci\tchromosome_length\tproportion_genome_annotated\tmean_length\tmedian_length\treads\tproportion_libraries_annotated\tmean_abundance\tmedian_abundance", file=outf)


	overall_file = f"{output_directory}/{dir_name}/stats.txt"
	with open(overall_file, 'w') as outf:
		outf.write('')

	overall_d = {}

	all_loci = []
	total_read_count = 0
	cl_i = 0


	# total_reads = sum(chrom_depth_c.values())
	read_equivalent = 1 / sum(chrom_depth_c.values()) * 1000000


	## iterating through chromosomes and reads

	print(f" {len(chromosomes)} chromosomes/scaffolds/contigs")
	print(f" {sum([l for c,l in chromosomes]):,} bp total")
	print(f" {sum(chrom_depth_c.values()):,} reads")
	print()

	loci = []

	def get_peaks(bam, rgs):
		""" Produces coverages for alignment using samtools depth. 
		Returns a dict() of Counter() objects, where the dict keys are chromosomes, counter keys are chromosome positions, and counter values are depth at positions.
		"""

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



	# def dep_get_kernel_peaks(bam, rgs):
	# 	""" Produces coverages for alignment based off a kernel density estimation. 

	# 	Sum of all reads within the user defined kernel bandwith are used to generate a coverage. This is meant to normalize the effect of read length on coverage.

	# 	Returns a dict() of Counter() objects, where the dict keys are chromosomes, counter keys are chromosome positions, and counter values are depth at positions.
	# 	"""


	# 	## Counting read depth by position.

	# 	half_window = math.floor(kernel_window/2)

	# 	kernel_c = {}
	# 	depth_c = dict()

	# 	perc = percentageClass(1, sum(chrom_depth_c.values()))


	# 	for chrom, chrom_length in chromosomes:
	# 		depth_c[chrom] = Counter()
	# 		kernel_c[chrom] = Counter()


	# 	reads = samtools_view(bam, rgs=rgs)

	# 	print(f"   reading position depths ..... 0%", end='\r', flush=True)

	# 	i = 0
	# 	for i, read in enumerate(reads):
	# 		_, length, _, pos, chrom, _, _, _ = read

	# 		pos += math.floor(length / 2)

	# 		depth_c[chrom][pos-half_window] += 1

	# 		perc_out = perc.get_percent(i)
	# 		if perc_out:
	# 			print(f"   reading position depths ..... {perc_out}%\t{i:,} reads", end='\r', flush=True)


	# 	print()


	# 	## Performing kernel density smoothing based on mean or median or --> sum <--?

	# 	window_dq = deque()

	# 	perc = percentageClass(1, sum([l for c,l in chromosomes]))
	# 	ii =0

	# 	for chrom, chrom_length in chromosomes:
	# 		kernel_c[chrom] = Counter()


	# 		for i in range(chrom_length):
	# 			ii += 1

	# 			window_dq.append(depth_c[chrom][i])

	# 			summary_value = round(sum(window_dq),2)

	# 			if summary_value > 0:
	# 				kernel_c[chrom][i] = summary_value
				



	# 			perc_out = perc.get_percent(ii)
	# 			if perc_out:
	# 				print(f"   calculating kernel density .. {perc_out}% \t{round(ii/1000000, 1)}M nt ", end='\r', flush=True)


	# 			if len(window_dq) > kernel_window:
	# 				window_dq.popleft()


	# 	print()

	# 	return(kernel_c)



	def get_kernel_peaks(bam, rgs):
		""" Produces coverages for alignment based off a kernel density estimation. 

		Sum of all reads within the user defined kernel bandwith are used to generate a coverage. This is meant to normalize the effect of read length on coverage.

		Returns a dict() of Counter() objects, where the dict keys are chromosomes, counter keys are chromosome positions, and counter values are depth at positions.
		"""


		## Counting read depth by position.

		half_window = math.floor(kernel_window/2)

		kernel_c = {}
		depth_c = dict()

		perc = percentageClass(1, sum(chrom_depth_c.values()))


		for chrom, chrom_length in chromosomes:
			depth_c[chrom] = Counter()
			kernel_c[chrom] = Counter()


		reads = samtools_view(bam, rgs=rgs)

		print(f"   calculating kernel density .. 0%", end='\r', flush=True)

		i = 0
		for i, read in enumerate(reads):
			_, length, _, pos, chrom, _, _, _ = read

			pos += math.floor(length / 2)

			# for r in range(pos-half_window, pos+half_window+1):
			# 	depth_c[chrom][r] += 1

			depth_c[chrom].update(list(range(pos-half_window, pos+half_window+1)))

			perc_out = perc.get_percent(i)
			if perc_out:
				print(f"   calculating kernel density .. {perc_out}%\t{i:,} reads", end='\r', flush=True)


		print()


		# ## Performing kernel density smoothing based on mean or median or --> sum <--?

		# window_dq = deque()

		# perc = percentageClass(1, sum([l for c,l in chromosomes]))
		# ii =0

		# for chrom, chrom_length in chromosomes:
		# 	kernel_c[chrom] = Counter()


		# 	for i in range(chrom_length):
		# 		ii += 1

		# 		window_dq.append(depth_c[chrom][i])

		# 		summary_value = round(sum(window_dq),2)

		# 		if summary_value > 0:
		# 			kernel_c[chrom][i] = summary_value
				



		# 		perc_out = perc.get_percent(ii)
		# 		if perc_out:
		# 			print(f"   calculating kernel density .. {perc_out}% \t{round(ii/1000000, 1)}M nt ", end='\r', flush=True)


		# 		if len(window_dq) > kernel_window:
		# 			window_dq.popleft()


		# print()

		return(depth_c)
	## Getting positional coverage accross the alignment


	start = time()
	if coverage_method == "depth":
		peak_c = get_peaks(
			bam=alignment_file, 
			rgs=annotation_readgroups)

	elif coverage_method == "kernel":
		peak_c = get_kernel_peaks(
			bam=alignment_file, 
			rgs=annotation_readgroups)

	elapsed = time() - start
	print(f"   {round(elapsed,2)}s elapsed")

	print()
	print('prog\tchrom\t\treg\tloc\t\tassess')

	for chrom_count, chrom_and_length in enumerate(chromosomes):


		chrom, chrom_length = chrom_and_length

		# print(f"{chrom_count+1} / {len(chromosomes)}")
		# print(f"chrom: {chrom}")
		# print(f"       {chrom_length:,} bp")
		# print(f"       {chrom_depth_c[chrom]:,} reads")



		status_line = f'{chrom_count+1}/{len(chromosomes)}\t'


		claim_d = {}
		loci = []


		def boolean_creep(direction, cursor, depth_threshold):
			"""
			Defines peak boundaries, based on a minimum threshold in a given direction.

			Peak prominence is tested over a user defined window (pad). When all positions in the window drop below the threshold, the edge is defined.

			Returns this coordinate.
			"""

			positions = deque()
			depths    = deque()

			while True:

				cursor += direction

				if cursor in claim_d:
					break

				positions.append(cursor)
				depths.append(peak_c[chrom][cursor])

				if len(positions) > pad:
					positions.popleft()
					depths.popleft()

					window_average = max(depths)

					if window_average <= depth_threshold:
						break


			boundary = cursor

			return(boundary)

		def creep(direction, cursor, depth_threshold):
			"""
			Defines peak boundaries, based on a minimum threshold in a given direction.

			Peak prominence is tested over a user-defined window (creep_dist). When the average all positions in the window drop below the threshold, an outer edge is defined. Next user-defined pad (pad) is added.

			Returns this coordinate.
			"""

			# print()
			# print(f"  --creep {direction}--")
			# print("  start position:  ", cursor)
			# print("      peak depth:  ", peak_c[cursor])
			# print("       threshold:  ", depth_threshold)

			positions = deque()
			depths    = deque()

			# positions = deque([x+(r*direction) for r in range(0, merge_dist)])
			# depths    = deque([peak_c[p] for p in positions])

			while True:

				cursor += direction

				if cursor in claim_d:
					break

				positions.append(cursor)
				depths.append(peak_c[chrom][cursor])

				if len(positions) > creep_dist:
					positions.popleft()
					depths.popleft()

					window_average = mean(depths)

					if window_average <= depth_threshold:
						break

			if len(positions) <= 1:
				return(cursor)


			# while True:
			# 	# print("positions", positions)
			# 	# print("depths", depths)
			# 	positions.pop()
			# 	# depths.pop()
				
			# 	if depths.pop() > depth_threshold or len(positions) == 1:
			# 		break


				# window_median = median(depths)

				# # print(window_median)

				# if window_median > depth_threshold or len(depths) == 1:
				# 	break

			# print(cursor)


			# boundary = math.ceil(median(positions)) + (direction * pad)
			boundary = positions[-1] + (direction * pad)

			return(boundary)

		def resolve_overlaps(lbound, center, rbound):


			for p in range(center, lbound-1, -1):

				if p in claim_d:
					lbound = p + 1
					break


			for p in range(center, rbound+1, 1):

				if p in claim_d:
					rbound = p - 1
					break

			return(lbound, rbound)

		candidate_peak_count = len([v for v in peak_c[chrom].values() if v / aligned_read_count * 1000000 > rpm_threshold])

		perc = percentageClass(1, candidate_peak_count)

		i = 0
		for center, depth in peak_c[chrom].most_common():


			perc_out = perc.get_percent(i)
			if perc_out:
				print(f"{chrom_count+1}/{len(chromosomes)}\t{chrom}\t{perc_out}%", end='\r', flush=True)
				# print(f"   forming peaks to regions .... {perc_out}%", end='\r', flush=True)
			i += 1


			if center not in claim_d:

				rpm = depth / aligned_read_count * 1000000
				depth_threshold = math.floor(depth * peak_threshold)


				if rpm > rpm_threshold:

					cl_i += 1

					if boundary_method == 'creep':
						lbound  = creep(direction=-1, cursor=center, depth_threshold=depth_threshold)
						rbound  = creep(direction= 1, cursor=center, depth_threshold=depth_threshold)

					elif boundary_method == 'bool':
						lbound  = boolean_creep(direction=-1, cursor=center, depth_threshold=depth_threshold)
						rbound  = boolean_creep(direction= 1, cursor=center, depth_threshold=depth_threshold)


					lbound, rbound = resolve_overlaps(lbound, center, rbound)

					claim = f"Cluster_{cl_i}"

					for r in range(lbound, rbound + 1):
						claim_d[r] = claim

					locus = [claim, chrom, lbound, rbound]
					loci.append(locus)



					## assessing regions which are actually depth 0
					# sc, ft, sizes, strands = test_by_samtools(locus)

					# if sc.depth == 0:




					# 	print()
					# 	print(f"{lbound:,} <- {center:,} -> {rbound:,}")
					# 	print(loci[-1])
					# 	print(rbound - lbound, "nt")
					# 	input()

		# print()






		## Sorting loci by position

		# print(f"   sorting and writing regions ......", flush=True)

		loci.sort(key=lambda x: x[2])






		with open(region_file, 'a') as outf:

			for l in loci:

				name, chrom, start, stop = l

				line = [chrom, 'yasma_peak','region', start, stop, '.', '.', '.',
			f'ID={name}']

				line = "\t".join(map(str,line))
				print(line, file=outf)





		## Performing clumping by reading whole chromosome alignments



		clump_set = set()




		def get_considered_loci(i):

			if i >= len(loci):
				return(False)

			out = [i]
			n = 1
			current_end = loci[i][3]
			while i+n < len(loci):

				# try:
				next_start = loci[i+n][2]
				# except:

				# 	print(loci[i+n])
				# 	sys.exit()

				if next_start < current_end + clump_dist:
					out.append(i+n)

				else:
					break

				n += 1

			return(out)






		# loci_copy = loci[:]


		# start = time()


		locus_i = 0
		considered_loci = get_considered_loci(locus_i)
		unclumped_loci_count = len(loci)


		print(f"{chrom_count+1}/{len(chromosomes)}\t{chrom}\t100%\t{unclumped_loci_count} -> {len(loci)}", end='\r', flush=True)
		# print(f"   clumping similar neighbors... {unclumped_loci_count} -> {len(loci)} loci    ", end='\r', flush=True)

		last_claim = 'start'

		# test_d = {}
		in_locus = False


		sizecall_d = {}
		strand_d   = {}
		seq_d      = {}
		n = False


		def final_check_and_increment(i):
			if loci[i][3] - loci[i][2] < min_locus_length:
				del loci[i]
			else:
				i += 1
			return(i)



		for read in samtools_view(alignment_file, locus=chrom, rgs=annotation_readgroups):

			## Breaks for the final region
			if not considered_loci:
				break


			## Processing sam output
			sam_strand, sam_length, _, sam_lbound, _, _, sam_seq, read_name = read
			sam_rbound = sam_lbound + sam_length


			## Identifies the current claim

			try:
				rclaim = claim_d[sam_rbound]
			except KeyError:
				rclaim = None

			try:
				lclaim = claim_d[sam_lbound]
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



			if sam_lbound > loci[considered_loci[-1]][3]:
				## Checking for any broken loci

				for c in considered_loci[::-1]:
					c_name = loci[c][0]

					if c_name not in strand_d:
						to_loc   = f"after_{c_name}"
						from_loc = f"after_{loci[c-1][0]}"

						try:
							sizecall_d[to_loc] += sizecall_d[from_loc]
							strand_d[to_loc] += strand_d[from_loc]
							seq_d[to_loc] += seq_d[from_loc]
						except KeyError:
							pass

						del loci[c]

						if show_warnings:
							print(f"Warning: bad locus {c_name} removed")


				considered_loci = get_considered_loci(locus_i)




			if not considered_loci:
				break

			# try:
			# 	sam_lbound > loci[considered_loci[-1]][3]
			# except:

			# 	print()
			# 	print("all loci:")
			# 	print(loci)
			# 	print()
			# 	print("considered_loci:")
			# 	print(considered_loci)
			# 	print()
			# 	print("response:")
			# 	print(loci[considered_loci[-1]])
			# 	print("Warning: unknown problem in locus merging...")
			# 	sys.exit()
			# 	break

			if sam_lbound > loci[considered_loci[-1]][3]:
				## if all good, moving forward



				if len(considered_loci) == 1:

					with open(merge_file, 'a') as outf:
						print('', file=outf)
						pprint([loci[c] for c in considered_loci], outf)
						print('      -> locus has no other regions in range', file=outf)

					locus_i = final_check_and_increment(locus_i)
					# print("\nonly_one_considered <- increment")

					considered_loci = get_considered_loci(locus_i)



				else:

					n = 0
					while True:
						with open(merge_file, 'a') as outf:
							print('', file=outf)
							pprint([loci[c] for c in considered_loci], outf)
							# input()

						none_merged = True
						
						current_locus = loci[locus_i]

						# print(considered_loci)
						for n in considered_loci[1:]:
							next_locus = loci[n]

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
										from_loc = loci[r+1][0]


										try:
											sizecall_d[to_loc] += sizecall_d[from_loc]
											strand_d[to_loc] += strand_d[from_loc]
											seq_d[to_loc] += seq_d[from_loc]
											print(f"  packing -> {to_loc} with {sizecall_d[from_loc].depth} reads from {from_loc}", file=outf)
										except KeyError:
											print("Warning: unknown KeyError in locus packing...")
											# sys.exit("error!")
											# print(f"\nWarning: region {next_locus} may not have counts")

											# del loci[locus_i]
											# considered_loci = get_considered_loci(locus_i) 
											pass

										to_loc   = current_locus[0]
										from_loc = f"after_{loci[r][0]}"

										try:
											sizecall_d[to_loc] += sizecall_d[from_loc]
											strand_d[to_loc] += strand_d[from_loc]
											seq_d[to_loc] += seq_d[from_loc]
											print(f"  packing -> {to_loc} with {sizecall_d[from_loc].depth} reads from {from_loc}", file=outf)
										except KeyError:
											# sys.exit("error!")
											# print(f"\nWarning: region {next_locus} may not have counts")
											# del loci[locus_i]
											# considered_loci = get_considered_loci(locus_i) 
											pass

									to_loc   = f"after_{current_locus[0]}"
									from_loc = f"after_{loci[n][0]}"

									try:
										sizecall_d[to_loc] += sizecall_d[from_loc]
										strand_d[to_loc] += strand_d[from_loc]
										seq_d[to_loc] += seq_d[from_loc]
										print(f"  packing -> {to_loc} with {sizecall_d[from_loc].depth} reads from {from_loc}", file=outf)
									except KeyError:
										# print(f"\nWarning: region {next_locus} may not have counts")

										# sys.exit("error!")
										# del loci[locus_i]
										# considered_loci = get_considered_loci(locus_i) 
										pass



								## Redefining current locus boundaries
								loci[locus_i][3] = loci[n][3]


								## Eliminating loci up to the one that was merged
								for r in range(locus_i, n)[::-1]:
									del loci[r+1]


								## Filling in claim_d
								for r in range(loci[locus_i][2],loci[locus_i][3]+1):
									claim_d[r] = loci[locus_i][0]


								break


						## updating considered loci after merging
						considered_loci = get_considered_loci(locus_i) 

						if none_merged:
							with open(merge_file, 'a') as outf:
								print("      -> no valid merges for considered_loci", file=outf)
							locus_i = final_check_and_increment(locus_i)
							# print("\nnone_merged <- increment")
							considered_loci = get_considered_loci(locus_i) 
							break

						if len(considered_loci) == 1:
							with open(merge_file, 'a') as outf:
								print('      -> new locus has no other regions in range', file=outf)
							locus_i = final_check_and_increment(locus_i)
							# print("\nlen(considered_loci) <- increment")
							considered_loci = get_considered_loci(locus_i) 
							break

						if sam_lbound < loci[considered_loci[-1]][3]:
							with open(merge_file, 'a') as outf:
								print("      -> loci passed read location", file=outf)
							last_claim = loci[considered_loci[-1]][0]
							break




				# print(f"   clumping similar neighbors... {unclumped_loci_count} -> {len(loci)} loci    ", end='\r', flush=True)

				print(f"{chrom_count+1}/{len(chromosomes)}\t{chrom}\t100%\t{unclumped_loci_count} -> {len(loci)}      ", end='\r', flush=True)


		# print()
		to_delete = []
		for i,locus in enumerate(loci):
			name = locus[0]

			if name not in strand_d:
				print("Serious problem - loci with no reads broke through:", name)

				to_delete.append(i)

		for i in to_delete[::-1]:
			del loci[i]




		## Assessing locus dimensions and making annotations

		# read_depths = []

		perc = percentageClass(increment=5, total=len(loci))


		last_stop = 0
		for i,locus in enumerate(loci):


			print_percentage = perc.get_percent(i)
			if print_percentage:

				print(f"{chrom_count+1}/{len(chromosomes)}\t{chrom}\t100%\t{unclumped_loci_count} -> {len(loci)}\t{print_percentage}%", end='\r', flush=True)
				# print(f"   assessing loci .............. {print_percentage}%", end="\r", flush=True)

			name, chrom, start, stop = locus
			coords = f"{chrom}:{start}-{stop}"

			# print(locus)
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

			if results_line:
				with open(results_file, 'a') as outf:
					print("\t".join(map(str, results_line)), file=outf)

				with open(gff_file, 'a') as outf:
					print("\t".join(map(str, gff_line)), file=outf)

				top_reads_save(read_c, reads_file, read_equivalent, name)

		if loci:
			print(f"{chrom_count+1}/{len(chromosomes)}\t{chrom}\t100%\t{unclumped_loci_count} -> {len(loci)}\t100%", end='\r', flush=True)

		all_loci += loci
		# sys.exit()

		# print()
		# print()


		locus_lengths = [l[3]-l[2] for l in loci]
		try:
			mean_length   = round(mean(locus_lengths),1)
		except StatisticsError:
			mean_length = None
		try:
			median_length = median(locus_lengths)
		except StatisticsError:
			median_length = None

		if len(loci) > 0:
			proportion_chromosome_annotated = round(sum(locus_lengths) / chrom_length, 4)
		else:
			proportion_chromosome_annotated = None

		# print(f"     length:")
		# print(f"       mean ---> {mean_length} bp")
		# print(f"       median -> {median_length} bp")
		# print(f"       proportion of chromosome annotated -> {proportion_chromosome_annotated}")

		read_depths = [sum(strand_d[l[0]].values()) for l in loci]
		try:
			mean_depth   = round(mean(read_depths),1)
		except StatisticsError:
			mean_depth = None
		try:
			median_depth = median(read_depths)
		except StatisticsError:
			median_depth = None
		if len(loci) > 0:
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
			out += [unclumped_loci_count, len(loci)]
			out += [chrom_length, proportion_chromosome_annotated, mean_length, median_length]
			out += [chrom_depth_c[chrom], proportion_libraries_annotated, mean_depth, median_depth]
			print("\t".join(map(str, out)), file=outf)

		# print()


		try:
			overall_d['region_count']  += unclumped_loci_count
			overall_d['loci_count']    += len(loci)
			overall_d['genome_length'] += chrom_length
			overall_d['locus_lengths'] += locus_lengths
			overall_d['total_depth']   = aligned_read_count
			overall_d['read_depths']   += read_depths

		except KeyError:
			overall_d['region_count']   = unclumped_loci_count
			overall_d['loci_count']     = len(loci)
			overall_d['genome_length']  = chrom_length
			overall_d['locus_lengths']  = locus_lengths
			overall_d['total_depth']    = aligned_read_count
			overall_d['read_depths']    = read_depths


		print()



	print()
	print()


	with open(overall_file, 'a') as outf:


		print('project\tannotation_name\tregion_count\tloci_count\tgenome_length\tproportion_genome_annotated\tmean_length\tmedian_length\ttotal_depth\tproportion_library_annotated\tmean_depth\tmedian_depth', file=outf)

		line = [
			project_name,
			annotation_name,
			overall_d['region_count'], 
			overall_d['loci_count'], 
			overall_d['genome_length']
		]

		if overall_d['loci_count'] == 0:
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

		if overall_d['loci_count'] == 0:
			line += ['NA', "NA", 'NA']
		else:
			line += [
				round(sum(overall_d['read_depths'])/overall_d['total_depth'], 4),
				round(mean(overall_d['read_depths']),1),
				median(overall_d['read_depths'])
			]

		print("\t".join(map(str, line)), file=outf)


		# print('region_count ..................', overall_d['region_count'], file=outf)
		# print('loci_count ....................', overall_d['loci_count'], file=outf)
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















