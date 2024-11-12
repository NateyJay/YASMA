

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

np.seterr(all='raise')


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

		if not sizes:
			return

		if type(sizes) == int:
			sizes = [sizes]

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
		self.header += ['Gap', 'skew', 'size_1n','size_1n_depth', 'size_2n','size_2n_depth', 'size_3n','size_3n_depth', 'sizecall', 'condition']



	def format(self, locus, seq_c, strand_c, sizecall, aligned_depth, last_stop, condition):

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
		result_line += [condition]


		if 'N' in sizecall.sizecall:
			feature_type = "OtherRNA"
		else:
			feature_type = f"RNA_{sizecall}"

		if start < 1:
			start = 1
		gff_line = [
			chrom, 'yasma_locus',feature_type, start, stop, '.', strand, '.',
			f'ID={name};sizecall={sizecall};condition={condition};depth={depth};rpm={rpm};fracTop={frac_top};complexity={complexity};skew={skew};majorRNA={major_rna}'
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







@cli.command(group='Annotation', help_priority=2)


@optgroup.group('\n  Basic options',
				help='')

@optgroup.option("-a", "--alignment_file", 
	required=False, 
	type=click.UNPROCESSED, callback=validate_path,
	help='Alignment file input (bam or cram).')

@optgroup.option('-ac', '--annotation_conditions', 
	required=False,
	multiple=True,
	default=None,
	help="List of conditions names which will be included in the annotation. Defaults to use all libraries, though this is likely not what you want if you have multiple groups.")

@optgroup.option("-o", "--output_directory", 
	# default=f"Annotation_{round(time())}", 
	required=True,
	type=click.Path(),
	help="Directory name for annotation output.")

@optgroup.option("-n", "--name", 
	# default=f"Annotation_{round(time())}", 
	required=False,
	type=str,
	help="Optional name for annotation. Useful if comparing annotations.")


@optgroup.option("-c", "--conditions", 
	required=False, 
	multiple=True,
	type=click.UNPROCESSED, callback=validate_condition,
	help='Values denoting condition groups (sets of replicate libraries) for projects with multiple tissues/treatments/genotypes. Can be entered here as space sparated duplexes, with the library base_name and condition groups delimited by a colon. E.g. SRR1111111:WT SRR1111112:WT SRR1111113:mut SRR1111114:mut')



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
	default=2,
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


@optgroup.option('--tradeoff_round', default=4, help='Significance rounding for tradeoff average. Defaults to 3 digits (e.g. 0.977 or 97.7%)')





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


@optgroup.option('--filter_skew/--no_filter_skew', default=False, help='filter highly skewed loci (default: False)')

@optgroup.option("--max_skew",
	default=0.90,
	type=float,
	help="Filter value for loci which are skewed toward only one sequence in abundance. By default (0.95), if more than 1 in 20 reads for a locus are a single sequence, they are excluded from the annotation.")

@optgroup.option('--filter_complexity/--no_filter_complexity', is_flag=True, default=False, help='filter low complexity loci (default: False)')
@optgroup.option("--min_complexity",
	default=10,
	type=int,
	help="Filter value for locus complexity. This is defined as the number of unique-reads / 1000 nt (default: 10).")

@optgroup.option('--filter_abundance/--no_filter_abundance', is_flag=True, default=True, help='filter low abundance loci (default: True). This is meant to remove loci which have not reached an absolute level of abundance.')
@optgroup.option("--min_abundance",
	default=50,
	type=int,
	help="Min reads in a locus")

@optgroup.option('--filter_abundance_density/--no_filter_abundance_density', is_flag=True, default=True, help='filter low abundance loci (default: True). This is meant to remove loci which have not reached an absolute level of abundance.')
@optgroup.option("--min_abundance_density",
	default=100,
	type=int,
	help="Min reads per 1000 nucleotides in a locus")


@optgroup.group('\n Other options',
				help='')

@optgroup.option('--force', is_flag=True, default=False, help='force resubsample')
@optgroup.option('--debug', is_flag=True, default=False, help='Debug flag')
@optgroup.option('--test_mode', is_flag=True, default=False, help='test_mode flag')
@optgroup.option('--override', is_flag=True, default=False, help='Overrides config file changes without prompting.')


def tradeoff(**params):
	'''Annotator using large coverage window and prec/sens tradeoff.'''

	rc = requirementClass()
	# rc.add_samtools()
	rc.check()

	ic = inputClass(params)
	ic.check(['alignment_file'])

	output_directory        = str(ic.output_directory)
	alignment_file          = ic.inputs["alignment_file"]
	conditions              = ic.inputs['conditions']
	project_name            = ic.inputs['project_name']
	annotation_conditions   = ic.inputs['annotation_conditions']

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



	### Organizing directories

	if annotation_name:
		dir_name = f"tradeoff_{annotation_name}"
	# elif annotation_conditions:
	# 	dir_name = f"tradeoff_{'_'.join(annotation_conditions)}"
	else:
		dir_name = 'tradeoff'

	Path(output_directory, dir_name).mkdir(parents=True, exist_ok=True)



	### Processing basic options
	if params['target_genome_perc'] and params['target_read_perc']:
		sys.exit("ERROR: cannot specify target read AND genome percentages (one is dependent on the other)")



	if debug: 
		show_warnings = True
	else:
		show_warnings = False



	### Starting log
	log_file = f"{output_directory}/{dir_name}/log.txt"
	sys.stdout = Logger(log_file)


	start_time = datetime.now()
	date_time = start_time.strftime("%Y/%m/%d, %H:%M:%S")
	print()
	print("Start time:",date_time)	

	stat_d = {} ## an object which records basic stats of the annotation
	# project	
	# annotation_name	
	# region_count	
	# locus_count	
	# genome_length	
	# proportion_genome_annotated	
	# mean_length	median_length	
	# total_depth	
	# proportion_library_annotated	
	# mean_depth	
	# median_depth


	### Writing annotation parameters to folder

	params['tool'] = 'tradeoff'
	parameters_file = Path(output_directory, dir_name, "params.json")
	with open(parameters_file, 'w') as outf:
		oparams = params
		oparams['alignment_file'] = str(params['alignment_file'])
		outf.write(json.dumps(params, indent=2))
		del oparams
		# sys.exit()


	### Getting alignment metadata

	chromosomes, libraries = get_chromosomes(alignment_file)



	### getting basic metrics, including a test_mode chromosome filter

	if params['test_mode']:
		# chromosomes = chromosomes[3:5]
		# chromosomes = chromosomes[20:30]
		# chromosomes = chromosomes[2:5]
		# chromosomes = chromosomes[:2]
		# chromosomes = chromosomes[:1]
		chromosomes = chromosomes[4:5]


	chromosome_max_lengths = {}
	for c,l in chromosomes:
		chromosome_max_lengths[c] = l


	genome_length = sum([l for c,l in chromosomes])


	# print(alignment_file)
	# annotation_readgroups = check_rgs(new_rgs, bam_rgs)
	# params['annotation_readgroups'] = annotation_readgroups



	### cleaning up conditions and libraries

	if not conditions:
		conditions = {'all' : libraries}


	if annotation_conditions:
		print(f" Annotation conditions: {annotation_conditions}")

	else:
		print(f" Annotation conditions: {annotation_conditions} (all libraries used)")
		annotation_conditions = list(conditions.keys())
	# print(" Annotation libraries:", annotation_readgroups)


	conditions = {c: conditions[c] for c in annotation_conditions}


	rev_conditions = {}
	for cond, lib_set in conditions.items():
		for lib in lib_set:
			rev_conditions[lib] = cond

	libraries  = []
	for lib_set in conditions.values():
		for lib in lib_set:
			libraries.append(lib)




	print()


	### getting alignment depths by library and chromosome

	chrom_depth_c = get_global_depth(alignment_file, aggregate_by=['rg','chrom'])

	keys = list(chrom_depth_c.keys())
	for key in keys:
		if key[0] in libraries:
			chrom_depth_c[key[1]] += chrom_depth_c[key]

		del chrom_depth_c[key]

	for key in list(chrom_depth_c.keys()):
		if key not in [c for c,l in chromosomes]:
			del chrom_depth_c[key]

	aligned_read_count = sum(chrom_depth_c.values())


	### getting Read Per Billion equivalents for all the libraries

	g = get_global_depth(alignment_file, aggregate_by=['rg'])

	rpb_d = {}
	for lib in libraries:
		condition = rev_conditions[lib]
		depth = g[lib]
		rpb = round(1000000000 / depth, 2)

		rpb_d[lib] = rpb

		print(f"    {condition} {lib} -> {depth:,} -> {rpb:,} scaled read")


	rpbs = np.zeros(shape=(len(annotation_conditions), max([len(conditions[c]) for c in annotation_conditions])), dtype='uint32')
	for i,c in enumerate(annotation_conditions):
		for j,s in enumerate(conditions[c]):
			try:
				rpbs[i,j] = round(rpb_d[s])
			except KeyError:
				rpbs[i,j] = 0


	# print(rpbs)
	# print(rpbs.shape)

	### calling subsampling if needed

	if params['subsample']:

		alignment_file = subsample(aligned_read_count, alignment_file, params)
		chrom_depth_c = get_global_depth(alignment_file, aggregate_by=['rg','chrom'])
		aligned_read_count = sum(chrom_depth_c.values())



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


	filter_file = Path(output_directory, dir_name, 'filtered_loci.txt')
	with open(filter_file, 'w') as outf:
		print('coords\tlength\tabd\tpass_abd\tabd_dens\tpass_abd_dens\tcomplexity\tpass_complexity\tskew\tpass_skew', sep='\t', file=outf)


	# stats_file = f"{output_directory}/{dir_name}/stats_by_ref.txt"
	# with open(stats_file, 'w') as outf:
	# 	print("project\tchromosome\tregions\tregions\tchromosome_length\tproportion_genome_annotated\tmean_length\tmedian_length\treads\tproportion_libraries_annotated\tmean_abundance\tmedian_abundance", file=outf)


	# overall_file = f"{output_directory}/{dir_name}/stats.txt"
	# with open(overall_file, 'w') as outf:
	# 	outf.write('')

	# overall_d = {}





	## iterating through chromosomes and reads

	print()
	print(f" {len(chromosomes)} chromosomes/scaffolds/contigs")
	print(f" {genome_length:,} bp total")
	print(f" {aligned_read_count:,} reads")
	print()


	stat_d['genome_length'] = genome_length
	stat_d['aligned_reads'] = aligned_read_count


	## Getting positional coverage accross the alignment


	start = time()



	def get_kernel_coverage():
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

		def get_max_lib_counts():
			lib_counts = []
			for srrs in conditions.values():
				lib_counts.append(len(srrs))
			return max(lib_counts)

		max_lib_counts = get_max_lib_counts()

		pos_size_d = dict()
		pos_d      = dict()

		print(" processing alignment...")

		for chrom, chrom_length in chromosomes:

			pos_size_d[chrom] = [None] * chrom_length
			pos_d[chrom] = np.zeros(shape=(chrom_length, len(annotation_conditions), max_lib_counts, 2), dtype='uint32')



		ec = elapsedClass()
		iterables = []
		for c,l in chromosomes:
			iterables.append(samtools_view(alignment_file, rgs=libraries, contig=c))#, read_minmax=(params['min_read_length'], params['max_read_length'])))

		reads = chain.from_iterable(iterables)

		# print(f"    encoding reads ............... 0%", end='\r', flush=True)
		perc = percentageClass(1, sum(chrom_depth_c.values()))
		perc.update()

		aligned_read_count = 0
		for i, read in enumerate(reads):
			aligned_read_count+=1
			strand, length, _, pos, chrom, lib, _, _ = read
			condition = rev_conditions[lib]


			pos += math.floor(length / 2)

			strand_i    = ["+", "-"].index(strand)
			condition_i = annotation_conditions.index(condition)
			lib_i       = conditions[condition].index(lib)


			pos_d[chrom][pos, condition_i, lib_i, strand_i] += 1

			try:
				pos_size_d[chrom][pos].append(length)
			except AttributeError:
				pos_size_d[chrom][pos] = [length]



			perc_out = perc.update()
			if perc_out:
				sys.stdout.write(f"    encoding reads ............... {perc_out}%\t {i+1:,} reads   \n", terminal_only=True)
				sys.stdout.flush()
				sys.stdout.overwrite_lines(1)


		# print()
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
			# ec = elapsedClass()

			gen_c = Counter() # a counter of kernel depths at positions in the genome 
			read_c = Counter() # a counter of the max kernel depth across the span of a read.



			cov_track = trackClass(Path(output_directory, dir_name, "coverage.bw"), chromosomes)
			ker_track = trackClass(Path(output_directory, dir_name, "kernel.bw"), chromosomes)


			# rpbs = np.array([round(d['rpb']) for d in lib_d.values()], dtype='float32')

			for chrom_i, chromosome_entry in enumerate(chromosomes):
				chrom, chrom_length = chromosome_entry

				sys.stdout.write(f"    computing coverage ........... {chrom_i+1}/{len(chromosomes)} {chrom}                \n", terminal_only=True)
				sys.stdout.flush()
				sys.stdout.overwrite_lines(1)

				# print("    rpm array...")

				# print("    calc coverage...")
				# print(pos_d[chrom].dtype)
				coverage = np.sum(pos_d[chrom], axis=3, dtype='uint32')

				if coverage.shape[0] < cov_window or coverage.shape[0] < ker_window:
					print(chrom, 'pass, too short!')

					coverage = np.zeros(shape=(coverage.shape[0]), dtype='uint32')
					kernel   = np.zeros(shape=(coverage.shape[0]), dtype='uint32')

				else:
					# print(coverage.dtype)
					# print(coverage.shape, "<- raw reads")
					coverage = np.multiply(coverage, rpbs)
					# print(coverage.dtype)
					# print(coverage.shape, "<- rpb normalized")
					coverage = np.median(coverage, axis=(2))
					# print(coverage.dtype)
					# print(coverage.shape, "<- median of libraries")
					coverage = coverage.astype('float32')
					# print(coverage.dtype)
					# print(coverage.shape, "<- condensint to f32")
					coverage = np.mean(coverage, axis=(1))
					# print(coverage.dtype)
					# print(coverage.shape, "<- mean of conditions")
					coverage = sliding_window_view(coverage, cov_window, axis=0)
					# print(coverage.dtype)
					# print(coverage.shape, "<- windows")
					coverage = np.mean(coverage, axis=1)
					# print(coverage.dtype)
					# print(coverage.shape, "<- mean kernel")
					coverage = np.concatenate(
						(np.array([coverage[0]] * half_cov_window), 
							coverage, 
							np.array([coverage[-1]] * (half_cov_window-1))), 
						axis=0)
					# print(coverage.shape, "<- padded")



					# print("    calc kernel...")
					kernel = sliding_window_view(coverage, ker_window, axis=0)
					# print(kernel.shape, "<- k windows")
					kernel = np.max(kernel, axis=1)
					# print(kernel.dtype)
					# print(kernel.shape, "<- max kernel")
					kernel = np.concatenate(
						(np.array([kernel[0]] * half_cov_window), 
							kernel, 
							np.array([kernel[0]] * (half_cov_window-1))), 
						axis=0)
					# print(kernel.shape, "<- padded")


					# print("    write bigwigs...")

				# print(coverage)
				for i, c in enumerate(coverage):
					cov_track.add(chrom, i+1, float(c))

				for i, k in enumerate(kernel):
					ker_track.add(chrom, i+1, float(k))

				# print("    tally...")

				# gen_c.update([k for k in kernel])

				summed = np.sum(pos_d[chrom], axis=(1,2,3), dtype='uint32')

				for i, k in enumerate(kernel):
					v = int(summed[i])

					k = round(k, 2)

					read_c[k] += v
					gen_c[k]  += 1


			print(f"    computing coverage ........... {chrom_i+1}/{len(chromosomes)} {chrom}                ", end='\n', flush=True)

			print()
			cov_track.close()
			ker_track.close()



			print()
			print(ec)
			print()

			return(gen_c, read_c)

		gen_c, read_c = numpy_method()


		### gen_c is a counter object
		# key: floor(rpm) value based on kernel
		# value: the number of ==genomic spaces== with that rpm value

		### read_c is a counter object
		# key: floor(rpm) value based on kernel
		# value: the number of ==reads== in total with that rpm value



		# print("gen_c")
		# pprint(gen_c.most_common(20))
		# print(sum(gen_c.values()), "<- calculated genome length")

		# print()
		# print("read_c")

		# pprint(read_c.most_common(20))
		# print(sum(read_c.values()), "<- calculated read count")


		def knee(gen_c, read_c):
			found_depths = list(gen_c.keys())
			found_depths.sort()

			total_genomic_space  = sum([l for c,l in chromosomes])
			total_read_space     = sum(read_c.values())


			thresholds      = []
			annotated_space = []
			p_gens          = []
			annotated_reads = []
			p_reads         = []
			averages        = []
			kdiffs          = []

			genp_thresholds  = []
			readp_thresholds = []


			for threshold in found_depths:

				total_genomic_space -= gen_c[threshold]
				total_read_space    -= read_c[threshold]

				p_read = total_read_space / aligned_read_count
				p_gen  = total_genomic_space / genome_length


				genp_thresholds.append((p_gen, threshold))
				readp_thresholds.append((p_read, threshold))

				average = mean([p_read, p_gen])

				vdist = p_read - p_gen
				pdist = math.sqrt( (average - p_gen)**2 * 2 )
				kdiff = vdist - pdist


				thresholds.append(threshold)
				annotated_space.append(total_genomic_space)
				annotated_reads.append(total_read_space)
				p_gens.append(p_gen)
				p_reads.append(p_read)
				averages.append(average)
				kdiffs.append(kdiff)


			peak_i = max(range(len(kdiffs)), key=kdiffs.__getitem__)


			with open(Path(output_directory, dir_name, 'thresholds.txt'), 'w') as outf:
				print('depth\tannotated_space\tp_genome\tannotated_reads\tp_reads\taverage_score\tkdiff\tpeak', file=outf)

				for i, threshold in enumerate(found_depths):

					space     = annotated_space[i]
					p_gen     = p_gens[i]
					reads     = annotated_reads[i]
					p_read    = p_reads[i]
					average   = averages[i]
					kdiff     = kdiffs[i]

					threshold = round(threshold, 3)


					if i == peak_i:
						peak = 1
						out = {
							'threshold' : threshold,
							'p_gen'     : p_gen,
							'genome'    : space,
							'p_read'    : p_read,
							'reads'     : reads,
							'average'   : average,
							'kdiff'     : kdiff
							}
					else:
						peak = 0

					print(threshold, space, p_gen, reads, p_read, average, kdiff, peak, sep='\t', file=outf)

			return(out, readp_thresholds, genp_thresholds)
		

		out, readp_thresholds, genp_thresholds = knee(gen_c, read_c)
		print(out)



		# def tally(gen_c, read_c):
		# 	found_depths = list(gen_c.keys())
		# 	found_depths.sort()

		# 	averages = []
		# 	table    = []

		# 	total_genomic_space  = sum([l for c,l in chromosomes])
		# 	total_possible_space = genome_length - gen_c[0]
		# 	total_read_space     = sum(read_c.values())
		# 	genp_thresholds      = []
		# 	adj_genp_thresholds  = []
		# 	readp_thresholds     = []

		# 	for rpm_threshold in found_depths:

		# 		total_genomic_space -= gen_c[rpm_threshold]
		# 		total_read_space    -= read_c[rpm_threshold]



		# 		# print(depth_threshold, total_genomic_space, sep='\t')

		# 		gen_score     = total_genomic_space / genome_length
		# 		adj_gen_score = total_genomic_space / total_possible_space
		# 		read_score    = total_read_space / aligned_read_count

		# 		genp_thresholds.append((gen_score, rpm_threshold))
		# 		adj_genp_thresholds.append((adj_gen_score, rpm_threshold))
		# 		readp_thresholds.append((read_score, rpm_threshold))

		# 		## unweighted avg
		# 		avg_score = ((1-gen_score) + read_score) / 2


		# 		## weight avg
		# 		weight_score = ((1-gen_score) * params['genome_weight'] + read_score * params['read_weight']) / sum([params['genome_weight'], params['read_weight']])



		# 		geom_score = math.sqrt(((1-gen_score) * read_score))


		# 		averages.append(round(weight_score, params['tradeoff_round']))

		# 		table.append([rpm_threshold, 
		# 			total_genomic_space, 
		# 			round(gen_score,4), 
		# 			round(adj_gen_score, 4),
		# 			total_read_space, 
		# 			round(total_read_space/aligned_read_count,4),
		# 			round(avg_score, 4), 
		# 			round(geom_score, 4), 
		# 			round(weight_score, 4)])


		# 	peak_index = averages.index(max(averages))

		# 	with open(Path(output_directory, dir_name, 'thresholds.txt'), 'w') as outf:
		# 		print('depth\tannotated_space\tp_genome\tadj_p_genome\tannotated_reads\tp_reads\taverage_score\tgeom_score\tweighted_avg\tpeak', file=outf)
		# 		for i,t in enumerate(table):

		# 			# print(i,t)
		# 			# input()

		# 			if i == peak_index:
		# 				peak = 1
		# 				out = {
		# 					'rpm_threshold' : t[0],
		# 					'gen_score' : t[2],
		# 					'adj_gen_score' : t[3],
		# 					'read_score' : t[5],
		# 					'avg_score' : t[6],
		# 					'weighted_avg' : t[8]
		# 					}
		# 			else:
		# 				peak = 0

		# 			print("\t".join(map(str, t)), peak, sep='\t', file=outf)
		# 			if total_genomic_space > genome_length:
		# 				sys.exit("problem!!")

		# 	# pprint(out)
		# 	return(out, readp_thresholds, genp_thresholds)
		

		# out, readp_thresholds, genp_thresholds = tally(gen_c, read_c)
		# print(out)
		# print()


		return(pos_d, pos_size_d, out, readp_thresholds, genp_thresholds)





	pos_d, pos_size_d, threshold_stats, readp_thresholds, genp_thresholds = get_kernel_coverage()

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
		print(f"    depth threshold: ......... {depth_threshold} rpb")
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
		print(f"    depth threshold: ......... {depth_threshold} rpb")
		print(f"    exp. genome proportion: .. {gen_score}")
		print(f" -> set read proportion: ..... {read_score}")

	else:
		print(f"Finding threshold through weighted tradeoff. Weight: [{params['read_weight']}] reads to [{params['genome_weight']}] genome")

		depth_threshold = threshold_stats['threshold']
		gen_score       = threshold_stats['p_gen']
		# adj_gen_score   = threshold_stats['adj_gen_score']
		read_score      = threshold_stats['p_read']

		print(" annotation parameters...")
		print(f"    depth threshold: ......... {round(depth_threshold,2):,} rpb")
		print(f"    exp. genome proportion: .. {gen_score}")
		print(f"    exp. read proportion: .... {read_score}")


	if threshold_stats['threshold'] == 0.0:
		sys.exit("Error: detected threshold for annotation is 0 reads per million (0 reads).\nThis will annotate 100%% of reads, leading to a highly unrepresentative sample. This might be caused by problems with the alignment, libraries (check these), or an internal problem with YASMA (please make an issue on github)")



	def get_regions(depth_threshold, chromosomes):

		regions = []
		def check_and_cash_region(in_region, reg_i, chrom, reg_start, reg_stop, chrom_length):
			if in_region:

				if reg_stop >= chrom_length:
					reg_stop = chrom_length - 1


				if np.sum(pos_d[chrom][reg_start:reg_stop,...]) == 0:
					return(reg_i)
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
				if depth > depth_threshold:
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

	stat_d['regions'] = len(all_regions)


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
		total_annotated_reads += np.sum(pos_d[chrom][start:(stop+1),])

	# print(total_annotated_reads)
	print(f"    {total_annotated_reads:,} reads ({ round(total_annotated_reads / aligned_read_count *100,1) }%) in regions")
	print(f"        expected: {round(100*read_score,1)}%")


	stat_d['total_region_space'] = total_region_space
	stat_d['total_region_reads'] = total_annotated_reads

	# sys.exit("holding here...")


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


			p = np.sum(pos_d[self.chrom][self.start:self.stop+1, ...], axis=(0,1,2))

			self.strands['+'], self.strands['-'] = p

			try:
				self.frac_top = self.strands['+'] / sum(self.strands.values())
			except ZeroDivisionError:
				print(self.region)

				sys.exit("ZeroDivisionError - regions")

			except FloatingPointError:
				print(self.region)
				print(p)
				sys.exit("FloatingPointError - regions")


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
				# print(window_start, self.start, window_end, self.stop)
				return(False, ['outofbounds'])


			window = list(range(window_start, window_end)) 
			window_size = len(window)

			w_strands = Counter()
			w_sizes   = sizeClass(minmax=read_minmax)
			w_depths  = 0

			p = np.multiply(np.sum(pos_d[chrom][window_start:window_end+1, ...], axis=(0,3)), rpbs)
			p = np.median(p, axis=1)
			p = np.mean(p, axis=0)

			w_depths = p


			p = np.sum(pos_d[self.chrom][window_start:window_end+1,], axis=(0,1,2))

			w_strands['+'] = p[0]
			w_strands['-'] = p[1]

			for w in range(window_start, window_end+1):
				try:
					sizes = pos_size_d[chrom][w]
				except KeyError:
					pass
				if sizes:
					w_sizes.update(sizes)


			if expand: 
				## trimming doesnt worry about hitting other loci
				if w in claim_d[chrom] and claim_d[chrom][w] != name:
					# print('claim break')
					return(False, ['claim'])

			try:
				w_fractop = w_strands['+'] / sum(w_strands.values())
			except ZeroDivisionError:
				w_fractop = 1
				return(False, ['no reads detected (error?)'])
			except FloatingPointError:
				return(False, ['no reads detected (error?)'])


			expand_fails = []

			if w_depths == 0:
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
					print("Warning: OOB")
					# sys.exit("OOB error")

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

			try:
				rc = reviseClass(region)
			except AttributeError:
				print(f"Warning: {region} encountered an error (empty)")
				# print(f"  to check: samtools view {alignment_file} {chrom}:{start}-{stop}")
				continue

			new_start, new_stop = rc.expand()

			if new_stop > chromosome_max_lengths[chrom]:
				new_stop = chromosome_max_lengths[chrom] - 1

			if new_start > new_stop:
				print(region)
				print(new_start, new_stop)
				print("expand failed!! stop before start")
				sys.exit("EF ERROR")

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


		total_revised_reads += np.sum(pos_d[chrom][start:stop+1,])

		# for r in range(start, stop+1):
		# 	try:
		# 		total_revised_reads += sum(pos_d[chrom][start:stop+1,])
		# 	except IndexError:
		# 		print(chrom, r, "<- index error 4")
		# 		pass


	def string_plus_white(s, length = 7):
		s = str(s)
		return s + " " * (length - len(s))



	print(f"    {revised_genomic_space:,} genomic nt ({round(revised_genomic_space/genome_length *100,1)}%) in revised regions")
	print(f"    {total_revised_reads:,} reads ({round(total_revised_reads/aligned_read_count *100,1)}%) in revised regions")
	diff = round(total_revised_reads/aligned_read_count *100,1) - round(total_annotated_reads/aligned_read_count *100,1)
	print(f"      +{round(diff,2)}% over unrevised regions")



	stat_d['total_revised_space'] = revised_genomic_space
	stat_d['total_revised_reads'] = total_revised_reads




	max_chrom_word_length = max([len(c) for c,l in chromosomes])

	class progressClass():
		def __init__(self, chrom, region_count=None):
			max_chrom_word_length = max([len(c) for c,l in chromosomes])

			self.header = "\t".join(['prog', "chrom"+(max_chrom_word_length-5)*" ",'regions','loci','assess'])

			self.chrom          = chrom
			self.region_count   = region_count
			self.locus_count    = region_count
			self.assessed_count = 0

			self.i = [c for c,l in chromosomes].index(chrom) + 1
			self.n = len(chromosomes)

		def add_white(self, s, length = 7):
			s = str(s)
			return s + " " * (length - len(s))

		def show(self, region_count=None, locus_count=None, assessed_count=None, terminal_only=True):

			if region_count:
				self.region_count = region_count

			if locus_count:
				self.locus_count = locus_count

			if assessed_count:
				self.assessed_count = assessed_count

			if self.assessed_count == 0 or self.locus_count == 0:
				assess_p = ''
			else:
				# assess_p = round(self.assessed_count / self.locus_count / 20) * 5
				# assess_p = f"{assess_p}%"
				# assess_p = str(self.assessed_count)
				assess_p = self.add_white(self.assessed_count)
	
			to_print = [
				f"{self.i}/{self.n}",
				self.chrom,
				self.add_white(self.region_count),
				self.add_white(self.locus_count),
				assess_p
			]
			to_print = "\t".join(to_print) + "\n"

			if terminal_only:
				sys.stdout.write("\x1b[1A\x1b[2K", terminal_only = True)
			sys.stdout.write(to_print, terminal_only = terminal_only)





	# def print_progress_string(i, n, chrom, input_loci, output_loci, assess=False, terminal_only=False):

	# 	chrom = chrom + (max_chrom_word_length - len(chrom)) * " "

	# 	if not assess:
	# 		assess = ''
	# 	else:
	# 		assess = f"\t{assess}%"

	# 	sys.stdout.write(f"{i+1}/{n}\t{chrom}\t{string_plus_white(input_loci)}\t{string_plus_white(output_loci)}{assess}  \r", 
	# 		terminal_only=terminal_only)
	# 	sys.stdout.flush()

	print()
	# print('prog', "chrom"+(max_chrom_word_length-5)*" ",'regions\tloci\tassess', sep='\t')
	print(progressClass(chrom).header)
	print()

	total_region_space = 0
	regions_name_i = 0
	total_annotated_reads = 0

	# Some filter counters
	complexity_filter = 0
	skew_filter       = 0
	abd_filter        = 0
	abd_dens_filter   = 0

	annotated_space = 0
	annotated_reads = 0


	locus_name_i = 0
	all_loci = []

	for chrom_count, chrom_and_length in enumerate(chromosomes):
		chrom, chrom_length = chrom_and_length

		status_line = f'{chrom_count+1}/{len(chromosomes)}\t'



		regions = [r for r in all_regions if r[1] == chrom]


		sizecall_d = {}
		strand_d   = {}
		seq_d      = {}



		def get_region_stats(chrom, start, stop):

			size= sizeClass(minmax=read_minmax)
			for w in range(start, start):
				try:
					size.update(pos_size_d[chrom][w])
				except IndexError:
					pass

			strand = Counter()
			p = np.sum(pos_d[chrom][start: stop, ...], axis=(0,1,2))
			strand['+'] = int(p[0])
			strand['-'] = int(p[1])

			return(size, strand)

		for i in range(len(regions)):

			curr_region = regions[i]
			curr_name = curr_region[0].replace("region_", "r")
			curr_start, curr_stop = curr_region[2:]

			sizecall_d[curr_name], strand_d[curr_name] = get_region_stats(chrom, curr_start, curr_stop)

			if i == len(regions)-1:
				break

			next_region = regions[i+1]
			betw_name = curr_region[0].replace("region_", "a")
			next_start, next_stop = next_region[2:]
			betw_start, betw_stop = curr_stop + 1, next_start - 1

			sizecall_d[betw_name], strand_d[betw_name] = get_region_stats(chrom, next_start, next_stop)




		def get_ft(strands):

			try:
				ft = strands["+"] / sum(strands.values())
			except ZeroDivisionError:
				ft = 1

			return(ft)

		def check_merge(curr_name, next_name, curr_stop, next_start):

			dist_test = next_start - curr_stop <= params['merge_dist']

			size_test = sizecall_d[curr_name] == sizecall_d[next_name]


			curr_ft = strand_d[curr_name]
			next_ft = strand_d[next_name]

			frac_test = abs(get_ft(curr_ft) - get_ft(next_ft)) < clump_strand_similarity

			return size_test and frac_test and dist_test

		def do_merge(curr_name, betw_name, next_name):
			strands = strand_d[curr_name] + strand_d[next_name] + strand_d[betw_name]
			sizes   = sizecall_d[curr_name] + sizecall_d[next_name] + sizecall_d[betw_name]

			return(strands, sizes)

		i = 0
		pc = progressClass(chrom, len(regions))

		locus_count = len(regions)
		pc.show(locus_count=locus_count)

		if len(regions) == 0:
			# print()
			pc.show(terminal_only=False)
			continue

		if len(regions) == 1:
			loci = [[f"locus_{locus_name_i}"] + regions[0][1:]]


		else:
			curr_name = regions[i][0].replace("region_", "r")
			betw_name = regions[i][0].replace("region_", "a")
			next_name = regions[i+1][0].replace("region_", "r")
			next_start = regions[i][2]

			loc_start = regions[i][2]
			loc_stop  = regions[i][3]

			loci = []
			
			# print_progress_string(chrom_count, len(chromosomes), chrom, len(regions), locus_count, terminal_only=True)
			# progressClass


			while True:

				if check_merge(curr_name, next_name, loc_stop, next_start):

					strand_d[curr_name], sizecall_d[curr_name] = do_merge(curr_name, betw_name, next_name)

					del strand_d[betw_name]
					del strand_d[next_name]
					del sizecall_d[betw_name]
					del sizecall_d[next_name]

					with open(merge_file, 'a') as outf:
						print(i, f"{curr_name} <<< {betw_name}, {next_name}", file=outf)
					locus_count -= 1

					loc_stop  = regions[i+1][3]

				else:

					locus_name_i += 1
					locus_name   = f"locus_{locus_name_i}"
					loci.append([locus_name, chrom, loc_start, loc_stop])

					strand_d[locus_name]   = strand_d.pop(curr_name)
					sizecall_d[locus_name] = sizecall_d.pop(curr_name)


					with open(merge_file, 'a') as outf:
						print(i, f"{locus_name} <<< {curr_name}", file=outf)

					curr_name = next_name
					loc_start = regions[i+1][2]
					loc_stop  = regions[i+1][3]

				i += 1


				pc.show(locus_count=locus_count-1)

				if i == len(regions) - 1:
					break

				betw_name = regions[i][0].replace("region_", "a")
				next_name = regions[i+1][0].replace("region_", "r")
				next_start = regions[i+1][2]



				# print_progress_string(chrom_count, len(chromosomes), chrom, len(regions), locus_count, terminal_only=True)
				# input(i)

				# sleep(0.05)

			all_loci += loci
			# continue

			# print()
			# sys.exit('\n')


		## Assessing locus dimensions and making annotations

		# perc = percentageClass(increment=5, total=len(regions))

		stat_d['unfiltered_loci'] = len(loci)

		last_stop = 0
		for i,locus in enumerate(loci):

			# print()
			# print()
			# print("############")
			# print(locus)

			rc = reviseClass(locus)
			locus[2], locus[3] = rc.trim()
			# print(locus, i)

			old_name = locus[0]
			regions_name_i += 1
			locus[0] = f"locus_{regions_name_i}"
			# print(locus[0])
			# name = locus[0]


			pc.show(assessed_count=i+1)

			# print_percentage = perc.get_percent(i)
			# if print_percentage:
			# 	print_progress_string(chrom_count, len(chromosomes), chrom, len(regions), len(loci), print_percentage, terminal_only=True)

			name, chrom, start, stop = locus
			coords = f"{chrom}:{start}-{stop}"



			def check_best_condition():
				locus_abd = pos_d[chrom][start:stop, ]
				p = np.sum(locus_abd, axis=(0,3))
				p = np.multiply(p, rpbs)
				p = np.median(p, axis=1)
				# print(p)
				best_ci = int(np.argmax(p))
				return(best_ci, locus_abd)


				# pass
			best_condition_i, locus_abd = check_best_condition()
			best_condition = list(conditions.keys())[best_condition_i]


			# print(conditions[best_condition])

			# rgs = [r+".t" for r in conditions[best_condition]]

			# for read in samtools_view(alignment_file, contig=chrom, start=start, stop=stop,#, rgs=rgs, 
				# boundary_rule = 'tight'):

			best_strand = Counter()
			best_size   = sizeClass()
			best_seq    = Counter()

			strand      = Counter()
			size        = sizeClass()
			seq         = Counter()



			for read in samtools_view(alignment_file, contig=chrom, start=start, stop=stop):

				sam_strand, sam_read_length, _, _, _, sam_lib, sam_seq, _ = read


				strand[sam_strand] += 1
				size.update(sam_read_length)
				seq[sam_seq] += 1

				if sam_lib in libraries:

					best_strand[sam_strand] += 1
					best_size.update(sam_read_length)
					best_seq[sam_seq] += 1

				# print(read)

			# print()


			abd = sum(seq.values())
			length = stop-start

			complexity = len(seq.keys()) / length * 1000
			skew       = seq.most_common(1)[0][1] / sum(seq.values())


			pass_complexity = complexity >= params['min_complexity']
			pass_skew       = skew <= params['max_skew']
			pass_abd        = sum(best_seq.values()) >= params['min_abundance']
			pass_abd_dens   = sum(best_seq.values()) / (stop-start) * 1000 >= params['min_abundance_density']


			pass_all_filters = 0

			if not pass_complexity:
				complexity_filter += 1
				pass_all_filters += int(params['filter_complexity'])

			if not pass_skew:
				skew_filter += 1
				pass_all_filters += int(params['filter_skew'])

			if not pass_abd:
				abd_filter += 1
				pass_all_filters += int(params['filter_abundance'])

			if not pass_abd_dens:
				abd_dens_filter += 1
				pass_all_filters += int(params['filter_abundance_density'])


			if pass_all_filters > 0:
				regions_name_i -= 1
				with open(filter_file, 'a') as outf:
					print(coords, length, abd, pass_abd, round(abd/length*1000,4), pass_abd_dens, round(complexity,4), pass_complexity, skew, pass_skew, sep='\t', file=outf)
				continue

			
			annotated_space += length
			annotated_reads += np.sum(locus_abd)

			# print(f"{length} ({annotated_space})\t{abd} ({annotated_reads})")
			# print()


			results_line, gff_line = assessClass().format(locus, best_seq, best_strand, best_size, sum(chrom_depth_c.values()), last_stop, best_condition)


			last_stop = stop

			with open(results_file, 'a') as outf:
				print("\t".join(map(str, results_line)), file=outf)

			with open(gff_file, 'a') as outf:
				print("\t".join(map(str, gff_line)), file=outf)

			top_reads_save(best_seq, reads_file, abd, name)


		pc.show(terminal_only=False)



	def bool_to_check(val):
		if val:
			return "[x]"
		else:
			return "[ ]"

	print()
	print()
	print(f"locus filters:  (x = filter activated)")
	print(f"  {bool_to_check(params['filter_skew'])} {skew_filter} loci are extremely skewed (> {params['max_skew']} prop. abundance is a single sequence)")
	print(f"  {bool_to_check(params['filter_complexity'])} {complexity_filter} loci have low complexity (< {params['min_complexity']} unique reads per 1000 nt)")
	print(f"  {bool_to_check(params['filter_abundance'])} {abd_filter} loci are below min abundance (< {params['min_abundance']} aligned reads)")
	print(f"  {bool_to_check(params['filter_abundance_density'])} {abd_dens_filter} loci are below min abundance density (< {params['min_abundance_density']} aligned reads / 1000 nt)")
	print()
	print(f"  {locus_name_i:,} loci passing activated filter(s)")
	print()
	print(f"final annotation metrics:")
	print(f"  {annotated_space:,} genomic nt ({round(annotated_space / genome_length * 100, 1)}%) in loci")
	print(f"  {annotated_reads:,} reads ({round(annotated_reads / aligned_read_count * 100, 1)}%) in loci")
	print()


	stat_d['total_final_space'] = annotated_space
	stat_d['total_final_reads'] = annotated_reads


	stat_d['filtered_loci'] = locus_name_i

	# with open(overall_file, 'a') as outf:


	# 	print('project\tannotation_name\tregion_count\tlocus_count\tgenome_length\tproportion_genome_annotated\tmean_length\tmedian_length\ttotal_depth\tproportion_library_annotated\tmean_depth\tmedian_depth', file=outf)

	# 	line = [
	# 		project_name,
	# 		annotation_name,
	# 		overall_d['region_count'], 
	# 		overall_d['regions_count'], 
	# 		overall_d['genome_length']
	# 	]

	# 	if overall_d['regions_count'] == 0:
	# 		line += ['NA', "NA", 'NA']
	# 	else:
	# 		line += [
	# 			round(sum(overall_d['locus_lengths'])/overall_d['genome_length'], 4),
	# 			round(mean(overall_d['locus_lengths']),1),
	# 			median(overall_d['locus_lengths'])
	# 		]

	# 	line += [
	# 		overall_d['total_depth']
	# 	]

	# 	if overall_d['regions_count'] == 0:
	# 		line += ['NA', "NA", 'NA']
	# 	else:
	# 		line += [
	# 			round(sum(overall_d['read_depths'])/overall_d['total_depth'], 4),
	# 			round(mean(overall_d['read_depths']),1),
	# 			median(overall_d['read_depths'])
	# 		]

	# 	print("\t".join(map(str, line)), file=outf)


	end_time = datetime.now()


	date_time = end_time.strftime("%Y/%m/%d, %H:%M:%S")
	elapsed = end_time - start_time
	print(f"Run completed: {date_time}  ({elapsed} elapsed)")	















