

import sys
import click
from pathlib import Path
from os.path import isfile, isdir
from collections import Counter#, deque
from pprint import pprint
from random import sample

import numpy as np
# from statistics import quantiles
import math

from .generics import *
from .cli import cli

from statistics import mean, median






class assessClass():
	'''produces a line assessment of a locus, similar to ShortStack3'''

	def __init__(self):

		self.header = ['Locus','Name','Length','Reads','RPM']
		self.header += ['UniqueReads','FracTop','Strand','MajorRNA','MajorRNAReads','Complexity']
		self.header += ['Gap', 'size_1n','size_1n_depth', 'size_2n','size_2n_depth', 'size_3n','size_3n_depth', 'dicercall']





	def format(self, name, chrom, start, stop, reads, aligned_depth, last_stop):

		def get_size_keys(size, n, min=15, max=30):

			if n == 1:
				return(size)

			out = []

			for r in range(n):
				key = []

				for k in range(n):
					key.append(str(size-n+r+1+k))

				out.append("_".join(key))
			return(out)

		size_key_d = dict()
		for n in [1,2,3]:
			size_key_d[n] = {}
			for size in range(15,31):
				keys = get_size_keys(size, n)
				size_key_d[n][size] = keys



		### Basic information

		input_len = len(reads)

		reads = [r for r in reads if not r[3] + r[1] < start or not r[3] > stop ]
		depth = len(reads)
		rpm = depth / aligned_depth * 1000000





		### ShortStack standard metrics

		seq_c       = Counter()
		strand_c    = Counter()
		size_d      = {1 : Counter(), 2 : Counter(), 3 : Counter()}

		for read in reads:
			sam_strand, sam_length, sam_size, sam_pos, sam_chrom, sam_rg, sam_read, sam_name = read



			if 15 <= sam_length <= 30:

				# print(sam_length)
				# print(self.get_size_keys(sam_length, 1))
				# print(self.get_size_keys(sam_length, 2))
				# print(self.get_size_keys(sam_length, 3))
				# print(self.get_size_keys(sam_length, 4))
				# sys.exit()



				seq_c[sam_read] += 1
				strand_c[sam_strand] += 1

				size_d[1][sam_length] += 1
				size_d[2].update(size_key_d[2][sam_length])
				size_d[3].update(size_key_d[3][sam_length])

		if sum(strand_c.values()) == 0:
			# this solution is a stop-gap. Apparently, in some species low-lambdas can lead to trace-loci being annotated, which only have reads outside of accepted ranges. 
			return(None,None)



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

		def get_most_common(c):
			mc = c.most_common(1)

			if mc == []:
				return("NA", 0)

			else:
				return(mc[0][0], mc[0][1])



		size_1_key, size_1_depth = get_most_common(size_d[1])
		size_2_key, size_2_depth = get_most_common(size_d[2])
		size_3_key, size_3_depth = get_most_common(size_d[3])


		if size_1_depth > depth * 0.5:
			dicercall = size_1_key

		elif size_2_depth > depth * 0.5 and depth > 30:
			dicercall = size_2_key

		elif size_3_depth > depth * 0.5 and depth > 60:
			dicercall = size_3_key

		else:
			dicercall = "N"


		frac_top   = round(frac_top,3)
		complexity = round(complexity,3)
		rpm        = round(rpm,3)

		result_line = [f"{chrom}:{start}-{stop}", name, stop-start, depth, rpm]
		result_line += [unique_reads, frac_top, strand, major_rna, major_rna_depth, complexity]
		result_line += [
			gap, 
			size_1_key, size_1_depth,
			size_2_key, size_2_depth,
			size_3_key, size_3_depth,
			dicercall
		]


		if dicercall == 'N':
			feature_type = "OtherRNA"
		else:
			# size_range = len(dicercall.split("_"))
			feature_type = f"RNA_{dicercall}"

		gff_line = [
			chrom, 'yasma_empirical',feature_type, start, stop, '.', strand, '.',
			f'ID={name};dicercall={dicercall};depth={depth};rpm={rpm};fracTop={frac_top};majorRNA={major_rna}'
		]


		return(result_line, gff_line)

class locusClass():
	'''keeps track of coordinates and a list of loci'''

	def __init__(self, threshold, merge_dist, chrom_length, chrom, depth_c):
		self.in_locus     = False
		self.threshold    = threshold
		self.merge_dist   = merge_dist
		self.chrom_length = chrom_length
		self.chrom        = chrom
		self.cluster_i    = 0

		self.loci = []

		self.depth_c = depth_c
		# self.alignment_file  = alignment_file
		# self.rgs             = rgs
		# self.total_reads     = total_reads
		# self.rpm_threshold   = rpm_threshold


	def add(self, depth, position):
		if depth > self.threshold:
			# print(depth, self.threshold, position, sep='\t')
			

			if not self.in_locus:
				self.start = position
				if self.start <= 0:
					self.start = 1


			self.in_locus = True
			self.stop = position


		else:

			if self.in_locus:
				if position > self.stop + self.merge_dist:
					self.cluster_i += 1

					if self.stop > self.chrom_length:
						self.stop = self.chrom_length

					# if simple_assessment(self):

					self.refine_boundaries()

					locus_name = f"pl_{self.cluster_i}"

					self.loci.append([locus_name, self.chrom, self.start, self.stop])
					self.in_locus = False


	def refine_boundaries(self):
		boundary_window = 20
		boundary_threshold = 0.1

		start, stop = self.start, self.stop

		depths = []
		for r in range(start, stop+1):
			depths.append(self.depth_c[r])

		positional_median = median(depths)
		print(depths)
		print(positional_median)


		start_in = start
		# print(start, "<- start in")
		while True:

			left_median = median([self.depth_c[r] for r in range(start, start + boundary_window + 1)])
			print(left_median)

			if left_median <= boundary_threshold * positional_median:
				break

			start -= 1
		# print(start, "<- start out")

		print(f"{start_in} <- {start} ({start-start_in} nt)")



		stop_in = stop
		# print(stop, "<- stop in")
		while True:

			right_median = median([self.depth_c[r] for r in range(stop - boundary_window, stop+1)])

			if right_median <= boundary_threshold * positional_median:
				break

			stop += 1
		# print(stop, "<- stop out")

		print(f"{stop} -> {stop_in} ({stop - stop_in} nt)")



		self.start = start
		self.stop  = stop

		input()







# def simple_assessment(self):
# 	coords = f"{self.chrom}:{self.start}-{self.stop}"
# 	sam_iter = samtools_view(self.alignment_file, rgs=self.rgs, locus=coords)

# 	for i,read in enumerate(sam_iter):
# 		sam_strand, sam_length, sam_size, sam_pos, sam_chrom, sam_rg, sam_read, sam_name = read
# 		# print(sam_name)


# 	rpm = (i+1) / self.total_reads * 1000000

# 	if rpm / self.rpm_threshold:
# 		return(False)




@cli.command(group='Annotation', help_priority=1)

@click.option("-a", "--alignment_file", 
	required=True, 
	type=click.Path(exists=True),
	help='Alignment file input (bam or cram).')

@click.option('-r', '--annotation_readgroups', 
	required=True,
	multiple=True,
	help="List of read groups (RGs, libraries) to be considered for the annotation. 'ALL' uses all readgroups for annotation, but often pertainent RGs will need to be specified individually.")

@click.option("-g", "--gene_annotation", 
	required=False, 
	type=click.Path(exists=True),
	help='Gene annotation file in gff3 format. Tested with NCBI annotation formats.')

@click.option("-o", "--output_directory", 
	default=f"Annotation_{round(time())}", 
	type=click.Path(),
	help="Directory name for annotation output.")

@click.option("-f", "--force",
	is_flag=True,
	help='Force remake of supporting files.')

@click.option("-s", "--sensitivity",
	default=0.95,
	help="Sensitivity threshold for discovering sRNA-producing regions. Must be 0 < S <= 1.0")

@click.option("--window",
	default=40,
	help="Window size (centered on position) for counting reads.")

@click.option("--merge_dist",
	default=150,
	help="Maximum gap size between valid regions to merge to a single locus. Default 150 nt.")

@click.option('--pad',
	default=10,
	help='Number of bases arbitrarily added to either end of a defined locus. Default 10 nt.')

@click.option('--peak_trim',
	default=0.05,
	help='Minimum depth, as a proportion of the maximum coverage depth in a locus, used for trimming low-depth edges from a locus. Default 0.05 percent of peak.')


def annotate(alignment_file, annotation_readgroups, gene_annotation, output_directory, force, sensitivity, window, merge_dist, pad, peak_trim):
	'''Annotator based on empircally-derived probability scores.'''

	output_directory = output_directory.rstrip("/")

	Path(output_directory+ "/Empirical/").mkdir(parents=True, exist_ok=True)

	log_file = f"{output_directory}/Empirical_log.txt"

	sys.stdout = Logger(log_file)


	cl_i = 0


	assert 0 < sensitivity <= 1.0, "sensitivity must be 0 < S <= 1.0"
	assert 0 < peak_trim < 1, "peak_trim must be between 0 and 1."

	half_window = int(window/2)

	chromosomes, bam_rgs = get_chromosomes(alignment_file, output_directory)
	annotation_readgroups = check_rgs(annotation_readgroups, bam_rgs)


	chrom_depth_c = get_global_depth(output_directory, alignment_file, aggregate_by=['rg','chrom'])

	pprint(chrom_depth_c)
	sys.exit()

	keys = list(chrom_depth_c.keys())
	for key in keys:
		if key[0] in annotation_readgroups:
			chrom_depth_c[key[1]] += chrom_depth_c[key]

		del chrom_depth_c[key]



	## preparing output files

	gff_file = f"{output_directory}/Annotation.gff3"

	with open(gff_file, 'w') as outf:
		print("##gff-version 3", file=outf)

		for chrom, chrom_length in chromosomes:
			print(f"##sequence-region   {chrom} 1 {chrom_length}", file=outf)


	results_file = f"{output_directory}/Results.txt"
	with open(results_file, 'w') as outf:
		print("\t".join(assessClass().header), file=outf)


	reads_file = f"{output_directory}/Reads.txt"
	with open(reads_file, 'w') as outf:
		print(TOP_READS_HEADER, file=outf)


	assess_file = f"{output_directory}/Assessment.txt"
	with open(assess_file, 'w') as outf:
		print('chrom\tthreshold\tcount\tmedian_length\tmean_length\tmedian_gap\tmean_gap', file=outf)


	## iterating through chromosomes and reads

	all_loci = []
	total_read_count = 0
	cl_i = 0


	total_reads = sum(chrom_depth_c.values())
	read_equivalent = 1 / sum(chrom_depth_c.values()) * 1000000

	rpm_threshold = 0.5


	for chrom_count, chrom_and_length in enumerate(chromosomes):



		loci = []
		chrom, chrom_length = chrom_and_length

		print()
		print()
		print(f"{chrom_count+1} / {len(chromosomes)}")
		print(f"chrom: {chrom}")
		print(f"       {chrom_length:,} bp")
		print(f"       {chrom_depth_c[chrom]:,} reads")


		depth_c  = Counter()
		detail_d = dict()

		read_position_d = dict()
		read_identity_d = dict()


		max_t = 36


		## Iterating through reads for a chromosome and counting window depths

		perc = percentageClass(1, chrom_depth_c[chrom])

		sam_iter = samtools_view(alignment_file, rgs=annotation_readgroups, locus=chrom)
		for i, read in enumerate(sam_iter):
			total_read_count += 1

			perc_out = perc.get_percent(i)
			if perc_out:
				print(f"   reading alignment... {perc_out}%", end='\r', flush=True)


			sam_strand, sam_length, sam_size, sam_pos, sam_chrom, sam_rg, sam_read, sam_name = read

			key = sam_pos + int(sam_length/2)



			# added_to_intergenic = False
			for r in range(sam_pos, sam_pos + sam_length + window):

				depth_c[r - half_window] += 1

				# try:
				# 	read_position_d[r].add(sam_name)
				# except:
				# 	read_position_d[r] = set([sam_name])




		print()

		# print("      ", round(threshold,4), 'reads / window (chromosome, whole)')


		print('  writing window depths to file...')
		with open(f"{output_directory}/{chrom}.txt", 'w') as outf:
			for r in range(chrom_length):
				print(r, depth_c[r], sep='\t', file=outf)



		sys.exit()
		print("   parsing to loci...")

		loci_d = {}
		for t in range(max_t):
			loci_d[t] = locusClass(t, merge_dist, chrom_length, chrom, depth_c)


		# in_locus = False

		# for t in range(max_t):
		# 	print(t)

		for t in [5]:
			for r in range(chrom_length):

				loci_d[t].add(depth=depth_c[r], position=r)




		for t in range(max_t):
			print(t, len(loci_d[t].loci), sep='\t')







		sys.exit()

		# print(f"       {len(loci):,} loci found")






		# perc = percentageClass(increment=5, total=len(loci))


		# last_stop = 0
		# for i,locus in enumerate(loci):


		# 	print_percentage = perc.get_percent(i)
		# 	if print_percentage:
		# 		print(f"   assessing loci... {print_percentage}%", end="\r", flush=True)

		# 	name, chrom, start, stop = locus
		# 	coords = f"{chrom}:{start}-{stop}"


		# 	reads = [r for r in samtools_view(alignment_file, locus=coords, rgs=annotation_readgroups)]

		# 	read_c = Counter()
		# 	for read in reads:
		# 		sam_strand, _, _, _, _, _, sam_read, _ = read

		# 		if sam_strand == '-':
		# 			sam_read = complement(sam_read[::-1])

		# 		read_c[sam_read] += 1



			
		# 	results_line, gff_line = assessClass().format(name, chrom, start, stop, reads, sum(chrom_depth_c.values()), last_stop)


		# 	last_stop = stop

		# 	if results_line:
		# 		with open(results_file, 'a') as outf:
		# 			print("\t".join(map(str, results_line)), file=outf)

		# 		with open(gff_file, 'a') as outf:
		# 			print("\t".join(map(str, gff_line)), file=outf)

		# 		top_reads_save(read_c, reads_file, read_equivalent, name)

		# all_loci += loci
		# # sys.exit()


	# print()
	# print()
	# print('chromosome', 'lambda', 'loci', sep='\t')
	# for c,l in chromosomes:
	# 	print(c, round(output_d['lambda'][c]), output_d['loci_count'][c], sep='\t')


	print()
	print(f"{len(all_loci):,} loci found in total")
		# sys.exit()
















