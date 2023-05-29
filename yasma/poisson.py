# sRNA library background assessment

import sys
import click
from pathlib import Path
from os.path import isfile, isdir
from collections import Counter#, deque
from pprint import pprint

import numpy as np
# from statistics import quantiles
import math

from .generics import *
from .cli import cli



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
			chrom, 'yasma_poisson',feature_type, start, stop, '.', strand, '.',
			f'ID={name};dicercall={dicercall};depth={depth};rpm={rpm};fracTop={frac_top};majorRNA={major_rna}'
		]


		return(result_line, gff_line)

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
	required=True, 
	type=click.Path(exists=True),
	help='Gene annotation file in gff3 format. Tested with NCBI annotation formats.')

@click.option("-o", "--output_directory", 
	default=f"Annotation_{round(time())}", 
	type=click.Path(),
	help="Directory name for annotation output.")

@click.option("-f", "--force",
	is_flag=True,
	help='Force remake of supporting files.')

@click.option("--window",
	default=40,
	help="Window size (centered on position) for counting reads, used as k in the poisson model. Default 40 nt.")

@click.option("--merge_dist",
	default=150,
	help="Maximum gap size between valid regions to merge to a single locus. Default 150 nt.")

@click.option('--pad',
	default=10,
	help='Number of bases arbitrarily added to either end of a defined locus. Default 10 nt.')

@click.option('--peak_trim',
	default=0.05,
	help='Minimum depth, as a proportion of the maximum coverage depth in a locus, used for trimming low-depth edges from a locus. Default 0.05 percent of peak.')


def poisson(alignment_file, annotation_readgroups, gene_annotation, output_directory, force, window, merge_dist, pad, peak_trim):
	'''Annotator based on poisson-derived probability scores.'''

	output_directory = output_directory.rstrip("/")

	Path(output_directory).mkdir(parents=True, exist_ok=True)

	# window = 40
	cl_i = 0
	# merge_dist = 200
	# pad = 20
	# peak_trim = 0.05


	assert 0 < peak_trim < 1, "peak_trim must be between 0 and 1."

	half_window = int(window/2)

	chromosomes, bam_rgs = get_chromosomes(alignment_file, output_directory)
	annotation_readgroups = check_rgs(annotation_readgroups, bam_rgs)


	chrom_depth_c = get_global_depth(output_directory, alignment_file, aggregate_by=['rg','chrom'])

	keys = list(chrom_depth_c.keys())
	for key in keys:
		if key[0] in annotation_readgroups:
			chrom_depth_c[key[1]] = chrom_depth_c[key]

		del chrom_depth_c[key]


	# pprint(chrom_depth_c)
	# sys.exit()


	read_equivalent = 1 / sum(chrom_depth_c.values()) * 1000000


	def calc_poisson(k, l):
		p = l**k * math.e**(-l) / math.factorial(k)
		return(p)


	gene_d = {}
	with open(gene_annotation, 'r') as f:
		for line in f:
			line = line.strip().split("\t")

			feature_type = line[2]

			if feature_type == "mRNA":

				chrom = line[0]
				start, stop = [int(l) for l in line[3:5]]

				try:
					gene_d[chrom].append((start, stop))
				except KeyError:
					gene_d[chrom] = [(start, stop)]



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




	## iterating through chromosomes and reads

	all_loci = []
	total_read_count = 0
	cl_i = 0


	for chrom_count, chrom in enumerate(chromosomes):

		loci = []
		chrom, chrom_length = chrom

		print()
		print()
		print(f"{chrom_count+1} / {len(chromosomes)}")
		print(f"chrom: {chrom}")
		print(f"       {chrom_length:,} bp")
		print(f"       {chrom_depth_c[chrom]:,} reads")


		depth_c  = Counter()
		detail_d = dict()

 		## Defining parts of the genome which are not from annotated mRNAs

		intergenic_length = 0
		intergenic_read_count = 0
		intergenic_positions = set()
		last_stop = 0
		for start, stop in gene_d[chrom]:
			intergenic_length += stop - start
			for r in range(start, stop):
				intergenic_positions.add(r)
			last_stop = stop


		## Iterating through reads for a chromosome and counting window depths

		perc = percentageClass(1, chrom_depth_c[chrom])

		sam_iter = samtools_view(alignment_file, rgs=annotation_readgroups, locus=chrom)
		for i, read in enumerate(sam_iter):
			total_read_count += 1

			perc_out = perc.get_percent(i)
			if perc_out:
				print(f"   reading alignment... {perc_out}%", end='\r', flush=True)


			# if (i+1) % 10000 == 0:
			# 	print(".", end='', flush=True)
			# 	if (i+1) % 100000 == 0:
			# 		print(" ", end='', flush=True)

			# 		if (i+1) % 1000000 == 0:
			# 			print(flush=True)
			# read_count += 1
			sam_strand, sam_length, sam_size, sam_pos, sam_chrom, sam_rg, sam_read, sam_name = read

			key = sam_pos + int(sam_length/2)
			try:
				detail_d[key][sam_length] += 1
				detail_d[key][sam_strand] += 1
			except KeyError:
				detail_d[key] = Counter()
				detail_d[key][sam_length] += 1
				detail_d[key][sam_strand] += 1


			added_to_intergenic = False
			for r in range(sam_pos, sam_pos + sam_length + window):

				depth_c[r - half_window] += 1

				if not added_to_intergenic:
					if r in intergenic_positions:
						intergenic_read_count += 1
						added_to_intergenic = True

		print()


		lambda_chrom = window / intergenic_length * intergenic_read_count

		# lambda_chrom_rpm = lambda_chrom / i * 1000000

		# lambda_chrom = window / chrom_length * read_count

		print("      ", round(lambda_chrom,4), 'lambda (chromosome, intergenic)')

		# print("      ", round(lambda_chrom_rpm,4), 'lambda_rpm')


		poisson_d = {}
		k = 0
		# threshold = 0.00001
		while True:
			k += 1
			try:
				p = calc_poisson(k,lambda_chrom)
			except OverflowError:
				break
			poisson_d[k] = p


		# pprint(poisson_d)

		print("   parsing to loci...")
		loci = []

		in_locus = False
		nucleated = False
		for r in range(chrom_length):

			k = depth_c[r]

			try:
				p = poisson_d[k]
			except KeyError:
				p = 0.0



			if k > lambda_chrom:
				if p <= 0.00001:
					nucleated = True


				if not in_locus:
					start = r
					if start <= 0:
						start = 1
					region_c = Counter()



				in_locus = True
				stop = r 

				try:
					region_c.update(detail_d[r])
				except KeyError:
					pass



			else:

				if in_locus and nucleated:
					if r > stop + merge_dist:
						cl_i += 1

						if stop > chrom_length:
							stop = chrom_length



						depth_profile = [depth_c[i] for i in range(start, stop+1)]
						trim_threshold = max(depth_profile) * peak_trim


						# pprint(depth_profile)
						# print(trim_threshold, 'trim_threshold')

						for from_left, depth_k in enumerate(depth_profile):
							if depth_k > trim_threshold:
								break

						for from_right, depth_k in enumerate(depth_profile[::-1]):
							if depth_k > trim_threshold:
								break

						# print(len(depth_profile), 'len(depth_profile)')

						# print(stop - start, "stop-start")

						# print(start, stop, 'initial')

						start += from_left
						stop  -= from_right
						# print(start, stop, 'trimmed (95%)')

						start -= pad
						stop  += pad



						locus_name = f"pl_{cl_i}"

						loci.append((locus_name, chrom, start, stop))#, frac_top, size_profile, total/region_length))
						in_locus = False
						nucleated = False

						# if locus_name == "Cl_418":
						# 	sys.exit()
				else:
					in_locus = False
					nucleated = False



		print(f"       {len(loci):,} loci found")


		perc = percentageClass(increment=5, total=len(loci))


		last_stop = 0
		for i,locus in enumerate(loci):


			print_percentage = perc.get_percent(i)
			if print_percentage:
				print(f"   assessing loci... {print_percentage}%", end="\r", flush=True)

			name, chrom, start, stop = locus
			coords = f"{chrom}:{start}-{stop}"


			reads = [r for r in samtools_view(alignment_file, locus=coords, rgs=annotation_readgroups)]

			read_c = Counter()
			for read in reads:
				sam_strand, _, _, _, _, _, sam_read, _ = read

				if sam_strand == '-':
					sam_read = complement(sam_read[::-1])

				read_c[sam_read] += 1



			
			results_line, gff_line = assessClass().format(name, chrom, start, stop, reads, sum(chrom_depth_c.values()), last_stop)

			last_stop = stop

			with open(results_file, 'a') as outf:
				print("\t".join(map(str, results_line)), file=outf)

			with open(gff_file, 'a') as outf:
				print("\t".join(map(str, gff_line)), file=outf)

			top_reads_save(read_c, reads_file, read_equivalent, name)

		all_loci += loci
		# sys.exit()

	print()
	print(f"{len(all_loci):,} loci found in total")
		# sys.exit()
















