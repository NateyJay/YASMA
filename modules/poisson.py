# sRNA library background assessment

import sys
import click
from pathlib import Path
from os.path import isfile, isdir
from collections import Counter#, deque
from pprint import pprint

# import numpy as np
# from statistics import quantiles
import math

from modules.generics import *




class assessClass():
	'''produces a line assessment of a locus, similar to ShortStack3'''

	def __init__(self):

		self.header = ['Locus','Name','Length','Reads','RPM']
		self.header += ['UniqueReads','FracTop','Strand','MajorRNA','MajorRNAReads','Complexity']
		self.header += ['ma_size','ma_size_depth', 'ma_size_pair', 'ma_size_pair_depth', 'dicercall']


	def format(self, name, chrom, start, stop, reads, aligned_depth):



		### Basic information

		reads = [r for r in reads if not r[3] + r[1] < start or not r[3] > stop ]
		depth = len(reads)
		rpm = depth / aligned_depth * 1000000

		out = [f"{chrom}:{start}-{stop}", name, stop-start, depth, rpm]



		### ShortStack standard metrics

		seq_c       = Counter()
		strand_c    = Counter()
		size_c      = Counter()
		size_pair_c = Counter()

		for read in reads:
			sam_strand, sam_length, sam_size, sam_pos, sam_chrom, sam_rg, sam_read, sam_name = read

			seq_c[sam_read] += 1
			strand_c[sam_strand] += 1
			size_c[sam_length] += 1

			if sam_length > 15:
				size_pair_c[f"{sam_length-1}_{sam_length}"] += 1

			if sam_length > 30:
				size_pair_c[f"{sam_length}_{sam_length+1}"] += 1



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

		out += [unique_reads, frac_top, strand, major_rna, major_rna_depth, complexity]




		### More derived metrics

		predominant_size       = size_c.most_common()[0][0]
		predominant_size_depth = size_c.most_common()[0][1]

		try:
			predominant_size_pair       = size_pair_c.most_common()[0][0]
			predominant_size_pair_depth = size_pair_c.most_common()[0][1]
		except IndexError:
			predominant_size_pair       = "NA"
			predominant_size_pair_depth = 0

		if predominant_size_depth > depth / 2:
			dicercall = predominant_size

		elif predominant_size_pair_depth > depth / 2:
			dicercall = predominant_size_pair

		else:
			dicercall = "N"

		out += [predominant_size, predominant_size_depth, predominant_size_pair, predominant_size_pair_depth, dicercall]

		# print(out)
		# print('here!!!')
		return(out)

@click.command()

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
	help="Directory name for annotation output")

@click.option("-f", "--force",
	is_flag=True,
	help='Force remake of supporting files')

def poisson(alignment_file, annotation_readgroups, gene_annotation, output_directory, force):

	output_directory = output_directory.rstrip("/")

	window = 40
	cl_i = 0
	merge_dist = 200
	pad = 20

	half_window = int(window/2)

	chromosomes, bam_rgs = get_chromosomes(alignment_file, output_directory)
	annotation_readgroups = check_rgs(annotation_readgroups, bam_rgs)


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

	output_file = f"{output_directory}/Poisson.gff3"

	with open(output_file, 'w') as outf:
		print("##gff-version 3", file=outf)

		for chrom, chrom_length in chromosomes:
			print(f"##sequence-region   {chrom} 1 {chrom_length}", file=outf)


	out_results = f"{output_directory}/Poisson.txt"
	with open(out_results, 'w') as outf:
		print("\t".join(assessClass().header), file=outf)




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

		# read_count = 0

		sam_iter = samtools_view(alignment_file, rgs=annotation_readgroups, locus=chrom)
		for i, read in enumerate(sam_iter):
			total_read_count += 1


			if (i+1) % 10000 == 0:
				print(".", end='', flush=True)
				if (i+1) % 100000 == 0:
					print(" ", end='', flush=True)

					if (i+1) % 1000000 == 0:
						print(flush=True)
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

		lambda_chrom_rpm = lambda_chrom / i * 1000000

		# lambda_chrom = window / chrom_length * read_count

		print("      ", round(lambda_chrom,4), 'lambda (chromosome, intergenic)')

		print("      ", round(lambda_chrom_rpm,4), 'lambda_rpm')


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
					start = r - pad
					if start <= 0:
						start = 1
					region_c = Counter()


				in_locus = True
				stop = r + pad

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


						# region_length = stop - start

						# frac_top = region_c['+'] / (region_c['-'] + region_c['+'])

						# total = sum([region_c[l] for l in range(15,31)])

						# size_profile = [
						# 	sum([region_c[l] for l in range(15,20)])/total,
						# 	sum([region_c[l] for l in range(20,26)])/total,
						# 	sum([region_c[l] for l in range(26,31)])/total
						# ]


						locus_name = f"Cl_{cl_i}"

						with open(output_file, 'a') as outf:
							print(chrom, 'poissonLocus','nc_RNA', start, stop, f'.\t.\t.\tID={locus_name}', sep='\t', file=outf)


						loci.append((locus_name, chrom, start, stop))#, frac_top, size_profile, total/region_length))
						in_locus = False
						nucleated = False
				else:
					in_locus = False
					nucleated = False




		print(f"       {len(loci):,} loci found")
		print("       assessing...")


		for locus in loci:

			name, chrom, start, stop = locus
			coords = f"{chrom}:{start}-{stop}"


			reads = [r for r in samtools_view(alignment_file, locus=coords, rgs=annotation_readgroups)]

			# assess = 
			
			with open(out_results, 'a') as outf:
				out = assessClass().format(name, chrom, start, stop, reads, total_read_count)
				out = "\t".join(map(str, out))
				print(out, file=outf)

		all_loci += loci


	print()
	print(f"{len(all_loci):,} loci found in total")
		# sys.exit()























