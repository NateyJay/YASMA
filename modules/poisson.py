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
	half_window = int(window/2)

	chromosomes, bam_rgs = get_chromosomes(alignment_file, output_directory)
	annotation_readgroups = check_rgs(annotation_readgroups, bam_rgs)


	cl_i = 0
	merge_dist = 200
	pad = 20


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



	with open(f"{output_directory}/poisson.gff3", 'w') as outf:
		print("##gff-version 3", file=outf)

		for chrom, chrom_length in chromosomes:
			print(f"##sequence-region   {chrom} 1 {chrom_length}", file=outf)

	for chrom_count, chrom in enumerate(chromosomes):
		chrom, chrom_length = chrom

		print()
		print()
		print(f"{chrom_count+1} / {len(chromosomes)}")
		print(f"chrom: {chrom}")
		print(f"       {chrom_length:,} bp")

		sam_iter = samtools_view(alignment_file, rgs=annotation_readgroups, locus=chrom)


		depth_c  = Counter()
		# strand_c = Counter()

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
		for i, read in enumerate(sam_iter):

			if (i+1) % 10000 == 0:
				print(".", end='', flush=True)
				if (i+1) % 100000 == 0:
					print(" ", end='', flush=True)

					if (i+1) % 1000000 == 0:
						print(flush=True)
			# read_count += 1
			sam_strand, sam_length, sam_size, sam_pos, sam_chrom, sam_rg, sam_read, sam_name = read

			added_to_intergenic = False
			for r in range(sam_pos, sam_pos + sam_length + window):

				depth_c[r - half_window] += 1

				if not added_to_intergenic:
					if r in intergenic_positions:
						intergenic_read_count += 1
						added_to_intergenic = True

		print()


		lambda_chrom = window / intergenic_length * intergenic_read_count
		# lambda_chrom = window / chrom_length * read_count

		print("      ", round(lambda_chrom,4), 'lambda (chromosome, intergenic)')

		def poisson(k, l):
			p = l**k * math.e**(-l) / math.factorial(k)
			return(p)

		poisson_d = {}
		k = 0
		threshold = 0.00001
		while True:
			k += 1
			try:
				p = poisson(k,lambda_chrom)
			except OverflowError:
				break
			poisson_d[k] = p



		# pprint(poisson_d)





		in_locus = False
		for r in range(chrom_length):
			k = depth_c[r]
			# print(k)
			try:
				p = poisson_d[k]
			except KeyError:
				p = 0.0
			# print(r, k, round(p,5), sep='\t')

			if p <= 0.0001 and k > lambda_chrom:
				if not in_locus:
					start = r - pad
				in_locus = True
				stop = r + pad 


			else:
				if in_locus:
					if r > stop + merge_dist:

						cl_i += 1

						with open("test.gff3", 'a') as outf:
							print(chrom, 'poissonLocus','nc_RNA', start, stop, f'.\t.\t.\tID=Cl_{cl_i}', sep='\t', file=outf)

						in_locus = False



				# NC_037310.1	smoothLocus	nc_RNA	360214	360303	.	-	.	ID=Cl_1;dicercall=21;frac_dicercall=0.674


















		# depth_c = Counter()
		# for read in sam_iter:

		# 	# print(read)
		# 	sam_strand, sam_length, sam_size, sam_pos, sam_chrom, sam_rg, sam_read = read


		# 	for r in range(sam_length):
		# 		depth_c[r + sam_pos] += 1



		# # for r in range(chrom_length):
		# # 	print(depth_c[r])

		# 	# if sam_pos > 200:
		# 	# 	pprint(depth_c)
		# 	# 	sys.exit()


		# gene_depths = []
		# for start, stop in gene_d[chrom]:
		# 	gene_depths += [depth_c[r] for r in range(start, stop+1)]

		# intergenic_depths = []
		# intergenic_regions = []
		# last_stop = 0
		# for start, stop in gene_d[chrom]:
		# 	intergenic_regions.append((last_stop, start))
		# 	last_stop = stop

		# intergenic_regions.append((last_stop, chrom_length))


		# for start, stop in intergenic_regions:
		# 	intergenic_depths += [depth_c[r] for r in range(start, stop+1)]


		# n_quantiles=100
		# q_gene = quantiles(gene_depths, n=n_quantiles)
		# q_inter = quantiles(intergenic_depths, n=n_quantiles)
		# q_all  = quantiles([depth_c[r] for r in range(chrom_length)], n=n_quantiles)

		# for r in range(n_quantiles-1):
		# 	print(round(r/n_quantiles,3), q_gene[r], q_inter[r], q_all[r], sep='\t')

		# sys.exit()




































