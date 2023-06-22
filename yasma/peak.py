

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
from scipy.stats import kruskal, kstest


class dcClass():
	def __init__(self, sizes=[]):
		self.size_d = {1 : Counter(), 2 : Counter(), 3 : Counter()}

		self.depth = 0

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

		self.size_key_d = dict()
		for n in [1,2,3]:
			self.size_key_d[n] = {}
			for size in range(15,31):
				keys = get_size_keys(size, n)
				self.size_key_d[n][size] = keys

		if len(sizes) > 0:
			self.update(sizes)



	def update(self, sizes):
		for size in sizes:
			self.depth += 1

			if 15 <= size <= 30:

				self.size_d[1][size] += 1
				self.size_d[2].update(self.size_key_d[2][size])
				self.size_d[3].update(self.size_key_d[3][size])

	def get(self):

		def get_most_common(c):
			mc = c.most_common(1)

			if mc == []:
				return("NA", 0)

			else:
				return(mc[0][0], mc[0][1])



		self.size_1_key, self.size_1_depth = get_most_common(self.size_d[1])
		self.size_2_key, self.size_2_depth = get_most_common(self.size_d[2])
		self.size_3_key, self.size_3_depth = get_most_common(self.size_d[3])


		if self.size_1_depth > self.depth * 0.5:
			dicercall = self.size_1_key

		elif self.size_2_depth > self.depth * 0.5 and self.depth > 30:
			dicercall = self.size_2_key

		elif self.size_3_depth > self.depth * 0.5 and self.depth > 60:
			dicercall = self.size_3_key

		else:
			dicercall = "N"

		self.dicercall = str(dicercall)

		return(dicercall)

	def __str__(self):
		return(str(self.get()))

	def __eq__(self, other):

		def process_dcall(dcall, s3k):

			if dcall == "N":
				return(set(["N"]))


			keys = set(str(s3k).split("_"))

			if dcall.count("_") == 2:
				keys.add("N")

			return(keys)

		self_keys  = process_dcall(self.dicercall,  self.size_3_key)
		other_keys = process_dcall(other.dicercall, other.size_3_key)

		common = self_keys.intersection(other_keys)

		if len(common) > 0:
			return True
		else:
			return False



class assessClass():
	'''produces a line assessment of a locus, similar to ShortStack3'''

	def __init__(self):

		self.header = ['Locus','Name','Length','Reads','RPM']
		self.header += ['UniqueReads','FracTop','Strand','MajorRNA','MajorRNAReads','Complexity']
		self.header += ['Gap', 'size_1n','size_1n_depth', 'size_2n','size_2n_depth', 'size_3n','size_3n_depth', 'dicercall']





	def format(self, name, chrom, start, stop, reads, aligned_depth, last_stop):


		# size_key_d = dict()
		# for n in [1,2,3]:
		# 	size_key_d[n] = {}
		# 	for size in range(15,31):
		# 		keys = get_size_keys(size, n)
		# 		size_key_d[n][size] = keys



		### Basic information

		input_len = len(reads)

		reads = [r for r in reads if not r[3] + r[1] < start or not r[3] > stop ]
		depth = len(reads)
		rpm = depth / aligned_depth * 1000000





		### ShortStack standard metrics

		seq_c       = Counter()
		strand_c    = Counter()
		# size_d      = {1 : Counter(), 2 : Counter(), 3 : Counter()}

		dicercall = dcClass()

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

				dicercall.update([sam_length])

				# size_d[1][sam_length] += 1
				# size_d[2].update(size_key_d[2][sam_length])
				# size_d[3].update(size_key_d[3][sam_length])

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


		



		frac_top   = round(frac_top,3)
		complexity = round(complexity,3)
		rpm        = round(rpm,3)

		dicercall.get()

		result_line = [f"{chrom}:{start}-{stop}", name, stop-start, depth, rpm]
		result_line += [unique_reads, frac_top, strand, major_rna, major_rna_depth, complexity]
		result_line += [
			gap, 
			dicercall.size_1_key, dicercall.size_1_depth,
			dicercall.size_2_key, dicercall.size_2_depth,
			dicercall.size_3_key, dicercall.size_3_depth,
			dicercall.dicercall
		]


		if dicercall.dicercall == 'N':
			feature_type = "OtherRNA"
		else:
			# size_range = len(dicercall.split("_"))
			feature_type = f"RNA_{dicercall}"

		gff_line = [
			chrom, 'yasma_empirical',feature_type, start, stop, '.', strand, '.',
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
	required=False, 
	type=click.Path(exists=True),
	help='Gene annotation file in gff3 format. Tested with NCBI annotation formats.')

@click.option("-o", "--output_directory", 
	# default=f"Annotation_{round(time())}", 
	required=True,
	type=click.Path(),
	help="Directory name for annotation output.")

@click.option("-f", "--force",
	is_flag=True,
	help='Force remake of supporting files.')

# @click.option("-s", "--sensitivity",
# 	default=0.95,
# 	help="Sensitivity threshold for discovering sRNA-producing regions. Must be 0 < S <= 1.0")

# @click.option("--window",
# 	default=40,
# 	help="Window size (centered on position) for counting reads.")

@click.option("--rpm_threshold",
	default=0.5,
	help="Depth threshold in reads per million for discovery of a locus peak.")

@click.option("--clump_dist",
	default=500,
	help="")

@click.option("--creep_dist",
	default=50,
	help="")

@click.option('--pad',
	default=10,
	help='Number of bases arbitrarily added to either end of a defined locus. Default 10 nt.')

@click.option('--peak_threshold',
	default=0.05,
	help='Minimum depth, as a proportion of the maximum coverage depth in a locus, used for trimming low-depth edges from a locus. Default 0.05 percent of peak.')

@click.option('--kernel_window',
	default=0.40,
	help='')


@click.option('--coverage_method',
	default='kernel',
	type=click.Choice(['kernel', 'depth']),
	help='')

@click.option('--boundary_method',
	default='creep',
	type=click.Choice(['creep', 'bool']),
	help='')



def peak(alignment_file, annotation_readgroups, gene_annotation, output_directory, force, rpm_threshold, clump_dist, creep_dist, pad, peak_threshold, kernel_window, coverage_method, boundary_method):
	'''Annotator based on empircally-derived probability scores.'''

	# window = 40
	# clump_dist = 500
	# creep_dist = 50


	output_directory = output_directory.rstrip("/")

	Path(output_directory+ "/Empirical/").mkdir(parents=True, exist_ok=True)

	log_file = f"{output_directory}/Empirical_log.txt"

	sys.stdout = Logger(log_file)


	cl_i = 0


	# assert 0 < sensitivity <= 1.0, "sensitivity must be 0 < S <= 1.0"
	assert 0 < peak_threshold < 1, "peak_trim must be between 0 and 1."

	# half_window = int(window/2)

	chromosomes, bam_rgs = get_chromosomes(alignment_file, output_directory)
	annotation_readgroups = check_rgs(annotation_readgroups, bam_rgs)

	# chromosomes = [c for c in chromosomes if c[0] == 'NC_037320.1']

	chrom_depth_c = get_global_depth(output_directory, alignment_file, aggregate_by=['rg','chrom'])

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






	## iterating through chromosomes and reads

	all_loci = []
	total_read_count = 0
	cl_i = 0


	total_reads = sum(chrom_depth_c.values())
	read_equivalent = 1 / sum(chrom_depth_c.values()) * 1000000



	for chrom_count, chrom_and_length in enumerate(chromosomes):



		loci = []
		chrom, chrom_length = chrom_and_length

		print()
		print()
		print(f"{chrom_count+1} / {len(chromosomes)}")
		print(f"chrom: {chrom}")
		print(f"       {chrom_length:,} bp")
		print(f"       {chrom_depth_c[chrom]:,} reads")


		def get_peaks(bam, rgs, chrom):

			depth_c = Counter()

			perc = percentageClass(1, chrom_length)


			c1 = ['samtools', 'view', '-h', '-F', '4']
			for rg in rgs:
				c1 += ['-r', rg]
			c1 += [bam, chrom]

			p1 = Popen(c1, stdout=PIPE, stderr=PIPE, encoding='utf-8')

			c2 = ['samtools', 'depth', '-a', '-']
			p2 = Popen(c2, stdin=p1.stdout, stdout=PIPE, stderr=PIPE, encoding='utf-8')


			
			for i,line in enumerate(p2.stdout):

				_, pos, depth = line.strip().split('\t')

				depth_c[int(pos)] = int(depth)


				perc_out = perc.get_percent(i)
				if perc_out:
					print(f"   reading position depths... {perc_out}%", end='\r', flush=True)

			p2.wait()

			print()

			return(depth_c)



		def get_kernel_peaks(bam, rgs, chrom):

			half_window = math.floor(kernel_window/2)

			depth_c = Counter()

			perc = percentageClass(1, chrom_depth_c[chrom])


			call = ['samtools', 'view', '-F', '4']
			for rg in rgs:
				call += ['-r', rg]
			call += [bam, chrom]

			p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')

			for i,line in enumerate(p.stdout):
				line = line.strip().split("\t")

				pos    = int(line[3])
				length = int(line[5].rstrip("M"))

				pos += math.floor(length / 2)
				depth_c[pos] += 1

				perc_out = perc.get_percent(i)
				if perc_out:
					print(f"   reading position depths ..... {perc_out}%", end='\r', flush=True)

			p.wait()
			print()


			## Performing kernel density smoothing based on mean or median or sum?

			kernel_c = Counter()
			window_dq = deque()

			perc = percentageClass(1, chrom_length)
			for i in range(chrom_length):

				window_dq.append(depth_c[i])

				kernel_c[i] = round(sum(window_dq),2)
				



				perc_out = perc.get_percent(i)
				if perc_out:
					print(f"   calculating kernel density .. {perc_out}%", end='\r', flush=True)

				if len(window_dq) > kernel_window:
					window_dq.popleft()

			print()

			return(kernel_c)






		## Evaluating peaks to form loci

		if coverage_method == "depth":
			peak_c = get_peaks(
				bam=alignment_file, 
				rgs=annotation_readgroups, 
				chrom=chrom)

		elif coverage_method == "kernel":
			peak_c = get_kernel_peaks(
				bam=alignment_file, 
				rgs=annotation_readgroups, 
				chrom=chrom)

		claim_d = {}
		loci = []


		def boolean_creep(direction, cursor, depth_threshold):

			positions = deque()
			depths    = deque()

			while True:

				cursor += direction

				if cursor in claim_d:
					break

				positions.append(cursor)
				depths.append(peak_c[cursor])

				if len(positions) > pad:
					positions.popleft()
					depths.popleft()

					window_average = max(depths)

					if window_average <= depth_threshold:
						break


			boundary = cursor

			return(boundary)

		def creep(direction, cursor, depth_threshold):

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
				depths.append(peak_c[cursor])

				if len(positions) > creep_dist:
					positions.popleft()
					depths.popleft()

					window_average = mean(depths)

					if window_average <= depth_threshold:
						break

			if len(positions) <= 1:
				return(cursor)


			while True:
				# print("positions", positions)
				# print("depths", depths)
				positions.pop()
				depths.pop()
				
				window_median = median(depths)

				# print(window_median)

				if window_median > depth_threshold or len(depths) == 1:
					break

			# print(cursor)


			boundary = math.ceil(median(positions)) + (direction * pad)

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


		candidate_peak_count = len([v for v in peak_c.values() if v / total_reads * 1000000 > rpm_threshold])





		perc = percentageClass(5, candidate_peak_count)

		i = 0

		for center, depth in peak_c.most_common():


			perc_out = perc.get_percent(i)
			if perc_out:
				print(f"   forming peaks to loci ....... {perc_out}%", end='\r', flush=True)
			i += 1


			if center not in claim_d:

				rpm = depth / total_reads * 1000000


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

					loci.append([claim, chrom, lbound, rbound])

					# print()
					# print(f"{lbound:,} <- {center:,} -> {rbound:,}")
					# print(loci[-1])
					# print(rbound - lbound, "nt")
					# input()

		print()






		## Sorting loci by position

		loci.sort(key=lambda x: x[2])





		## Clumping neighbor loci

		def get_basic_dimensions(l):
			name, chrom, start, stop = l
			coords = f"{chrom}:{start}-{stop}"


			reads = [r for r in samtools_view(alignment_file, locus=coords, rgs=annotation_readgroups)]

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





		i = 0

		unclumped_loci_count = len(loci)
		clump_events = 0
				
		print(f"   clumping similar neighbors... {unclumped_loci_count} -> {len(loci)} loci", end='\r', flush=True)

		while i < (len(loci)-1):

			loc1 = loci[i]

			if loc1[0] == "Cluster_36":
				print(loci[i-1])
				print(loci[i-2])
				print(loci[i-3])

			j = 0 
			clump_j = 0
			while i+j+1 < len(loci):

				loc2 = loci[i+j+1]



				if loc2[2] - loc1[3] > clump_dist:


					for r in range(loc1[2], loc1[3]+1):
						claim_d[r] = loc1[0]

					for r in range(clump_j+1, 0,-1):
						del loci[i+r]

					break


				else:


					ftop1, len1 = get_basic_dimensions(loc1)
					ftop2, len2 = get_basic_dimensions(loc2)


					# c1 = Counter(len1)
					# c2 = Counter(len2)

					# most_common_p1 = c1.most_common(1)[0][1] / sum(c1.values())
					# most_common_p2 = c2.most_common(1)[0][1] / sum(c2.values())



					dc1 = dcClass(len1)
					dc2 = dcClass(len2)



					dc1.get()
					dc2.get()


					if loc1[0] == "Cluster_36":
						print()
						print()
						print("i:", i)
						print(loc1, "->", loc2)
						print(dc1)
						print(dc2)
						print(dc1 == dc2)
						print(ftop1, ftop2)
						print(abs(ftop1 - ftop2) < 0.5)
						print("merge?", abs(ftop1 - ftop2) < 0.5 and dc1 == dc2)

						print()
						print()
						print(loci[i+2])
						input()


					# if 0.5 < len_error < 2:
					if abs(ftop1 - ftop2) < 0.5 and dc1 == dc2:

						# print(loc1, dep1, ftop1)
						# print(loc2, dep2, ftop2)
						# print('merge!!')
						# print()


						print(f"   clumping similar neighbors .. {unclumped_loci_count} -> ", end='', flush=True)


						loc1[3] = loc2[3]
						loci[i] = loc1




						print(f"{len(loci)} loci", end='\r', flush=True)

						clump_events += 1
						clump_j = j

				j += 1
			i += 1

		print()



		# for i in range(1356816, 1360359):

		# 	try:
		# 		d = peak_c[i]
		# 	except KeyError:
		# 		d = 0

		# 	try:
		# 		c = claim_d[i]
		# 	except KeyError:
		# 		c = None

		# 	rpm = round(read_equivalent * d, 2)

		# 	print(i, d, rpm, rpm >= rpm_threshold, c, sep="\t")


		# print(rpm_threshold)

		# sys.exit()




		## Assessing locus dimensions and making annotations

		read_depths = []

		perc = percentageClass(increment=5, total=len(loci))


		last_stop = 0
		for i,locus in enumerate(loci):


			print_percentage = perc.get_percent(i)
			if print_percentage:
				print(f"   assessing loci .............. {print_percentage}%", end="\r", flush=True)

			name, chrom, start, stop = locus
			coords = f"{chrom}:{start}-{stop}"


			reads = [r for r in samtools_view(alignment_file, locus=coords, rgs=annotation_readgroups)]

			read_depths.append(len(reads))

			read_c = Counter()
			for read in reads:
				sam_strand, _, _, _, _, _, sam_read, _ = read

				if sam_strand == '-':
					sam_read = complement(sam_read[::-1])

				read_c[sam_read] += 1



			
			results_line, gff_line = assessClass().format(name, chrom, start, stop, reads, sum(chrom_depth_c.values()), last_stop)


			last_stop = stop

			if results_line:
				with open(results_file, 'a') as outf:
					print("\t".join(map(str, results_line)), file=outf)

				with open(gff_file, 'a') as outf:
					print("\t".join(map(str, gff_line)), file=outf)

				top_reads_save(read_c, reads_file, read_equivalent, name)

		all_loci += loci
		# sys.exit()

		print()
		print()

		print(f"   {len(loci):,} loci found")
		print(f"      from {clump_events} clumping events")

		mean_length   = mean([l[3]-l[2] for l in loci])
		median_length = median([l[3]-l[2] for l in loci])

		print(f"     length:")
		print(f"       mean ---> {round(mean_length,1)} bp")
		print(f"       median -> {median_length} bp")


		mean_depth   = mean(read_depths)
		median_depth = median(read_depths)

		print(f"     abundance:")
		print(f"       mean ---> {round(mean_depth,1)} reads ({round(mean_depth * read_equivalent, 2)} rpm)")
		print(f"       median -> {median_depth} reads ({round(median_depth * read_equivalent, 2)} rpm)")

	print()
















