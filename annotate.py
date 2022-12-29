#!/usr/bin/env python3 

import sys
import os
from pprint import pprint

import click

from subprocess import PIPE, Popen, call
from pathlib import Path

from os.path import isfile, isdir

from time import time, sleep
from collections import Counter, deque
from itertools import count, chain

from statistics import median, mean



def get_chromosomes(file,output_directory):
	chromosomes = []
	rgs = []
	# call = f"samtools view -@ 4 -H {file}"
	call = ['samtools','view','-H', file]
	# print(call)

	p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')
	out, err = p.communicate()

	for o in out.strip().split("\n"):
		o = o.strip().split('\t')
		if o[0] == "@SQ":
			name = o[1].split(":")[-1]
			length = int(o[2].split(":")[-1])
			chromosomes.append((name,length))
		if o[0] == "@RG":
			rgs.append(o[1].split(":")[-1])


	with open(f"./{output_directory}/chrom.sizes.txt", 'w') as outf:
		for chrom, size in chromosomes:
			print(chrom, size, sep='\t', file=outf)

	return(chromosomes, rgs)




def check_rgs(annotation_readgroups, bam_rgs):

	if annotation_readgroups[0].lower() == 'all':
		annotation_readgroups = bam_rgs
	else:
		for rg in annotation_readgroups:
			if rg not in bam_rgs:
				sys.exit(f"Error: provided readgroup '{rg}' not found within bamfile header:\n{bam_rgs}")

	annotation_readgroups = set(annotation_readgroups)
	return(annotation_readgroups)
	# def get(self):


class Logger(object):
	def __init__(self, file_name):
		self.terminal = sys.stdout
		self.file_name = file_name
		self.log = open(file_name, "w")
		# with open(file_name, "w") as outf:
		# 	outf.write("")

	def clear_ansi(self, message):
		return(message.replace("\033[1m", "").replace("\033[0m",""))

	def write(self, message):
		self.terminal.write(message)
		# with open(self.file_name, 'a') as outf:
		# 	outf.write(message)  
		self.log.write(self.clear_ansi(message))

	def flush(self):
		self.terminal.flush()
		self.log.flush()



def inf_counter():
	i = 1
	while True:
		yield(i)
		i += 1





@click.group()
def cli():
	pass


# @click.command()
# @click.option('--count', default=1, help='Number of greetings.')
# @click.option('--name', prompt='Your name',
#               help='The person to greet.')
# def hello(count, name):
#     """Simple program that greets NAME for a total of COUNT times."""
#     for x in range(count):
#         click.echo(f"Hello {name}!")

@click.command()

@click.option("-a", "--alignment_file", 
	required=True, 
	type=click.Path(exists=True),
	help='Alignment file input (bam or cram).')

@click.option('-r', '--annotation_readgroups', 
	required=True,
	multiple=True,
	help="List of read groups (RGs, libraries) to be considered for the annotation. 'ALL' uses all readgroups for annotation, but often pertainent RGs will need to be specified individually.")

@click.option("-o", "--output_directory", 
	default=f"Annotation_{round(time())}", 
	type=click.Path(),
	help="Directory name for annotation output")

@click.option("-f", "--force",
	is_flag=True,
	help='Force remake of supporting files')

def precheck(alignment_file, annotation_readgroups, output_directory, force):
	"""Runs precheck to identify likely dicer sizes."""

	Path(output_directory).mkdir(parents=True, exist_ok=True)


	log_file = f"{output_directory}/Log_precheck.txt"

	message = f"log_file is already exists ({log_file}). The annotator will not over-write by default (use --force to override). Be warned: this will trigger the overwrite of some files in this folder!"
	assert not isfile(log_file) or force, message

	sys.stdout = Logger(log_file)


	chromosomes, bam_rgs = get_chromosomes(alignment_file, output_directory)
	annotation_readgroups = check_rgs(annotation_readgroups, bam_rgs)


	def get_rg_depth():
		depth_file = f"./{output_directory}/readgroup_depth.txt"
		# if isfile(depth_file): #and not force:
		# 	with open(depth_file, 'r') as f:
		# 		line = f.readline()
		# 		line = line.split("\t")

		# 		if line[0] == file:
		# 			# print(f'\nread depth from {depth_file}...')
		# 			return(int(line[1]))



		print('reading annotation RG depth...')

		c = Counter()

		call = ['samtools', 'view', '-F', '4']

		for r in annotation_readgroups:
			call += ['-r', r]

		call += [alignment_file]
		# print(call)

		p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')

		depth = 0
		for line in p.stdout:
			line = line.strip().split("\t")

			length = int(line[5][:-1])

			c.update([length])

			depth += 1

			# if depth > 1000000:
			# 	p.kill()
			# 	break

		p.wait()

		highest = max(c.values())

		print("length\tprop\tprop_highest\tabundance")
		for r in range(15,31):
			prop = round(c[r] / depth, 4)
			prop_highest = round(c[r] / highest, 4)
			print(r, prop, prop_highest, f"{c[r]:,}", sep='\t')

		with open(depth_file, 'w') as outf:
			print(alignment_file, depth, sep='\t', file=outf)

	get_rg_depth()






@click.command()
@click.option("-a", "--alignment_file", 
	required=True, 
	type=click.Path(exists=True),
	help='Alignment file input (bam or cram).')

@click.option('-r', '--annotation_readgroups', 
	required=True,
	multiple=True,
	help="List of read groups (RGs, libraries) to be considered for the annotation. 'ALL' uses all readgroups for annotation, but often pertainent RGs will need to be specified individually.")

@click.option("-d", "--dicercall", 
	default=[20,21,22,23,24],
	multiple=True, 
	help='List of sRNA lengths that derive from dicer.')

@click.option("-o", "--output_directory", 
	default=f"Annotation_{round(time())}", 
	type=click.Path(),
	help="Directory name for annotation output")

@click.option("-f", "--force",
	is_flag=True,
	help='Force remake of supporting files')

@click.option("--partial_wigs",
	is_flag=True,
	help='Only make wiggle files associated with essential functions (ignoring size and strand specific coverages. (May improve speed)')

@click.option("--window",
	default=100,
	help="Window size (centered on position) for determining DCR vs non-DCR read ratio (counting overlapping reads).")

@click.option("--merge_dist",
	default=150,
	help="Maximum gap size between valid regions to merge to a single locus.")

@click.option('--pad',
	default=10,
	help='Number of bases arbitrarily added to either end of a defined locus.')

@click.option('--rpm_cutoff',
	default=1.0,
	help='RPM depth threshold for DCR-sized reads to be considered as a valid region.')

@click.option('--extension_ratio',
	default=0.5,
	help='Fraction of RPM threshold to be considered for extending a locus boundaries')

@click.option('--dicer_ratio',
	default=3.0,
	help='Ratio of dicer to non-dicer reads to be considered for a valid region')

def annotate(alignment_file, annotation_readgroups, dicercall, output_directory, force, partial_wigs, window, merge_dist, pad, rpm_cutoff, extension_ratio, dicer_ratio):
	"""Main annotation suite."""
	print('run annotation')




	assert isfile(alignment_file), f"{alignment_file} does not exist"


	Path(output_directory).mkdir(parents=True, exist_ok=True)
	Path(f'./{output_directory}/coverages').mkdir(parents=True, exist_ok=True)




	log_file = f"{output_directory}/Log_annotation.txt"
	message = f"log_file is already exists ({log_file}). The annotator will not over-write by default (use --force to override). Be warned: this will trigger the overwrite of some files in this folder!"
	assert not isfile(log_file) or force, message

	if force and isfile(log_file):
		print("force flag included, overwrite possible", end='', flush=True)
		counter = 0
		while counter < 5:
			print('.', end='', flush=True)
			sleep(1)
			counter += 1
		

	sys.stdout = Logger(log_file)



	dicercall = [int(d) for d in dicercall]
	dcr_range = set([r for r in range(min(dicercall), max(dicercall) + 1)])
	non_range = set([r for r in range(15,30) if r not in dcr_range])


	assert window % 2 == 0, "Window must be an even number!"
	half_window = int(window / 2)



	# possible names:
	# DicerLocus
	# smallDicer


	print()
	print()
	print("\033[1m-- annotator v0.3x --\033[0m")

	print()
	print()
	print(f"\033[1m[Prerequisites]\033[0m")


	def check_reqs():
		tool_responses = {
		'samtools version' : 'Samtools compilation details:',
		# 'gt --version' : 'gt (GenomeTools)',
		# 'bgzip --version' : 'bgzip (htslib)',
		# 'tabix --version' : 'tabix (htslib)',
		'wigToBigWig' : 'wigToBigWig v 2.8',
		}


		fails = []

		for tool, response in tool_responses.items():
			p = Popen(tool.split(), stdout=PIPE, stderr=PIPE, encoding='utf-8')

			out, err = p.communicate()

			merged = out + err

			tool = tool.split()[0]


			# print(out)
			# print(err)
			if response in merged:
				pass_str = "[x]"
			else:
				pass_str = "[ ]"
				fails.append(tool)




			print(" ", pass_str, tool)
			# sys.exit()

		# do_not_prepare_gff = False
		do_not_make_bigwig = False

		if 'samtools' in fails:
			sys.exit("Error: samtools not found in PATH (required)")

		# for tool in ['gt','bgzip','tabix']:
		# 	if tool in fails:
		# 		do_not_prepare_gff = True
		# 		break

		if 'wigToBigWig' in fails:
			do_not_make_bigwig = True

		# if do_not_prepare_gff:
		# 	print("Warning: will not prepare indexed gff for jbrowse due to missing reqs")
		if do_not_make_bigwig:
			print("Warning: will not prepare bigwig files due to missing reqs")

		return(do_not_make_bigwig)


	do_not_make_bigwig = check_reqs()

	# rpm_cutoff = round(rpm_cutoff / window, 6)

	chromosomes, bam_rgs = get_chromosomes(alignment_file, output_directory)
	annotation_readgroups = check_rgs(annotation_readgroups, bam_rgs)





	## initiating output files
	gff_file   = f"{output_directory}/Annotation.gff3"

	with open(gff_file, 'w') as outf:
		print("##gff-version 3", file=outf)

		for c, l in chromosomes:
			print(f"##sequence-region   {c} 1 {l}", file=outf)


	count_file = f'{output_directory}/Counts.txt'
	with open(count_file, 'w') as outf:
		print("cluster", 'ann_depth', 'tot_depth', "\t".join(bam_rgs), sep='\t', file=outf)



	results_file = f"{output_directory}/Results.txt"
	with open(results_file, 'w') as outf:
		print("#name\tlocus\tlength\tgap\tdepth\trpm\tdepth:length\tfrac_top\tstrand\tfrac_dicer\tdcr_reads\tnon_reads\tdicercall\tfrac_dicercall\t" + "\t".join(map(str, dcr_range)), file=outf)


	reads_file = f"{output_directory}/TopReads.txt"
	with open(reads_file, 'w') as outf:
		print("cluster\tseq\trank\tdepth\trpm\tlocus_prop", file=outf)





	print()
	print(f"\033[1m[General settings]\033[0m")
	print(f"             alignment_file: {alignment_file}")
	print(f"      output_directory: {output_directory}")
	print(f" annotation_readgroups: {list(annotation_readgroups)}")
	print(f"           dicer_sizes: {list(dcr_range)}")
	print(f"                 force: {force}")
	print(f"          partial_wigs: {partial_wigs}")
	print(f"              log_file: {log_file}")
	print()


	print(f"\033[1m[Annotation settings]\033[0m")
	print(f"     window: {window}")
	print(f" merge_dist: {merge_dist}")
	print(f"        pad: {pad}")



	def get_library_depth(output_directory, file):
		depth_file = f"./{output_directory}/library_depth.txt"
		if isfile(depth_file): #and not force:
			with open(depth_file, 'r') as f:
				line = f.readline()
				line = line.split("\t")

				if line[0] == file:
					# print(f'\nread depth from {depth_file}...')
					return(int(line[1]))



		print('reading annotation RG depth...')

		call = ['samtools', 'view', '-F', '4']

		for rg in annotation_readgroups:
			call += ['-r', rg]

		call += [file]

		p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')

		depth = 0
		for line in p.stdout:
			depth += 1

		p.wait()


		with open(depth_file, 'w') as outf:
			print(file, depth, sep='\t', file=outf)

		return(depth)

	library_depth = get_library_depth(output_directory, alignment_file)


	read_equivalent = 1 / library_depth * 1000000
	depth_cutoff = library_depth / rpm_cutoff / 1000000
	ext_cutoff = rpm_cutoff * extension_ratio


	print()
	print('\033[1m[Depth settings]\033[0m')
	print(f'    ann_rg_depth: {library_depth:,} reads')
	print(f'          1 read: {round(read_equivalent,5)} rpm')
	print(f"      rpm_cutoff: {rpm_cutoff} rpm -> {round(depth_cutoff,2)} reads")
	print(f"      ext_cutoff: {ext_cutoff} rpm > {round(depth_cutoff*extension_ratio,2)} reads")
	print(f" extension_ratio: {extension_ratio}")
	print(f"     dicer_ratio: {dicer_ratio}")





	class wiggleClass():
		def __init__(self, file):
			self.file = f"./{output_directory}/coverages/{file}.wig"
			self.outf = open(self.file, 'w')
			self.reset()


		def reset(self):
			self.val = 0
			self.start_pos = 1


		def add(self, val, pos, chrom):


			if val != self.val:
				span = pos - self.start_pos

				if span > 0:

					print(f"variableStep chrom={chrom} span={span}", file=self.outf)
					print(f"{self.start_pos} {self.val}", file=self.outf)

					self.val = val
					self.start_pos = pos

		def convert(self, cleanup=False):

			self.outf.close()

			wig = self.file

			bigwig = wig.replace(".wig", ".bigwig")

			print(f"  {wig} -> {bigwig}", flush=True)

			call = f"wigToBigWig {wig} ./{output_directory}/chrom.sizes.txt {bigwig}"

			p = Popen(call.split(), stdout=PIPE, stderr=PIPE, encoding='utf-8')

			out, err= p.communicate()

			if out.strip() + err.strip() != "":

				print(out)
				print(err)

			if cleanup:
				os.remove(wig)

	coverage_names = ['dcr','non']

	if not partial_wigs:
		for s in ["+","-"]:
			for l in list(dcr_range) + ['non']:
				coverage_names.append(f"{l}{s}")



	wig_d = {c : wiggleClass(c) for c in coverage_names + ['rpm_passing', 'ratio_passing', 'passing_all']}


	cluster_counter = inf_counter()


	def samtools_view(bam):

		if not isfile(f"{bam}.bai"):
			# call = f"samtools index {bam}"
			call= ['samtools','index',bam]
			p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')
			out,err=p.communicate()
			print(out)
			print(err)


		# call = f"samtools view -@ 4 -F 4 {bam}"
		call = ['samtools', 'view', '-F', '4', bam]
		# print(call)
		p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')

		for i,line in enumerate(p.stdout):

			line = line.strip().split()

			# read_id, flag, sam_chrom, sam_pos, _, length, _, _, _, _,_,_,_,_,_,_,_,_,rg= line

			read_id = line[0]
			flag = line[1]
			seq = line[9]

			if flag == "16":
				strand = '-'
			elif flag == '0':
				strand = "+"
			else:
				strand = False

			# print(line)
			length = int(line[5].rstrip("M"))
			# sam_pos = int(sam_pos)

			# length = len(line[9])

			sam_pos = int(line[3])
			sam_chrom = line[2]

			rg = line[18].lstrip("RG:Z:")


			if length in dcr_range:
				size = 'dcr'
			elif length in non_range:
				size = 'non'
			else:
				size = False



			yield(strand, length, size, sam_pos, sam_chrom, rg, seq)

		p.wait()

	sam_iter = samtools_view(alignment_file)
	read = next(sam_iter)
	sam_strand, sam_length, sam_size, sam_pos, sam_chrom, sam_rg, sam_read = read


	total_locus_count = 0
	chrom_count = 0

	class locusClass():
		def __init__(self):
			self.reads = deque([])


			self.in_locus = False
			self.nucleated = False
			self.last_hit_pos = 0
			self.last_end = 0
			self.start = False
			self.stop  = False

		def hit(self, pos, hit_type):
			self.last_hit_pos = pos

			if not self.in_locus:
				self.start = pos

			self.in_locus = True


			if hit_type == 'nuc':
				self.nucleated = True


		def add(self, read):
			self.reads[-1].append(read)

		def check(self, pos):
			# print(self.reads)

			# if pos > 1392000:
			# 	print(pos, self.last_hit_pos)

			if self.in_locus:

				if pos - self.last_hit_pos > merge_dist:
					self.in_locus = False

					self.stop = self.last_hit_pos

					if self.nucleated:
						self.nucleated = False
						return(True)

			else:

				while True:

					if len(self.reads) == 0:
						break


					if len(self.reads[0]) == 0:
						self.reads.popleft()
					elif self.reads[0][0][0] + 35 < pos - merge_dist:
						self.reads.popleft()

					else:
						break
						
			self.reads.append([])






		def summarize(self, chrom):

			# start = self.reads[0][1]
			start = self.start# - self.pad
			stop  = self.stop# + self.pad


			# pprint(self.reads)
			# print(self.reads.values())

			reads = chain.from_iterable(self.reads)


			reads = [r for r in reads if r[0] + r[1] >= start and r[0] <= stop]

			len_c    = Counter()
			size_c   = Counter()
			strand_c = Counter()
			rg_c     = Counter()
			read_c   = Counter()

			read_starts = []
			read_stops  = []


			for r in reads:
				sam_pos, sam_length, sam_size, sam_strand, sam_rg, sam_read = r

				read_starts.append(sam_pos)
				read_stops.append(sam_pos + sam_length)

				len_c.update([sam_length])
				size_c.update([sam_size])
				strand_c.update([sam_strand])
				rg_c.update([sam_rg])
				read_c.update([sam_read])


			# print(start, stop)
			start = min(read_starts) - pad
			stop  = max(read_stops) + pad


			name = f"Cl_{next(cluster_counter)}"


			if len(reads) == 0:
				print(f"WARNING: {name} detected no reads. Likely an error. Skipping.")
				return(0,0)

			# pprint(len_c)
			# pprint(read_c)


			dist_to_last = start - self.last_end
			self.last_end = stop

			coords = f"{chrom}:{start}..{stop}"
			length = stop - start
			n_reads = len(reads)
			frac_top = round(strand_c["+"] / n_reads,3)
			frac_dicercall = round(size_c['dcr'] / n_reads, 3)
			rpm = round(n_reads / library_depth * 1000000, 2)


			cum_count = 0
			top_reads = read_c.most_common(100)
			with open(reads_file, 'a') as outf:
				for rank, read in enumerate(top_reads):

					seq, dep = read
					rpm = round(dep * read_equivalent, 4)

					cum_count += dep

					loc_prop = round(cum_count / n_reads, 4)

					print(name, seq, rank, dep, rpm, loc_prop, file=outf, sep='\t')

					if loc_prop >= 0.3:
						break




			predominant_length, predominant_length_depth = len_c.most_common(1)[0]
			predominant_length_depth = round(predominant_length_depth/n_reads,3)
			# print(predominant_length, predominant_length_depth)

			if frac_top >= 0.8:
				strand = '+'
			elif frac_top <= 0.2:
				strand = "-"
			else:
				strand = '.'

			depth_by_length = round(n_reads / length, 3)


			to_print = [name, coords, length, dist_to_last, n_reads, rpm, depth_by_length, frac_top, strand]
			to_print += [frac_dicercall, size_c['dcr'],  size_c['non']]
			to_print += [predominant_length, predominant_length_depth]
			to_print += [len_c[d] for d in dcr_range]
			to_print = "\t".join(map(str, to_print))


			with open(results_file, 'a') as outf:
				print(to_print, sep='\t', file=outf)

			# print(name, coords, length, dist_to_last, n_reads, rpm, depth_by_length, frac_top, strand, 
			# 		frac_dicercall, size_c['dcr'],  size_c['non'], 
			# 		predominant_length, round(predominant_length_depth/n_reads,3), "\t".join([str(len_c[d]) for d in dcr_range]), sep='\t')


			# sys.exit()

			with open(gff_file, 'a') as outf:
				print(f"{chrom}\tsmoothLocus\tnc_RNA\t{start}\t{stop}\t.\t.\t.\tID={name};dicercall={predominant_length};frac_dicercall={predominant_length_depth}", file=outf)


			to_print = [name]

			to_print.append(sum([rg_c[rg] for rg in annotation_readgroups]))
			to_print.append(sum([rg_c[rg] for rg in bam_rgs]))
			to_print += [rg_c[rg] for rg in bam_rgs]
			to_print = "\t".join(map(str, to_print))

			with open(count_file, 'a') as outf:
				print(to_print, file=outf)


			# self.reads = {}

			return(length, dist_to_last)

	class coverageClass():
		def __init__(self, bandwidth=0):
			self.ls = deque()
			self.bandwidth = bandwidth

		def get(self):
			try:
				d = self.ls.popleft()
			except IndexError:
				d = 0

			return(d)

		def add(self, length):
			for r in range(length + self.bandwidth + 1):
				try:
					self.ls[r] += 1
				except IndexError:
					self.ls.append(1)



	def test1(rpm):
		if rpm == None:
			return('-')

		if rpm >= rpm_cutoff:
			return('n')
		elif rpm >= ext_cutoff:
			return('e')
		else:
			return('-')

	def test2(dcr, non):
		if dcr == None or non == None:
			return('-')

		if dcr >= non * dicer_ratio and dcr > 0:
			return('x')
		else:
			return('-')


	for chrom, chrom_length in chromosomes:
		chrom_count += 1


		print()
		print()
		print(f"{chrom_count} / {len(chromosomes)}")
		print(f"chrom: {chrom}")
		print(f"       {chrom_length:,} bp")
		pos = 0

		print("  ", end='')

		locus_lengths = []
		locus_gaps = []



		coverages = {c : coverageClass() for c in coverage_names}
		coverage_buffer = {c : deque([0]*half_window) for c in coverage_names}

		window_coverages = {'dcr' : coverageClass(window), 'non' : coverageClass(window)}
		locus = locusClass()#merge_dist, chrom, pad, cluster_counter, library_depth, read_equivalent, dcr_range, annotation_readgroups, bam_rgs)



		while sam_chrom == chrom:

			# if pos == 1393300:
			# 	sys.exit()

			# if pos == 1000000:
			# 	sys.exit("timeup!")


			corrected_pos = pos - half_window


			if pos == sam_pos:

				if sam_size:

					if sam_rg in annotation_readgroups:

						coverages[sam_size].add(sam_length)
						window_coverages[sam_size].add(sam_length)

						if not partial_wigs:
							if sam_size == 'non':
								coverages[f'{sam_size}{sam_strand}'].add(sam_length)
							else:
								coverages[f'{sam_length}{sam_strand}'].add(sam_length)



					locus.add((sam_pos, sam_length, sam_size, sam_strand, sam_rg, sam_read))

				try:
					read = next(sam_iter)
					sam_strand, sam_length, sam_size, sam_pos, sam_chrom, sam_rg, sam_read = read
					# print(pos, sam_id, sam_size, sep='\t')
				except StopIteration:
					break



			elif pos < sam_pos:

				read_count = {}
				# dens_rpm = {}
				# dens = {}
				win_cov = {}
				rpms = {}



				for size in coverage_names:

					cov = coverages[size].get()
					coverage_buffer[size].append(cov)
					coverage_buffer[size].popleft()

					cov = coverage_buffer[size][0]

					if size[-1] == "-":
						cov = cov * -1




					wig_d[size].add(round(cov * read_equivalent,4), corrected_pos, chrom)

					if size in ['dcr','non']:
						win_cov[size] = window_coverages[size].get()
						rpms[size] = round(win_cov[size] * read_equivalent, 4)


				t1 = test1(round(coverage_buffer['dcr'][0] * read_equivalent, 4))
				t2 = test2(win_cov['dcr'], win_cov['non'])

				tests = f"{t1}{t2}"


				# print(tests)

				if t1 == 'n':
					wig_d['rpm_passing'].add(1, corrected_pos, chrom)

				elif t1 == 'e':
					wig_d['rpm_passing'].add(0.3, corrected_pos, chrom)

				else:
					wig_d['rpm_passing'].add(0, corrected_pos, chrom)


				if not win_cov['dcr'] or win_cov['dcr'] == 0:
					ratio = 0
				else:
					try:
						ratio = round(win_cov['dcr'] / win_cov['non'], 2)
					except ZeroDivisionError:
						ratio = dicer_ratio


				wig_d['ratio_passing'].add(ratio, corrected_pos, chrom)


				if tests == "nx":
					locus.hit(corrected_pos, 'nuc')
					wig_d['passing_all'].add(1, corrected_pos, chrom)

				elif tests == "n-" or t1 == 'e':
					locus.hit(corrected_pos, 'ext')
					wig_d['passing_all'].add(0.3, corrected_pos, chrom)

				else:
					wig_d['passing_all'].add(0, corrected_pos, chrom)



				# if tests and "-" not in tests:
				# if rds['dcr'] and rds['dcr'] > 0 and pos-half_window > 0:
				# if pos >

				# if pos > 1392000:
				# 	if "n" in tests or 'e' in tests:
				# 		print(chrom, pos-half_window, 
				# 			"||", coverage_buffer['dcr'][0], round(win_cov['dcr'], 4), 
				# 			"||", coverage_buffer['non'][0], round(win_cov['non'], 4),
				# 			'||', tests,
				# 			sep='\t')




				if locus.check(corrected_pos):
					length, gap = locus.summarize(chrom)

					locus_lengths.append(length)
					locus_gaps.append(gap)
					total_locus_count += 1


				pos += 1


				if pos % 100000 == 0:
					print(".", end='', flush=True)
				if pos % 1000000 == 0:
					print(" ", end='', flush=True)




		for key in wig_d.keys():
			wig_d[key].add(0, pos, chrom)
			wig_d[key].reset()
		# for size in ['dcr','non']:
		# 	wig_densities[size].add(0, pos)
		# wig_rpm_pass.add(0, pos)
		# wig_pass.add(0, pos)

		locus_count = len(locus_gaps)
		med_length = median(locus_lengths)
		med_gap = median(locus_gaps)
		print()
		print(f"  {locus_count:,} loci found")
		print(f"  {med_length:,} median length")
		print(f"  {med_gap:,} median gap")

		# break


	print()
	print(f"{total_locus_count:,} loci found in total")

	print()

	if not do_not_make_bigwig:
		print("converting wigs to bigwigs...")
		for key in wig_d.keys():
			wig_d[key].convert()
	else:
		print("Not making bigwig files due to missing req...")



	def prepare_gff(gff_input):

		sorted_input = gff_input.replace(".gff3", ".sorted.gff3")
		zipped_input = sorted_input.replace(".gff3", ".gff3.gz")


		# print(gff_input)
		# print(sorted_input)
		# print(zipped_input)



		print("  sorting...")
		c2 = f"gt gff3 -retainids -sortlines -tidy {gff_input}"
		with open(sorted_input, 'w') as f:
			c2 = Popen(c2.split(), stdout=f)

		print("  zipping...")
		c3 = f"bgzip -f {sorted_input}"
		call(c3.split())

		print("  indexing...")
		c4 = f"tabix -f -p gff {zipped_input}"
		call(c4.split())


	# print()
	# if not do_not_prepare_gff:
	# 	print("indexing GFF for jbrowse...")
	# 	prepare_gff(gff_file)
	# else:
	# 	print("Not preparing indexed gff due to missing reqs...")




@click.command()
@click.option("-a", "--alignment_file", 
	required=True, 
	type=click.Path(exists=True),
	help='Alignment file input (bam or cram).')

@click.option('-r', '--annotation_readgroups', 
	required=True,
	multiple=True,
	help="List of read groups (RGs, libraries) to be considered for the annotation. 'ALL' uses all readgroups for annotation, but often pertainent RGs will need to be specified individually.")

@click.option("-o", "--output_directory", 
	default=f"Annotation_{round(time())}", 
	type=click.Path(),
	help="Directory name for annotation output")

def fold(alignment_file, annotation_readgroups, output_directory):



	results_file = f"{output_directory}/Results.txt"

	assert isfile(results_file), f"results_file {results_file} not found... (Have you run annotation with this directory?)"

	with open(results_file, 'r') as f:

		for line in f:
			line = line.strip().split()
			print(line)
			# sys.exit()







cli.add_command(precheck)
cli.add_command(annotate)
cli.add_command(fold)

if __name__ == '__main__':
	cli()








