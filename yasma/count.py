# Simple quantification of features

import os
import click
from subprocess import PIPE, Popen, call
from pathlib import Path
from os.path import isfile, isdir
from pprint import pprint
from multiprocessing import Lock, Process, Queue, current_process, Pool
# from tqdm.contrib.concurrent import process_map
from tqdm import tqdm

from .generics import *
from .cli import cli


def call_count(job):

	name, locus = job


	start = int(locus.split(":")[-1].split("-")[0])
	stop  = int(locus.split(":")[-1].split("-")[1])

	c = Counter()

	sam_iter = samtools_view(alignment_file, locus=locus)

	for read in sam_iter:
		sam_strand, sam_length, _, sam_pos, _, sam_rg, sam_seq = read


		if sam_pos >= start and sam_pos + sam_length <= stop:

			c.update([sam_rg])


	line = [name, locus, sum(c.values())]
	for rg in rgs:
		line.append(c[rg])

	# pprint(c)

	lock.acquire()
	# print(".", end='', flush=True)

	with open(output_file, 'a') as outf:
		print("\t".join(map(str, line)), file=outf)


	# pbar.update(1)
	lock.release()



def init(l, r, a, o, ):
	global lock
	global rgs
	global alignment_file
	global output_file
	lock = l
	rgs = r
	alignment_file = a
	output_file = o

@cli.command(group='Calculation', help_priority=3)

@click.option("-a", "--alignment_file", 
	required=False, 
	type=click.UNPROCESSED, callback=validate_path,
	help='Alignment file input (bam or cram).')

@click.option("-o", "--output_directory", 
	required=False,
	type=click.Path(),
	help="Directory name for annotation output")

@click.option("-l", "--locus_files", 
	required=False, 
	type=click.Path(exists=True),
	multiple=True,
	help='list of loci in 2-column format (name, chrom:start-stop). Multiples allowed. Defaults to the results file.')

@click.option("-g", "--gff_files", 
	required=False, 
	type=click.Path(exists=True),
	multiple=True,
	help='gff of genes to be evaluated. Looks for "mRNA" features as would be described in an NCBI annotation. Multiples allowed. Defaults to None.')

@click.option("--ignore_zeroes",
	is_flag=True,
	help="Don't save zero depth entries to deepcounts.txt. Saves file space, but can make it harder to analyze.")


def count(** params):
	"""Gets counts for all readgroups, loci, strand, and sizes."""

	rc = requirementClass()
	rc.add_samtools()
	rc.check()

	ic = inputClass(params)
	ic.check(['alignment_file'])

	output_directory     = str(ic.output_directory)
	alignment_file       = ic.inputs['alignment_file']

	locus_files          = params['locus_files']
	gff_files            = params['gff_files']
	ignore_zeroes        = params['ignore_zeroes']



	Path(output_directory+ "/counts/").mkdir(parents=True, exist_ok=True)


	counts_file = f"{output_directory}/counts/counts.txt"
	deep_counts_file = f"{output_directory}/counts/deepcounts.txt"

	if len(locus_files) == 0:
		locus_files = [f"{output_directory}/peak/loci.txt"]


	c = Counter()

	print('annotations:')
	loci = []
	for file in locus_files:
		with open(file, 'r') as f:
			header = f.readline()
			for line in f:
				c[file] += 1
				line = line.strip().split("\t")[:2]
				line[1] = line[1].replace("..", "-")
				loci.append((line[1], line[0]))

	for file in gff_files:
		with open(file, 'r') as f:
			for line in f:
				if not line.startswith("#"):
					line = line.strip().split("\t")
					if line[2] == 'gene':
						c[file] += 1
						coords = f"{line[0]}:{line[3]}-{line[4]}"
						name = line[8].split(";")[0].split("=")[1]

						loci.append((name, coords))

	for key, val in c.items():
		print(" ", key, "->", val, 'loci')

	print('')
	print('processing annotations...')


	chroms, rgs = get_chromosomes(alignment_file)

	# print('  -> making empty genome structure')
	# lookup = {}
	# for c,l in chroms:
	# 	lookup[c] = [None] * l
	# 	# lookup[c] = {}


	# print('  -> populating')
	# pbar = tqdm(total = len(loci), ncols=90)
	# for name, locus in loci:
	# 	pbar.update()
	# 	# print(locus)

	# 	chrom, start, stop = parse_locus(locus)
	# 	# print(chrom, start, stop)

	# 	for r in range(start, stop + 1):
	# 		try:
	# 			lookup[chrom][r].add(name)
	# 		except:
	# 			lookup[chrom][r] = set([name])
	# pbar.close()




	# depth_c = get_rg_depth(output_directory, alignment_file)

	# print(sum(depth_c.values()))
	# pprint(depth_c)
	# sys.exit()


	with open(counts_file, 'w') as outf:
		print('name','locus', "\t".join(rgs), sep='\t', file=outf)

	with open(deep_counts_file, 'w') as outf:
		print('name','rg','length','strand','count', sep='\t', file=outf)




	perc = percentageClass(1,len(loci))

	for i, locus in enumerate(loci):
		name, locus = locus

		c = Counter()
		deep_c = Counter()

		p_count = perc.get_percent(i)
		if p_count:
			print(f"   quantifying... {p_count}%", end='\r')

		chrom, start, stop = parse_locus(locus)

		for read in samtools_view(alignment_file, locus=locus):
			sam_strand, sam_length, sam_size, sam_pos, sam_chrom, sam_rg, sam_seq, sam_read_id = read


			c[sam_rg] += 1
			deep_c[(sam_rg, sam_length, sam_strand)] += 1


		line = [name, locus]

		line += [c[rg] for rg in rgs]

		with open(counts_file, 'a') as outf:

			print("\t".join(map(str, line)), file=outf)

		with open(deep_counts_file, 'a') as outf:
			for rg in rgs:
				for length in range(15,31):
					for strand in {"+", "-"}:
						count = deep_c[(rg, length, strand)]

						if count == 0 and ignore_zeroes:
							pass
						else:
							print(name, rg, length, strand, count, sep='\t', file=outf)


	print()

	# print()
	# print("processing alignments...")
	# print(f"  {sum(depth_c.values()):,} total")

	# pbar = tqdm(total = sum(depth_c.values()), ncols=90)

	# sam_iter = samtools_view(alignment_file)

	# for i,read in enumerate(sam_iter):

	# 	# if i % 1000000 == 0:
	# 	# 	print(read)
	# 	# 	break

	# 	sam_strand, sam_length, _, sam_pos, sam_chrom, sam_rg, sam_seq = read

	# 	if 15 <= sam_length <= 30:
	# 		pbar.update()


	# 		try:
	# 			start_loc = lookup[sam_chrom][sam_pos]
	# 		except IndexError:
	# 			start_loc = set([None])

	# 		try:
	# 			stop_loc  = lookup[sam_chrom][sam_pos + sam_length + 1]
	# 		except IndexError:
	# 			stop_loc = set([None])

	# 		# print(start_loc, stop_loc)

	# 		if start_loc and stop_loc:
	# 			names = list(start_loc & stop_loc)
	# 			names = [n for n in names if n]

	# 			# if len(names) == 0:
	# 			# 	names = ['NoLocus']

	# 		else:
	# 			# names = ['NoLocus']
	# 			names = []

	# 		for name in names:
	# 			c.update([(name, sam_rg)])
	# 			deep_c.update([(name, sam_rg, sam_length, sam_strand)])





		# if sum(c.values()) > 1000000:
		# 	pprint(c)
		# 	sys.exit()

	# pbar.close()

	# pprint(deep_c)
	# sys.exit()


	# print()
	# print('writing...')

	# with open(counts_file, 'w') as outf:
	# 	print('name','locus', "\t".join(rgs), sep='\t', file=outf)

	# with open(deep_counts_file, 'w') as outf:
	# 	print('name','rg','length','strand','count', sep='\t', file=outf)



	# loci.append(("NoLocus", "NA"))



	# pbar = tqdm(total = len(loci), ncols=90)
	# for name, locus in loci:
	# 	pbar.update()

	# 	line = [name, locus]

	# 	line += [c[(name, rg)] for rg in rgs]

	# 	with open(counts_file, 'a') as outf:
	# 		print("\t".join(map(str, line)), file=outf)


	# 	for rg in rgs:

	# 		for length in range(15,31):
	# 			for strand in {"+", "-"}:
	# 				count = deep_c[(name, rg, length, strand)]


	# 				with open(deep_counts_file, 'a') as outf:
	# 					print(name, rg, length, strand, count, sep='\t', file=outf)
	# 				# if count > 0:
	# 				# 	sys.exit()

	# pbar.close()
		# sys.exit()


	# with open(output_file, 'w') as outf:
	# 	header = "name\tlocus\ttotal"
	# 	for r in rgs:
	# 		header += f"\t{r}"
	# 	print(header, file=outf)




	# lock = Lock()

	# with Pool(initializer=init, 
	# 	initargs=(lock, rgs, alignment_file, output_file, ), 
	# 	processes=threads) as p:

	# 	r = list(tqdm(p.imap(call_count, loci), total=len(loci)))




