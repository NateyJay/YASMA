# Simple quantification of features

import os
import click
from click_option_group import optgroup


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

	contig = locus.split(":")[0]
	start = int(locus.split(":")[-1].split("-")[0])
	stop  = int(locus.split(":")[-1].split("-")[1])

	c = Counter()

	sam_iter = samtools_view(alignment_file, contig=contig, start=start, stop=stop)

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


@optgroup.group('\n  Basic options',
				help='')

@optgroup.option("-a", "--alignment_file", 
	required=False, 
	type=click.UNPROCESSED, callback=validate_path,
	help='Alignment file input (bam or cram).')

@optgroup.option("-o", "--output_directory", 
	required=True,
	type=click.Path(),
	help="Directory name for annotation output")

@optgroup.option("-c", "--conditions", 
	required=False, 
	multiple=True,
	type=click.UNPROCESSED, callback=validate_condition,
	help='Values denoting condition groups (sets of replicate libraries) for projects with multiple tissues/treatments/genotypes. Can be entered here as space sparated duplexes, with the library base_name and condition groups delimited by a colon. E.g. SRR1111111:WT SRR1111112:WT SRR1111113:mut SRR1111114:mut')


@optgroup.option("-a", "--annotation_file", 
	required=False, 
	type=click.UNPROCESSED, callback=validate_path,
	default = "tradeoff/loci.gff3",
	multiple=False,
	help='Locus annotation in gff, gff3, bed, or tabular format. Tabular requires contig:start-stop and locus_namein the first two columns (tab-delimited, "#" escape char). Defaults to tradeoff/loci.gff3, but that may not be prefered if alignment includes conditions.')

@optgroup.option("-n", "--name", 
	required=False,
	type=str,
	help="Optional name. Useful if comparing annotations.")



@optgroup.option("--include_zeros",
	is_flag=False,
	help="Include to save 0-depth rows in the deep counts. By default, these are excluded to save space (except for one entry to make sure downstream analyses will include un-found entries)")


def count(** params):
	"""Gets counts for all readgroups, loci, strand, and sizes."""

	rc = requirementClass()
	rc.check()

	ic = inputClass(params)
	ic.check(['alignment_file'])

	output_directory     = str(ic.output_directory)
	alignment_file       = ic.inputs['alignment_file']
	conditions           = ic.inputs['conditions']

	annotation_file      = Path(params['annotation_file'])
	include_zeros        = params['include_zeros']
	name                 = params['name']

	Path(output_directory, "counts").mkdir(parents=True, exist_ok=True)

	if name:
		counts_file      = Path(output_directory, 'counts', f'{name}_counts.txt')
		deep_counts_file = Path(output_directory, 'counts', f'{name}_deepcounts.txt') 

	else:
		counts_file      = Path(output_directory, 'counts', 'counts.txt')
		deep_counts_file = Path(output_directory, 'counts', 'deepcounts.txt') 


	c = Counter()


	def process_annotation(file):
		with open(file, 'r') as f:
			if file.suffix == '.txt':
				f.readline()

			for i,line in enumerate(f):
				if not line.startswith("#"):
					line = line.strip().split("\t")

					if file.suffix == ".txt":
						coords, name = line[:2]
						coords = coords.replace("..", '-')

					elif file.suffix == '.gff' or file.suffix == '.gff3':
						coords = f"{line[0]}:{line[3]}-{line[4]}"
						name   = line[8].split(";")[0].split("=")[-1]

					elif file.suffix == '.bed':
						name   = f"bed_{i}"
						coords = f"{line[0]}:{line[1]}-{line[2]}"

					else:
						print(f'file.suffix "{file.suffix}" not expected in annotation file...')
						sys.exit()

					yield (name, coords)




	print('annotations:')

	loci = []


	for name, coords in process_annotation(annotation_file):
		# c[file] += 1
		loci.append((name, coords))



	# for file in locus_files:
	# 	with open(file, 'r') as f:
	# 		header = f.readline()
	# 		for line in f:
	# 			c[file] += 1
	# 			line = line.strip().split("\t")[:2]
	# 			line[1] = line[1].replace("..", "-")
	# 			loci.append((line[1], line[0]))

	# for file in gff_files:
	# 	with open(file, 'r') as f:
	# 		for line in f:
	# 			if not line.startswith("#"):
	# 				line = line.strip().split("\t")
	# 				if line[2] == 'gene':
	# 					c[file] += 1
	# 					coords = f"{line[0]}:{line[3]}-{line[4]}"
	# 					name = line[8].split(";")[0].split("=")[1]

	# 					loci.append((name, coords))

	# for key, val in c.items():
	# 	print(" ", key, "->", val, 'loci')

	print('')
	print('processing annotation...')


	chroms, rgs = get_chromosomes(alignment_file)

	pprint(conditions)

	rev_conditions = {}
	for k,vs in conditions.items():
		for v in vs:
			rev_conditions[v] = k


	cond_rgs = [f"{rev_conditions[r]}.{r}" for r in rgs]

	with open(counts_file, 'w') as outf:
		print('name','locus', "\t".join(cond_rgs), sep='\t', file=outf)

	with open(deep_counts_file, 'w') as outf:
		print('name', 'condition', 'rg','length','strand','count', sep='\t', file=outf)


	missed_rg     = set(rgs)
	missed_length = set(range(15,31))
	missed_strand = set(["-","+"])

	perc = percentageClass(1,len(loci))

	for i, locus in enumerate(loci):

		name, locus = locus

		c = Counter()
		deep_c = Counter()

		p_count = perc.get_percent(i)
		if p_count:
			print(f"   quantifying... {p_count}%", end='\r')

		chrom, start, stop = parse_locus(locus)

		for read in samtools_view(alignment_file, contig=chrom, start=start, stop=stop):
			sam_strand, sam_length, sam_size, sam_pos, sam_chrom, sam_rg, sam_seq, sam_read_id = read


			c[sam_rg] += 1
			deep_c[(sam_rg, sam_length, sam_strand)] += 1



		line = [name, locus]

		line += [c[rg] for rg in rgs]

		with open(counts_file, 'a') as outf:

			print("\t".join(map(str, line)), file=outf)

		with open(deep_counts_file, 'a') as outf:
			for rg in rgs:
				try:
					cond = rev_conditions[rg]
				except KeyError:
					cond = 'None'

				for length in range(15,31):
					for strand in {"+", "-"}:
						count = deep_c[(rg, length, strand)]

						if count > 0 or include_zeros:
							print(name, cond, rg, length, strand, count, sep='\t', file=outf)

							try:
								missed_strand.remove(strand)
							except KeyError:
								pass

							try:
								missed_rg.remove(rg)
							except KeyError:
								pass
	
							try:
								missed_length.remove(length)
							except KeyError:
								pass
						
					

	with open(deep_counts_file, 'a') as outf:
		for rg in missed_rg:
			cond = rev_conditions[rg]
			print(name, cond, rg, '15', '-', '0', sep='\t', file=outf)
		for strand in missed_strand:
			cond = rev_conditions[rg]
			print(name, cond, rg, '15', strand, '0', sep='\t', file=outf)
		for length in missed_length:
			print(name, cond, rg, length, '-', '0', sep='\t', file=outf)






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




