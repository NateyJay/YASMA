# Locus quantification of all features for an input alignment

import os
import click
from subprocess import PIPE, Popen, call
from pathlib import Path
from os.path import isfile, isdir
from pprint import pprint
from multiprocessing import Lock, Process, Queue, current_process, Pool
from tqdm import tqdm


from modules.generics import *


def call_deepcount(job):

	name, locus = job

	# lock.acquire()
	# print(name, 'started', sep='\t')
	# lock.release()

	chrom = locus.split(":")[0]
	start = int(locus.split(":")[-1].split("-")[0])
	stop  = int(locus.split(":")[-1].split("-")[1])

	if start == 0:
		start = 1

	locus = f"{chrom}:{start}-{stop}"

	print(locus)
	c = Counter()

	sam_iter = samtools_view(alignment_file, locus=locus)

	for read in sam_iter:
		sam_strand, sam_length, _, sam_pos, _, sam_rg, sam_seq = read


		if sam_pos >= start and sam_pos + sam_length <= stop:
			print(read)

			c.update([(sam_rg, sam_length, sam_strand)])

	lock.acquire()
	# print(".", end='', flush=True)
	# print(name)

	for rg in rgs:
		for length in range(15,31):
			for strand in {"+", "-"}:

				# pprint(c)
				count = c[(rg, length, strand)]

				with open(output_file, 'a') as outf:
					print(name, rg, length, strand, count, sep='\t', file=outf)
	lock.release()

	# lock.acquire()
	# print(name, 'finished', sep='\t')
	# lock.release()


def init(l, r, a, o, ):
	global lock
	global rgs
	global alignment_file
	global output_file
	lock = l
	rgs = r
	alignment_file = a
	output_file = o

@click.command()

@click.option("-a", "--alignment_file", 
	required=True, 
	type=click.Path(exists=True),
	help='Alignment file input (bam or cram).')

@click.option("-o", "--output_directory", 
	default=f"Annotation_{round(time())}", 
	type=click.Path(),
	help="Directory name for annotation output")

@click.option("-f", "--force",
	is_flag=True,
	help='Force remake of supporting files')
@click.option("-t", "--threads",
	default=8,
	help="number of simultaneous threads for samtools_view.")


def deepcount(alignment_file, output_directory, force, threads):
	"""Gets counts for all readgroups, loci, strand, and sizes"""

	annotation_file = f"{output_directory}/Results.txt"
	output_file = f"{output_directory}/DeepCount.txt"

	with open(output_file, 'w') as outf:
		print("name\trg\tlength\tstrand\tabundance", file=outf)

	chroms, rgs = get_chromosomes(alignment_file, output_directory)
	jobs = []


	with open(annotation_file, 'r') as f:
		header = f.readline()
		# for i,h in enumerate(header.split('\t')):
		# 	print(i,h)
		# sys.exit()


		for line in f:
			line = line.strip().split("\t")

			name  = line[0]
			locus = line[1].replace("..","-")

			jobs.append((name, locus))

			

	# pprint(jobs)
	# sys.exit()

	jobs = jobs[:1]

	lock = Lock()
	pbar = tqdm(total=len(jobs))
	
	with Pool(initializer=init, 
		initargs=(lock, rgs, alignment_file, output_file,), 
		processes=threads) as p:


		for r in p.imap(call_deepcount, jobs):
			pbar.update()

	pbar.close()
		# r = list(tqdm(p.imap(call_deepcount, jobs), total=len(jobs)))

	# pool.map(call_deepcount, jobs)

	# pool.close()
	# pool.join()



