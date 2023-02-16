# Target predictions

import os
import click
from subprocess import PIPE, Popen, call
from pathlib import Path
from os.path import isfile, isdir
from pprint import pprint
from tqdm import tqdm
from multiprocessing import Lock, Process, Queue, current_process, Pool

from .generics import *
from .cli import cli



def call_job(job):

	s_name, s_seq = job

	lock.acquire()
	print(f"init\t{s_name}\t{s_seq}")
	lock.release()


	query_file = f"{s_name}.fa"

	with open(query_file, 'w') as outf:
		print(">", s_name, sep='', file=outf)
		print(s_seq, file=outf)


	z = len(s_seq) 


	call = ['RNAplex', '-t', cdna_file, '-q', query_file, '-f', rnaplex_f, '-z', str(z)]

	lock.acquire()
	print(f"call\t{s_name}", " ".join(call), sep='\t')
	lock.release()

	# sys.exit()
	p = Popen(call, stdout=PIPE, encoding='utf-8')

	out = []
	while True:
		
		mRNA_line = p.stdout.readline().lstrip(">").rstrip()
		if mRNA_line == "":
			break

		if mRNA_line != "(null)":
			m_name = mRNA_line
		p.stdout.readline()

		plex_line = p.stdout.readline().split()

		line = [m_name, s_name] + plex_line

		out.append(line)



	lock.acquire()
	print(f"out\t{s_name}\t{len(out)} hits found")

	with open(output_file, 'a') as outf:
		for line in out:
			print("\t".join(line), file=outf)
	lock.release()

	os.remove(query_file)


	lock.acquire()
	print(f"done\t{s_name}")
	lock.release()



def init(l, c, o, f):
	global lock
	global cdna_file
	global output_file
	global rnaplex_f

	lock = l
	cdna_file = c
	output_file = o
	rnaplex_f = f






@cli.command(group='Calculation', help_priority=6)

@click.option("-c", "--cDNA_file", 
	required=True, 
	type=click.Path(exists=True),
	help='File containing cDNA sequences in fasta format.')


@click.option("-o", "--output_directory", 
	required=True, 
	type=click.Path(),
	help="Directory name for annotation output.")


@click.option("-f", "--force",
	is_flag=True,
	help='Force remake of supporting files.')


@click.option("-m", "--max_rank",
	default=0,
	help='Maximum rank in abundance of a read within a locus (base 0).')


@click.option("--rnaplex_f",
	default='2',
	type=click.Choice(['0','1','2']),
	help='Speed variable for RNA plex, allowing options 0,2, and 1; ordered from slowest to fastest.')

@click.option("-t", "--threads",
	default=2,
	help='Number of thread processes for RNAplex. Each thread corresponds to an sRNA.')



def target(cdna_file, output_directory, force, max_rank, rnaplex_f, threads):
	"""Uses kmers and RNAplex to identify candidate targets in cDNA sequences."""



	output_directory = output_directory.rstrip("/")

	if not isdir(output_directory):
		sys.exit(f"output directory does not exist: {output_directory}")

	read_file = f"{output_directory}/Reads.txt"

	if not isfile(read_file):
		sys.exit(f"Read file {read_file} does not exist. Are you sure the output folder contains an annotation?")

	# results_file = f"{output_directory}/Results.txt"

	# if not isfile(results_file):
	# 	sys.exit(f"results file ({results_file}) not found")

	output_file = f"{output_directory}/Targeting.txt"

	with open(output_file, 'w') as outf:
		outf.write("")


	def get_kmers(seq, k):
		# print(seq)
		seq = deque(seq)

		kmer = deque()
		kmers = []

		while True:

			try:
				kmer[k-1]
				kmers.append("".join(kmer))

				if not seq:
					break
				kmer.popleft()

			except IndexError:
				pass
			kmer += seq.popleft()


		# for i,kmer in enumerate(kmers):
		# 	print(i * " " + kmer)
		# sys.exit()

		return(kmers)



	sRNAs = []
	with open(read_file, 'r') as f:
		header = f.readline()
		for line in f:
			line = line.strip().split("\t")

			cluster, seq, rank, _, _, _ = line
			rank = int(rank)

			if rank <= max_rank:
				sRNAs.append((cluster, seq))


			# sys.exit()



	mRNAs = []
	mRNA_lengths = {}

	with open(cdna_file, 'r') as f:
		header = f.readline()

		if header[0] != ">":
			print('ERROR: First line of cDNA file does not look like fasta:')
			print(header)
			sys.exit()

		seq = ''
		header = header[1:].strip()


		for line in f:
			line = line.strip()
			if line[0] == ">":
				seq = seq.upper().replace("T","U")

				mRNAs.append((header, seq))
				mRNA_lengths[header] = len(seq)

				seq = ''
				header = line[1:]


			else:
				seq += line


	# lock = Lock()
	# pool = Pool(initializer=init, 
	# 	initargs=(lock, cdna_file, output_file, rnaplex_f,), 
	# 	processes=threads)


	# pool.map(call_job, sRNAs)

	# pool.close()
	# pool.join()







	kmer_size = 6
	# kmer_window = 30
	kmer_window_expansion = 5
	kmer_threshold = 6

	kmer_d = {}

	window_d = {}

	# print(len(mRNAs))

	pbar = tqdm(total=len(mRNAs))
	for m_i, mRNA in enumerate(mRNAs):
		m_name, m_seq = mRNA
		# print(m_name)#, m_seq)
		pbar.update()

		complement_mRNA = m_seq[::-1]

		for k_i, kmer in enumerate(get_kmers(complement_mRNA, k=kmer_size)):

			# print(complement_mRNA)
			# print(kmer, m_name, k_i)

			try:
				kmer_d[kmer]
			except KeyError:
				kmer_d[kmer] = dict()

			try:
				kmer_d[kmer][m_name]
			except KeyError:
				kmer_d[kmer][m_name] = Counter()

			kmer_d[kmer][m_name][k_i] += 1

		if m_i > 1000:
			break


	pbar.close()
	# # pprint(kmer_d)

	# # sys.exit()

	theoretical_max = len(sRNAs) * m_i
	print(f'{theoretical_max:,} maximum number of queries by transcript')


	window_count = 0
	count = 0
	pbar = tqdm(total=len(sRNAs))
	for s_name, s_seq in sRNAs:
		pbar.update()

		s_kmers = get_kmers(s_seq, k=kmer_size)



		pos_d = dict()

		for s_kmer in s_kmers:

			try:
				keys = kmer_d[s_kmer].keys()
			except KeyError:
				keys = []

			for m_name in keys:

				try:
					pos_d[m_name]
				except KeyError:
					pos_d[m_name] = Counter()

				pos_d[m_name].update(kmer_d[s_kmer][m_name])

		kmer_window = len(s_seq)+kmer_window_expansion


		for m_name, c in pos_d.items():

			# pprint(c)

			# if m_name == 'Bcin03g05620.1' and s_name == 'pl_115':

			flag = ""
			for p in c.keys():

				window_hits = [j for j in c.keys() if p <= j <= p+kmer_window]


				if len(window_hits) > 6:

					count += 1
					# print(count, m_name)
					flag = "<---- " + str(len(window_hits))


					# print(s_name)
					# print("  ",m_name, c.keys(), flag)

					# print(window_hits)
					window_max = max(window_hits)
					window_min = min(window_hits)
					window_range = window_max - window_min
					offset = round((kmer_window - window_range)/2)
					window_max += offset
					window_min -= offset


					window = (m_name, mRNA_lengths[m_name] - window_max, mRNA_lengths[m_name] - window_min)
					# print(window)

					try:
						window_d[s_name].append(window)
					except KeyError:
						window_d[s_name] = [window]

					window_count += 1

					# sys.exit()



		# sys.exit()

	pbar.close()


	pprint(window_d)
	print(f'{len(window_d):,} sRNAs with windows')
	print(f'{window_count:,} windows in total')

	fold_reduction = theoretical_max / window_count
	print(f'{fold_reduction} fold reduction')

	sys.exit()

	# 		# sys.exit()



	# 	# print(len(c))
	# 	# pprint(c.most_common(10))

	# 	# print(len([v for v in c.values() if v > 3]))

	# 	# sys.exit()



	# sys.exit('done')




