# Target predictions

import os
import click
from subprocess import PIPE, Popen, call
from pathlib import Path
from os.path import isfile, isdir
from pprint import pprint
from tqdm import tqdm

from .generics import *
from .cli import cli

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


@click.option("--max_rank",
	default=0,
	help='Maximum rank position of reads in library to be tested for targets (base-0).')


@click.option("--kmer_size",
	default=6,
	help='k-mer size for zeroing-in on targets.')


@click.option("--kmer_window",
	default=30,
	help='Window size used for candidate targets.')


def target(cdna_file, output_directory, force, max_rank, kmer_size, kmer_window):
	"""Uses kmers and RNAplex to identify candidate targets in cDNA sequences."""


	half_kmer = round(kmer_size/2)

	output_directory = output_directory.rstrip("/")

	if not isdir(output_directory):
		sys.exit(f"output directory does not exist: {output_directory}")

	read_file = f"{output_directory}/Reads.txt"

	if not isfile(read_file):
		sys.exit(f"Read file {read_file} does not exist. Are you sure the output folder contains an annotation?")

	# results_file = f"{output_directory}/Results.txt"

	# if not isfile(results_file):
	# 	sys.exit(f"results file ({results_file}) not found")




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

				seq = ''
				header = line[1:]


			else:
				seq += line



	kmer_d = {}

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
	# pprint(kmer_d)

	# sys.exit()


	count = 0
	pbar = tqdm(total=len(sRNAs))
	for s_name, s_seq in sRNAs:
		pbar.update()

		s_kmers = get_kmers(s_seq, k=kmer_size)



		pos_d = dict()

		for s_kmer in s_kmers:

			# print(s_kmer)

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

				# print(kmer_d[kmer][m_name])
				# print(pos_d)
				# sys.exit()

			# pprint(pos_d)
			# sys.exit()

		# 	print(" ", len(pos_d))

		# print(len(pos_d))
		# sys.exit()


		for m_name, c in pos_d.items():

			# pprint(c)

			# if m_name == 'Bcin03g05620.1' and s_name == 'pl_115':

			flag = ""
			for p in c.keys():

				# pprint(c)
				# sys.exit()

				window_hits = [j for j in c.keys() if p <= j <= p+kmer_window]

				# print(window_hits)
				# sys.exit()


				if len(window_hits) > 6:

					count += 1
					print(count, m_name)
					flag = "<---- " + str(len(window_hits))


					print(s_name)
					print("  ",m_name, c.keys(), flag)

		# sys.exit()

	pbar.close()

			# sys.exit()



		# print(len(c))
		# pprint(c.most_common(10))

		# print(len([v for v in c.values() if v > 3]))

		# sys.exit()



	sys.exit('done')




