# Target predictions

import os
import click
from subprocess import PIPE, Popen, call
from pathlib import Path
from os.path import isfile, isdir
from pprint import pprint
from tqdm import tqdm
from multiprocessing import Lock, Process, Queue, current_process, Pool
import re

from .generics import *
from .cli import cli

from itertools import groupby


def plex_by_file(m_file, s_file):

	z = 30

	call = ['RNAplex', '-f', '2', '-z', str(z), '-t', m_file, '-q', s_file]
	p = Popen(call, stdout=PIPE, stdin=PIPE, encoding='utf-8')


	while True:

		m_name = p.stdout.readline().strip().lstrip(">")
		if m_name == '':
			break
		s_name = p.stdout.readline().strip().lstrip(">")
		out    = p.stdout.readline().strip().split()
		del out[2]

		out = [m_name, s_name] + out

		yield(out)







def plex_by_stdin(q_name, q_seq, t_name, t_seq, verbose=False):

	if verbose:
		print(q_name, t_name)

	z = len(q_seq) + 10

	call = ['RNAplex', '-f', '2', '-z', str(z)]
	p = Popen(call, stdout=PIPE, stdin=PIPE, encoding='utf-8')

	if verbose:
		print()
		print("[plex] target", t_name, t_seq)
		print("[plex] query ", q_name, q_seq)


	inp = f">{t_name}\n{t_seq}\n>{q_name}\n{q_seq}"
	# print(inp)

	print(f'command approximation:\necho "{inp}" | {" ".join(call)}')
	sys.exit()


	out, err = p.communicate(inp)
	out = out.strip().split("\n")
	if len(out) > 4:
		sys.exit()
	if verbose:
		print()
		for o in out:
			print("[plex]", o)
		print()
	out = out[-1].split()
	del out[2]
	
	# if verbose:
	# 	print('[plex]', out)
	# 	print()
	return(out)



def clean_plex(hit, s_full_seq, m_full_seq, verbose=False):

	m_name = hit[0]
	s_name = hit[1]

	if verbose:
		print("###############")
		print(m_name)
		print(s_name)

	s_from, s_to = [int(h)-1 for h in hit[4].split(",")]
	m_from, m_to = [int(h)-1 for h in hit[3].split(",")]


	# m_window_seq = m_full_seq[window_start : window_stop + 1]




	dots = hit[2].split("&")

	if verbose:
		print()
		print("Inputs:")
		print(">mRNA")
		print(m_full_seq)
		print(">sRNA")
		print(s_full_seq)




		print()
		print("Plex output:")
		print(m_full_seq[m_from : m_to+1])
		print(dots[0])
		print(dots[1][::-1])
		# print(s_full_seq[s_from : s_to])
		print(s_full_seq[s_from : s_to + 1][::-1])


	### Filling in gaps

	left_gap = len(s_full_seq) - s_to - 1
	right_gap = s_from + 1

	# s_seq = s_full_seq[s_left-1:s_right]
	s_seq = s_full_seq

	m_i_left = m_from - left_gap 
	if m_i_left < 0:
		return(False,False)

	m_i_right = m_to + right_gap

	m_seq = m_full_seq[m_i_left : m_i_right]



	m_deq = list(m_seq)
	s_deq = list(s_seq[::-1])
	m_dot = list(left_gap * "?" + dots[0] + (right_gap-1) * "?")
	s_dot = list(left_gap * "?" + dots[1][::-1] + (right_gap-1) * "?")

	if len(m_deq) != len(m_dot):
		return(False, False)


	# print(left_gap * "-" + dots[0] + (right_gap-1) * "-")
	# print(left_gap * "-" + dots[1][::-1] + (right_gap-1) * "-")
	# print(s_seq[::-1])


	if verbose:
		print()
		print("Modeled positions:")
		print("".join(m_deq))
		print("".join(m_dot))
		print("".join(s_dot))
		print("".join(s_deq))
		print()


	pair_d = {
		"UA":"|", "AU":"|", 
		"UG":":", "GU":":",
		"GC":"|", "CG":"|",
		"CA":".", "AC":".",
		"CC":".", "GG":".", 
		"UU":".", "AA":".",
		"GA":".", "AG":".",
		"UC":".", "CU":"."
	}

	dot = ''

	i = 0
	while True:
		try:
			m = m_dot[i].replace("(","|")
		except IndexError:
			break

		try:
			s = s_dot[i].replace(")","|")
		except:
			s = " "

		# print(m,s)

		if s == ' ':
			dot += " "
			s_deq.insert(i, " ")

		elif m != s and (m == "|" or s == "|"):
			if m == "|":
				m_deq.insert(i, "-")
				m_dot.insert(i, "-")
			elif s == "|":
				s_deq.insert(i, "-")
				s_dot.insert(i, "-")
			else:
				sys.exit("huh?")

			dot += "-"

		else:
			dot += pair_d[m_deq[i]+s_deq[i]]

		# print(m_dot[i].replace("(","|"),s_dot[i].replace(")","|"))
		i += 1

	m_deq = "".join(m_deq)
	s_deq = "".join(s_deq)


	dot = re.sub("^[\\-]+", " ", dot)
	s_deq = re.sub("^[\\-]+", " ", s_deq)
		
	if verbose:
		print()
		print("Clean duplex:")
		print(m_deq)
		# print("".join(m_dot))
		print(dot)
		# print("".join(s_dot))
		print(s_deq)
	# sys.exit()


	c = Counter()
	allenscore = 0

	pair_types = {
		"." : 'unpaired',
		"|" : 'paired-watcri',
		":" : 'paired-gu',
		"-" : 'bulge',
		" " : 'absent'
	}

	asc_d = {
		"|" : 0,
		"." : 1,
		"-" : 1,
		":" : 0.5,
		" " : 0
	}

	yasma_d = {
		"|" : 0,
		"." : 1,
		"-" : 1.5,
		":" : 0.5,
		" " : 0
	}

	yasma_score =0

	for i,d in enumerate(dot[::-1]):
		c[pair_types[d]] += 1

		c['YasmaScore'] += yasma_d[d]


		if i == 0:
			asc = 0

		else:
			asc = asc_d[d]

			if i <= 12:
				asc += asc

		# print(i, d, asc)
		c['AllenScore']  += asc


	c['InteractionLength'] = len(dot.strip())



	# for p, name in pair_types:
	# 	type_d[name] = dot.count(p)




	if verbose:
		# print([(k, sum(1 for i in g)) for k,g in groupby(list(dot))])
		print()
		print("Scoring:")
		pprint(c)
		# sys.exit()



	duplex = (m_deq, dot, s_deq)

	# pprint(duplex)
	return(duplex, c)





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


	s_dict = {}

	s_file_name = "queries.fa"
	with open(s_file_name, 'w') as outf:
		sRNAs = []
		with open(read_file, 'r') as f:
			header = f.readline()
			for line in f:
				line = line.strip().split("\t")

				cluster, seq, rank, _, _, _ = line
				rank = int(rank)

				if rank <= max_rank:
					sRNAs.append((cluster, seq))

					print(">" + cluster, file=outf)
					print(seq, file=outf)

					s_dict[cluster] = seq


			# sys.exit()


	kmer_size = 6
	# kmer_window = 30
	kmer_window_expansion = 30
	kmer_threshold = 6

	kmer_d = {}

	window_d = {}

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
				# mRNA_lengths[header] = len(seq)

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



	output_file = "target_output.txt"
	with open(output_file, 'w') as outf:

		header = ['m_name', 's_name', 'm_seq', 'interaction', 's_seq']
		score_keys = ['InteractionLength', 'AllenScore', 'YasmaScore', 'paired-watcri', 'paired-gu', 'unpaired', 'bulge']

		header += score_keys

		print("\t".join(header), file=outf)


	# print(len(mRNAs))

	pbar = tqdm(total=len(mRNAs))
	for m_i, mRNA in enumerate(mRNAs):
		kmer_d = {}

		m_name, m_seq = mRNA
		# print(m_name)#, m_seq)
		pbar.update()

		m_file_name = f"{m_name}_target.fa"
		with open(m_file_name, 'w') as outf:
			print(">" + m_name, file=outf)
			print(m_seq, file=outf)






		# complement_mRNA = m_seq[::-1]

		# for k_i, kmer in enumerate(get_kmers(complement_mRNA, k=kmer_size)):

		# 	try:
		# 		kmer_d[kmer]
		# 	except KeyError:
		# 		kmer_d[kmer] = Counter()

		# 	kmer_d[kmer][k_i] += 1


		# for s_name, s_seq in sRNAs:

		# print(s_name)
		verbose = False

		for hit in plex_by_file(m_file_name, s_file_name):

			m_name, s_name = hit[:2]

		# hit = plex_by_stdin(q_name=s_name, q_seq=s_seq, t_name=m_name, t_seq=m_seq, verbose=verbose)
			duplex, scoring_c = clean_plex(hit, s_dict[s_name], m_seq, verbose=verbose)

			if duplex:
				m, dot, s = duplex
				out = [m_name, s_name, m, dot, s]

				for key in score_keys:
					out.append(scoring_c[key])

				with open(output_file, 'a') as outf:
					print("\t".join(map(str, out)), file=outf)

		os.remove(m_file_name)
		# input()

		# for d in duplex:
		# 	print(d)
		# pprint(scoring_c)

		# if scoring_c and scoring_c['AllenScore'] < 6:
		# 	sys.exit()

		# if m_i > 10:
		# 	break

	pbar.close()






