# generic functions

import sys
# import os
# from pprint import pprint

# import click

from subprocess import PIPE, Popen, call
# from pathlib import Path

from os.path import isfile, isdir

from time import time, sleep
from collections import Counter, deque
from math import floor
# from itertools import count, chain

# from statistics import median, mean

# from Levenshtein import distance

# from timeit import timeit

from math import log10



TOP_READS_HEADER = "cluster\tseq\trank\tdepth\trpm\tlocus_prop"

def top_reads_save(read_c, file, read_equivalent, name):
	cum_count = 0
	top_reads = read_c.most_common(100)
	for rank, read in enumerate(top_reads):

		seq, dep = read
		rpm = round(dep * read_equivalent, 4)

		cum_count += dep

		loc_prop = round(cum_count / sum(read_c.values()), 4)

		with open(file, 'a') as outf:
			print(name, seq, rank, dep, rpm, loc_prop, file=outf, sep='\t')

			if loc_prop >= 0.3:
				break

def parse_locus(locus):
	locus = locus.replace("..", "-").strip()

	chrom = locus.split(":")[0]
	start = int(locus.split(":")[1].split("-")[0])
	stop  = int(locus.split(":")[1].split("-")[1])

	return(chrom, start, stop)

def samtools_faidx(locus, strand, genome_file):

	call = ['samtools', 'faidx', genome_file, locus]

	p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')

	out, err = p.communicate()

	if err != "":
		sys.exit(err)

	out = "".join(out.split("\n")[1:]).strip().upper()

	out = out.replace("T","U")
	# print(out)

	if strand == '-':
		out = out[::-1]
		out = complement(out)
	return(out)


def abundance_to_rgb(abd):

	if abd == 0:
		return((1,1,1))

	log = log10(abd)


	# RED
	if log <= 1:
		r = .4
	elif log < 2:
		r = log - 1
	else:
		r = 1

	# GREEN
	if log <= 1:
		g = .4
	elif log < 2:
		g = 1
	elif log <= 3:
		g = 3 - log
	else:
		g = 0

	# BLUE
	if log <= 1:
		b = 1
	elif log < 2:
		b = 2 - log
	elif log <= 3:
		b = 0
	elif log < 4:
		b = log - 3
	else:
		b = 1

	rgb = (round(c,2) for c in [r,g,b])

	return(rgb)


class percentageClass():
	'''used to make percent timers on known iterations'''

	def __init__(self, increment, total):
		self.increment = increment
		self.total = total

		self.points = [floor(p * increment / 100 * total) for p in range(int(100/increment))]
		self.points.append(total-1)

		# print(self.points)
		# print(total)
		# sys.exit()

	def get_percent(self, i):
		try:
			perc = self.points.index(i) * self.increment
		except ValueError:
			perc = False

		return(perc)


def get_global_depth(output_directory, alignment_file, force=False, aggregate_by=['rg','chrom','length']):
	depth_file = f"./{output_directory}/GlobalDepth.txt"

	header = ['rg','chrom','length','abundance']

	if not isfile(depth_file) or force:
		call = ['samtools', 'view', '-F', '4', alignment_file]

		c = Counter()

		p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')

		print(f"{depth_file} not found.")
		print("Reading alignment to find global depth dimensions...")
		print('  "." = 1M reads')

		rgs     = set()
		chroms  = set()
		lengths = set()

		for i,line in enumerate(p.stdout):
			if (i+1) % 1000000 == 0:
				print(".", end='', flush=True)
				if (i+1) % 10000000 == 0:
					print(" ", end='', flush=True)
					if (i+1) % 100000000 == 0:
						print("\n", end='', flush=True)

			line = line.strip().split("\t")

			rg     = line[18][5:]
			length = int(line[5][:-1])
			chrom  = line[2]

			key = (rg,length,chrom)

			rgs.add(rg)
			lengths.add(length)
			chroms.add(chrom)

			c[key] += 1

		print()
		# print(c.keys())
		with open(depth_file, 'w') as outf:
			print("\t".join(header), file=outf)
			for rg in rgs:
				for chrom in chroms:
					for length in lengths:
						print(rg, chrom, length, c[(rg,length, chrom)], sep='\t', file=outf)

	out_c = Counter()

	indicies = [i for i,h in enumerate(header) if h in aggregate_by]

	multiples = len(indicies) > 1
		

	with open(depth_file, 'r') as f:
		head_line = f.readline()
		for line in f:
			line = line.strip().split('\t')

			if multiples:
				key = tuple([line[i] for i in indicies])
			else:
				key = line[indicies[0]]

			freq = int(line[-1])
			if freq:
				out_c[key] += freq

	return(out_c)







def get_rg_depth(output_directory, alignment_file, force=False):
	depth_file = f"./{output_directory}/RGdepth.txt"

	if isfile(depth_file) and not force: #and not force:
		with open(depth_file, 'r') as f:

			line = f.readline()
			c = Counter()
			
			for line in f:
				line = line.strip().split("\t")


				c[(line[0], int(line[1]))] += int(line[3])



		return(c)


	print()
	print('reading annotation RG depth...')

	c = Counter()
	rg_c = Counter()

	call = ['samtools', 'view', '-F', '4', alignment_file]



	p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')

	for i,line in enumerate(p.stdout):
		line = line.strip().split("\t")

		rg = line[18][5:]
		length = int(line[5][:-1])

		c.update([(rg, length)])
		rg_c.update([rg])




		# if i > 1000000:
		# 	p.kill()
		# 	break

	p.wait()


	with open(depth_file, 'w') as outf:
		print("rg\tlength\tprop\tabundance", file=outf)
		for rg in rg_c.keys():
			for r in range(15,31):
				prop = round(c[(rg,r)] / rg_c[rg], 4)
				# prop_highest = round(c[r] / highest, 4)
				print(rg, r, prop, c[(rg,r)], sep='\t', file=outf)
	return(c)

		# with open(depth_file, 'w') as outf:
		# 	print(alignment_file, depth, sep='\t', file=outf)

def samtools_depth(bam, annotation_readgroups, locus=False):

	if bam.endswith('.bam'):
		index_file = f"{bam}.bai"
	elif bam.endswith('.cram'):
		index_file = f"{bam}.crai"

	if not isfile(index_file):
		# call = f"samtools index {bam}"
		call= ['samtools','index',bam]
		p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')
		out,err = p.communicate()

	# call = f"samtools view -@ 4 -F 4 {bam}"
	call = ['samtools', 'view', '-h','-F', '4']

	for rg in annotation_readgroups:
		call += ['-r', rg]

	call += [bam]


	if locus:
		call.append(locus)

	p1 = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')
	p2 = Popen(['samtools','depth','-a','-'], stdin=p1.stdout, stdout=PIPE)

	for line in p2.stdout:
		line = line.strip().split()

		depth = int(line[-1])

		yield(depth)

	p2.wait()





def samtools_view(bam, dcr_range=False, non_range=False, locus=False, rgs=[]):

	if bam.endswith('.bam'):
		index_file = f"{bam}.bai"
	elif bam.endswith('.cram'):
		index_file = f"{bam}.crai"

	if not isfile(index_file):
		# call = f"samtools index {bam}"
		call= ['samtools','index',bam]
		p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')
		out,err=p.communicate()
		# print(out)
		# print(err)

		# print("WHY AM I INDEXING???")


	# call = f"samtools view -@ 4 -F 4 {bam}"
	call = ['samtools', 'view', '-F', '4']
	
	for rg in rgs:
		call += ['-r', rg]

	call += [bam]


	if locus:
		call.append(locus)
		
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

		seq = seq.replace("T", "U")

		rg = line[18].lstrip("RG:Z:")

		if dcr_range and non_range:
			if length in dcr_range:
				size = 'dcr'
			elif length in non_range:
				size = 'non'
			else:
				size = False
		else:
			size = False



		yield(strand, length, size, sam_pos, sam_chrom, rg, seq, read_id)

	p.wait()

def get_chromosomes(file, output_directory=False):
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


	if output_directory:
		with open(f"./{output_directory}/ChromSizes.txt", 'w') as outf:
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


		# if not message.endswith("\r"):
		self.log.write(self.clear_ansi(message))
			# print(self.clear_ansi(message), end='', file=self.log)

	def flush(self):
		self.terminal.flush()
		self.log.flush()



def inf_counter():
	i = 1
	while True:
		yield(i)
		i += 1



def RNAfold(seq):
	assert len(seq) > 0, f"hairpin is length 0\n{seq}"

	call = ['RNAfold', '--noPS']

	p = Popen(call,
				  stdout=PIPE,
				stderr=PIPE,
				stdin=PIPE,
				encoding='utf-8')
	out, err = p.communicate(f">{time()}\n{seq}")


	mfe = float(out.strip().split()[-1].strip(")").strip("("))
	fold = out.strip().split("\n")[2].split()[0]


	pairing = []
	stack = []

	for i,f in enumerate(fold):
		if f == ".":
			pairing.append('.')

		elif f == "(":
			stack.append(i)
			pairing.append(i)

		elif f == ")":
			pairing.append(stack.pop())


	corr_pairing = []

	for i,p in enumerate(pairing):
		if p == '.':
			cp = p
		else:
			cp = [cor_i for cor_i, cor_p in enumerate(pairing) if cor_p == p and cor_i != i][0]
		corr_pairing.append(cp)

	# pairing = corr_pairing[:]

	return(fold, mfe, corr_pairing)



def complement(s):
	d = {"U" : "A", 
	"A":"U", "G":"C", "C":"G", "N":"N"}

	s = "".join([d[letter] for letter in s])
	return(s)




# def samtools_view(bam, dcr_range=False, non_range=False, locus=False):

# 	if not isfile(f"{bam}.bai"):
# 		# call = f"samtools index {bam}"
# 		call= ['samtools','index',bam]
# 		p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')
# 		out,err=p.communicate()
# 		# print(out)
# 		# print(err)


# 	# call = f"samtools view -@ 4 -F 4 {bam}"
# 	call = ['samtools', 'view', '-F', '4', bam]


# 	if locus:
# 		call.append(locus)
		
# 	# print(call)
# 	p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')

# 	for i,line in enumerate(p.stdout):

# 		line = line.strip().split()

# 		# read_id, flag, sam_chrom, sam_pos, _, length, _, _, _, _,_,_,_,_,_,_,_,_,rg= line

# 		read_id = line[0]
# 		flag = line[1]
# 		seq = line[9]

# 		if flag == "16":
# 			strand = '-'
# 		elif flag == '0':
# 			strand = "+"
# 		else:
# 			strand = False

# 		# print(line)
# 		length = int(line[5].rstrip("M"))
# 		# sam_pos = int(sam_pos)

# 		# length = len(line[9])

# 		sam_pos = int(line[3])
# 		sam_chrom = line[2]

# 		seq = seq.replace("T", "U")

# 		rg = line[18].lstrip("RG:Z:")

# 		if dcr_range and non_range:
# 			if length in dcr_range:
# 				size = 'dcr'
# 			elif length in non_range:
# 				size = 'non'
# 			else:
# 				size = False
# 		else:
# 			size = False



# 		yield(strand, length, size, sam_pos, sam_chrom, rg, seq)

# 	p.wait()

# def get_chromosomes(file,output_directory):
# 	chromosomes = []
# 	rgs = []
# 	# call = f"samtools view -@ 4 -H {file}"

# 	call = ['samtools','view','-H', file]
# 	# print(call)

# 	p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')
# 	out, err = p.communicate()

# 	for o in out.strip().split("\n"):
# 		o = o.strip().split('\t')
# 		if o[0] == "@SQ":
# 			name = o[1].split(":")[-1]
# 			length = int(o[2].split(":")[-1])
# 			chromosomes.append((name,length))
# 		if o[0] == "@RG":
# 			rgs.append(o[1].split(":")[-1])


# 	with open(f"./{output_directory}/ChromSizes.txt", 'w') as outf:
# 		for chrom, size in chromosomes:
# 			print(chrom, size, sep='\t', file=outf)

# 	return(chromosomes, rgs)




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


# class Logger(object):
# 	def __init__(self, file_name):
# 		self.terminal = sys.stdout
# 		self.file_name = file_name
# 		self.log = open(file_name, "w")
# 		# with open(file_name, "w") as outf:
# 		# 	outf.write("")

# 	def clear_ansi(self, message):
# 		return(message.replace("\033[1m", "").replace("\033[0m",""))

# 	def write(self, message):
# 		self.terminal.write(message)
# 		# with open(self.file_name, 'a') as outf:
# 		# 	outf.write(message)  
# 		self.log.write(self.clear_ansi(message))

# 	def flush(self):
# 		self.terminal.flush()
# 		self.log.flush()



def inf_counter():
	i = 1
	while True:
		yield(i)
		i += 1








