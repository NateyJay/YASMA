# generic functions

import sys
import click
from os import stat
import os
from git import Repo
from pprint import pprint

# import click

from subprocess import PIPE, Popen, call
# from pathlib import Path

from pathlib import Path, PurePath
from os.path import isfile, isdir

from time import time, sleep
from collections import Counter, deque
from math import floor
# from itertools import count, chain

from statistics import median, mean

# from Levenshtein import distance

# from timeit import timeit
import re

from math import log10

import pyBigWig

import json
from random import sample, seed, shuffle
import pysam


# import git




ENCODING='cp850'
# ENCODING=ENCODING

# def make_bash_safe_path(path):

# 	print(path)
# 	print(os.path.relpath(path))
# 	sys.exit()

# 	string = str(path)
# 	string = string.replace(" ", "\ ")
# 	return(string)



class requirementClass():
	def __init__(self):
		self.reqs = []


	def add_cutadapt(self):
		try:
			p = Popen(['cutadapt', '--version'], stdout=PIPE, stderr=PIPE, encoding=ENCODING)
			out,err = p.communicate()
			found=True
			version = out.strip()

		except FileNotFoundError:
			found=False
			version=''

		self.reqs.append(('cutadapt', found, version))

	def add_samtools(self):
		try:
			p = Popen(['samtools'], stdout=PIPE, stderr=PIPE, encoding=ENCODING)
			out,err = p.communicate()
			err = err.split("\n")
			# print(err)
			found=True
			version = err[2].split()[1]

		except FileNotFoundError:
			found=False
			version=''

		self.reqs.append(('samtools', found, version))

	def add_shortstack(self):
		try:
			p = Popen(['ShortStack', '-v'], stdout=PIPE, stderr=PIPE, encoding=ENCODING)
			out,err = p.communicate()
			out = out.split("\n")
			# print(out)
			found=True
			version = out[0].split()[-1]

		except FileNotFoundError:
			found=False
			version=''

		self.reqs.append(('ShortStack', found, version))


	def add_bowtie(self):
		try:
			p = Popen(['bowtie', '--version'], stdout=PIPE, stderr=PIPE, encoding=ENCODING)
			out,err = p.communicate()
			out = out.split("\n")
			found=True
			version = out[0].split()[-1]

		except FileNotFoundError:
			found=False
			version=''

		self.reqs.append(('bowtie', found, version))

	def add_rnafold(self):
		try:
			p = Popen(['RNAfold', '--version'], stdout=PIPE, stderr=PIPE, encoding=ENCODING)
			out,err = p.communicate()
			out = out.split("\n")
			found=True
			version = out[0].split()[-1]

		except FileNotFoundError:
			found=False
			version=''

		self.reqs.append(('RNAfold', found, version))

	def add_RNAfold(self):

		try:
			p = Popen(['RNAfold', '--version'], stdout=PIPE, stderr=PIPE, encoding=ENCODING)
			out,err = p.communicate()
			out = out.split("\n")
			found=True
			version = out[0].split()[-1]

		except FileNotFoundError:
			found=False
			version=''

		self.reqs.append(('RNAfold', found, version))

	def add_bedtools(self):

		try:
			p = Popen(['bedtools', '--version'], stdout=PIPE, stderr=PIPE, encoding=ENCODING)
			out,err = p.communicate()
			out = out.split("\n")
			found=True
			version = out[0].split()[-1]

		except FileNotFoundError:
			found=False
			version=''

		self.reqs.append(('bedtools', found, version))

	def check(self):

		fail=False

		print(color.BOLD + "Requirements:" + color.END)

		for tool, found, version in self.reqs:
			if not found:
				fail = True

			if found:
				print(color.BOLD + '[x]' + color.END, tool, "->", version)
			else:
				print(color.BOLD + '[ ]' + color.END, tool)

		
		if fail:
			sys.exit("Error: requirements not met!")

		print()

class inputClass():

	def __init__(self, params):

		try:
			if params['override']:
				self.override = True
			else:
				self.override = False
		except KeyError:
			self.override = False

		output_directory = params['output_directory']
		output_directory = output_directory.rstrip("/")

		self.output_directory = Path(output_directory)
		

		if self.output_directory == Path("."):
			self.output_directory = Path().cwd()

		else:
			self.output_directory.mkdir(parents=True, exist_ok=True)


		project_name = self.output_directory.name

		self.file = Path(output_directory, "inputs.json")


		self.inputs = {'project_name' : output_directory}


		self.input_list = [
			"untrimmed_libraries",
			"trimmed_libraries",
			"adapter",
			"alignment_file",
			'annotation_readgroups',
			'genome_file',
			'jbrowse_directory',
			'gene_annotation_file'
			]


		self.paths = [
			"untrimmed_libraries",
			"trimmed_libraries",
			"alignment_file",
			'genome_file',
			'jbrowse_directory',
			'gene_annotation_file'
			]

		for i in self.input_list:
			self.inputs[i] = None


		if isfile(self.file):
			try:
				self.read()
			except json.decoder.JSONDecodeError:
				print("DECODER ERROR!")
				pass

		self.parse(params)

		print("params:")
		pprint(params)
		print("inputs:")
		pprint(self.inputs)
		sys.exit()

		self.check_paired_end()



		self.write()
	


	def read(self):
		with open(self.file, 'r') as f:
			self.inputs = json.load(f)

		self.decode_inputs()


	def write(self):

		self.encode_inputs()
		with open(self.file, 'w') as outf:
			outf.write(json.dumps(self.inputs, indent=2))
		self.decode_inputs()



	def encode_inputs(self):

		od = self.output_directory
		print("OD:", od)

		def encode_path(path):


			path = str(path.absolute())

			if "/"+str(od)+"/" in path:
				path = "<OD>/" + path.split(str(od)+"/")[-1]

			return(path)


		for p in ["untrimmed_libraries", "trimmed_libraries"]:
			if self.inputs[p]:
				for i in range(len(self.inputs[p])):
					self.inputs[p][i] = encode_path(self.inputs[p][i])

		for p in ["alignment_file", 'genome_file', 'jbrowse_directory', 'gene_annotation_file']:
			if self.inputs[p]:
				self.inputs[p] = encode_path(self.inputs[p])






	def decode_inputs(self):
		od = self.output_directory


		def decode_path(path):
			# print(path)
			if "<OD>/" in path:
				path = path.split("<OD>/")[-1]
				path = Path(od, path)
				# print(path)
			elif "<OUTPUT_DIRECTORY>/" in path:
				path = path.split("<OUTPUT_DIRECTORY>/")[-1]
				path = Path(od, path)
				# print(path)
			else:
				path = Path(path)


			# path = Path(os.path.relpath(path))
			return(path)

		for p in ["untrimmed_libraries", "trimmed_libraries"]:
			if self.inputs[p]:
				for i in range(len(self.inputs[p])):
					self.inputs[p][i] = decode_path(self.inputs[p][i])


		for p in ["alignment_file", 'genome_file', 'jbrowse_directory', 'gene_annotation_file']:
			if self.inputs[p]:

				self.inputs[p] = decode_path(self.inputs[p])




	def parse(self, params):
		# self.params = params




		for p in ['annotation_readgroups', 'trimmed_libraries','untrimmed_libraries']:
			if p in params:
				try:
					params[p] = list(params[p])
				except TypeError:
					pass

		for option in self.input_list:

			try:
				value = params[option]
			except KeyError:
				value = None

			if option in self.paths and value:


				if isinstance(value, list):
					for i in range(len(value)):
						value[i] = Path(value[i])

				else:
					value = Path(value)

			self.add(option, value)



	def add(self, option, value):

		try:
			saved_value = self.inputs[option]
		except KeyError:
			saved_value = ''

		if not saved_value:
			self.inputs[option] = value

		elif not value:
			pass

		# elif saved_value != value and self.override:
		else:
			print(f"  Override!")
			print(f"  Replace: ... '{self.inputs[option]}'")
			print(f"  with: ...... '{value}'")
			self.inputs[option] = value

		# elif saved_value != value:
		# 	print(f"  Warning: input for option '{color.BOLD}{option}{color.END}' does not match logged value")

		# 	print(f"  Replace: ... '{self.inputs[option]}'")
		# 	print(f"  with: ...... '{value}'")

		# 	res = input("   (y)es or (n)o?\n")

		# 	if res == 'y':
		# 		self.inputs[option] = value



	def get(self):
		out = []

		for i in self.input_list:
			out.append(self.inputs[i])

		return(out)

	def check(self, required_options):

		if 'output_directory' not in required_options:
			required_options = ['project_name'] + required_options

		self.required_options = required_options

		print(color.BOLD + "Required options:" + color.END)
		pass_check = True
		offset = 25

		for option in required_options:
			value = self.inputs[option]

			print(color.BOLD + '[', end='')
			if not value:
				pass_check=False
				print(" ", end='')

			else:
				print(f"x", end='')

			print('] ' + color.END, end='')

			print(f"{option}:", "." * (offset-len(option)), f"{color.BOLD}{value}{color.END}")


			# 	pass_check = False
			# 	print(f'Error: required option {option} not supplied or in logged inputs...')
		print()
		if not pass_check:


			sys.exit("Error: one or more essential options are not provided in 'inputs.json' or program call")

		print(color.BOLD + "Other options:" + color.END)

		options = list(self.inputs.keys())
		options += [o for o in self.inputs.keys() if o not in options]
		options = [o for o in options if o not in required_options]

		for option in options:

			if option in self.inputs:
				value = self.inputs[option]
			else:
				value = self.params[option]

			print(f"    {option}:", "." * (offset-len(option)), f"{color.BOLD}{value}{color.END}")

		print()



	def check_chromosomes(self):

		genome_chromosomes = set()
		if self.inputs['genome_file']:
			genome_file = self.inputs['genome_file']
			genome_index_file = genome_file.parent / (genome_file.name + ".fai")

			if not isfile(genome_file):
				sys.exit(f"Error: genome_file not found {genome_file}")


			if not isfile(genome_index_file):
				call = ['samtools', 'faidx', str(genome_file)]
				print(" ".join(call))
				p = Popen(call, stdout=PIPE, stderr=PIPE)
				p.wait()

			with open(genome_index_file, 'r') as f:

				for line in f:
					genome_chromosomes.add(line.split()[0])

			# print(genome_chromosomes)


		gene_annotation_chromosomes = set()
		if self.inputs['gene_annotation_file']:

			if not isfile(self.inputs['gene_annotation_file']):
				sys.exit(f"Error: gene_annotation_file not found {self.inputs['gene_annotation_file']}")

			with open(self.inputs['gene_annotation_file']) as f:

				for line in f:
					if not line.startswith("#"):
						gene_annotation_chromosomes.add(line.split()[0])

			gene_annotation_chromosomes = list(gene_annotation_chromosomes)
			# print(gene_annotation_chromosomes)


		alignment_chromosomes = set()
		if self.inputs['alignment_file']:

			if not isfile(self.inputs['alignment_file']):
				sys.exit(f"Error: alignment_file not found {self.inputs['alignment_file']}")

			p = Popen(['samtools', 'view', '-H', self.inputs['alignment_file']],
				stdout=PIPE, stderr=PIPE, encoding=ENCODING)

			for line in p.stdout:
				if line.startswith("@SQ"):
					alignment_chromosomes.add(line.split("\t")[1][3:])

			p.wait()

			# print(alignment_chromosomes)
		all_chroms = set()
		all_chroms.update(genome_chromosomes)
		all_chroms.update(gene_annotation_chromosomes)
		all_chroms.update(alignment_chromosomes)

		def is_found(chrom, chrom_set):
			if chrom in chrom_set:
				return(" x")
			elif len(chrom_set) == 0:
				return(None)
			else:
				return("  ")

		print("Checking chromosome/scaffold/contig overlaps between inputs...")
		print()
		print('genome', 'align', 'gene', 'chrom/scaffold', sep='\t')
		print("--------------------------------------")
		for chrom in all_chroms:
			print(is_found(chrom, genome_chromosomes), is_found(chrom, alignment_chromosomes), is_found(chrom, gene_annotation_chromosomes), chrom, sep='\t')


	def check_paired_end(self):
		print("check_paired_end input:")
		pprint(self.inputs['untrimmed_libraries'])

		if self.inputs['untrimmed_libraries']:
			libs = self.inputs['untrimmed_libraries']

			first_pairs = []

			possible_basenames = set()

			for lib in libs:
				input_lib = lib

				while lib.suffix:
					lib = lib.with_suffix('')


				if lib.stem.endswith(("_1","_2","_3")):

					possible_basenames.add(lib.stem[:-2])

				if lib.stem.endswith("_1"):
					first_pairs.append(input_lib)


			if len(libs) != len(possible_basenames):
				# print("Likely paired ends!!!")
				# print("Use these libraries for the first pairs")

				self.inputs['untrimmed_libraries'] = first_pairs

		pprint(self.inputs)
			
def validate_glob_path(ctx, param, value):

	paths = list(value)
	paths = [p.strip() for p in paths]
	paths = " ".join(paths)
	paths = paths.strip().split(" ")

	if paths == ['']:
		# print(param)
		# raise click.UsageError("Error: Missing or empty option '-l'")
		return(None)

	full_paths = []
	for path in paths:
		full_paths.append(str(Path(path).absolute()))

		if not isfile(path) and not isdir(path):
			raise click.BadParameter(f"path not found: {path}")

	full_paths = tuple(full_paths)
	return(full_paths)		


def validate_path(ctx, param, value):


	if not value:
		return(None)

	path = value.strip()

	if not isfile(path) and not isdir(path):
		raise click.BadParameter(f"path not found: {path}")

	full_path = str(Path(path).absolute())

	return(full_path)		


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
	# print(" ".join(call))
	p = Popen(call, stdout=PIPE, stderr=PIPE, encoding=ENCODING)

	out, err = p.communicate()

	if err != "":
		print(f"WARNING: {err}")
		# sys.exit(err)

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

		self.removable = list(self.points)

		self.running = 0

		# print(self.points)
		# print(total)
		# sys.exit()

	def check(self, value):
		try:
			if value >= self.removable[0]:
				perc = self.points.index(self.removable.pop(0)) * self.increment
			else:
				perc = None
			return(perc)
		except IndexError:
			return(None)

	def get_percent(self, i):
		try:
			perc = self.points.index(i) * self.increment
		except ValueError:
			perc = False

		return(perc)

	def update(self, i=1):
		self.running += i

		return(self.get_percent(self.running))



# def get_window_depths(alignment_file, chromosome, window_length):
# 	i = 0
# 	window = deque()
# 	d_out = []

# 	sam_iter = samtools_view(alignment_file, locus=chromosome)

# 	read = next(sam_iter)

# 	while True:

# 		_, length, _, pos, chromosome, _, _, _ = read

# 		if pos > i:
# 			i += 1

# 			try:
# 				d = window.popleft()
# 			except IndexError:
# 				d = 0

# 			d_out.append(d)


# 		elif pos == i:

# 			for r in range(length + window_length):

# 				try:
# 					window[r] += 1
# 				except IndexError:
# 					window.append(1)

# 			try:
# 				read = next(sam_iter)
# 			except StopIteration:
# 				return(d_out)
# 			# if i > 1000000:
# 			# 	return(d_out)
# 		elif pos < i:
# 			print(pos, i)
# 			sys.exit("this shouldn't happen")


# def get_lambda(alignment_file, chromosomes, window_length, output_directory):

# 	lambda_file = f"./{output_directory}/Lambdas.txt"
# 	lambda_d = {}


# 	if isfile(lambda_file):
# 		with open(lambda_file, 'r') as f:
# 			for line in f:
# 				c,m = line.strip().split('\t')
# 				m = float(m)
# 				lambda_d[c] = m

# 	else:
# 		window_d = {}


# 		for c, l in chromosomes:

# 			window_d[c] = get_window_depths(alignment_file, c, window_length)


# 			with open(f"./{output_directory}/{c}.dist.txt", 'w') as outf:
# 				for d in window_d[c]:
# 					print(d, file=outf)


# 		print('\n')
# 		print('chrom','median','mean', sep='\t')
# 		for c,l in chromosomes:

# 			window_d[c] = sample(window_d[c],10000)

# 			print(c, median(window_d[c]), round(mean(window_d[c]),4), sep='\t')
			
# 			lambda_d[c] = mean(window_d[c])


# 		with open(lambda_file, 'w') as outf:
# 			for c,m in lambda_d.items():
# 				print(c,m, sep='\t', file=outf)


# 	return(lambda_d)


	# for i,line in enumerate(p.stdout):

	# 	if i % 1000000 == 0:
	# 		print(".", end='', flush=True)
	# 	if i % 10000000 == 0 and i > 0:
	# 		print(" ", end='', flush=True)

	# 	line = line.strip().split("\t")
	# 	print(line)

	# 	flag = line[1]
	# 	chrom = line[2]
	# 	pos = int(line[3])
	# 	length = int(line[5].rstrip("M"))


def read_loci(params):


	results_file = f"{params['output_directory']}/tradeoff/loci.txt"
	with open(results_file, 'r') as f:
		header = f.readline().strip().split("\t")
		header = [h.lower() for h in header]

		for line in f:
			line = line.strip().split("\t")
			d = dict(zip(header, line))

			yield d


def get_global_depth(alignment_file, force=False, aggregate_by=['rg','chrom','length']):

	depth_file = alignment_file.with_suffix(".depth.txt")

	header = ['rg','chrom','length','abundance']


	if not isfile(depth_file) or stat(depth_file).st_size < 50 or force:
		# call = ['samtools', 'view', '-F', '4', str(alignment_file)]
		# print(" ".join(call))
		c = Counter()

		# p = Popen(call, stdout=PIPE, encoding=ENCODING)

		print(f"  {depth_file} not found.")
		print("Reading alignment to find global depth dimensions...")
		print('  "." = 1M reads')

		rgs     = set()
		chroms  = set()
		lengths = set()



		for i, sam_out in enumerate(samtools_view(alignment_file)):
		# for i,line in enumerate(p.stdout):
			if (i+1) % 1000000 == 0:
				print(".", end='', flush=True)
				if (i+1) % 10000000 == 0:
					print(" ", end='', flush=True)
					if (i+1) % 100000000 == 0:
						print("\n", end='', flush=True)

			# line = line.strip().split("\t")

			# rg     = line[18][5:]
			# length = int(line[5][:-1])
			# chrom  = line[2]
			_, length, _, _, chrom, rg, _, _ = sam_out

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
			# print(line)

			if multiples:
				key = tuple([line[i] for i in indicies])
			else:
				key = line[indicies[0]]
			freq = int(line[-1])
			# print(key, freq)

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



	p = Popen(call, stdout=PIPE, stderr=PIPE, encoding=ENCODING)

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
		p = Popen(call, stdout=PIPE, stderr=PIPE, encoding=ENCODING)
		out,err = p.communicate()

	# call = f"samtools view -@ 4 -F 4 {bam}"
	call = ['samtools', 'view', '-h','-F', '4']

	for rg in annotation_readgroups:
		call += ['-r', rg]

	call += [bam]


	if locus:
		call.append(locus)

	p1 = Popen(call, stdout=PIPE, stderr=PIPE, encoding=ENCODING)
	p2 = Popen(['samtools','depth','-a','-'], stdin=p1.stdout, stdout=PIPE)

	for line in p2.stdout:
		line = line.strip().split()

		depth = int(line[-1])

		yield(depth)

	p2.wait()



def counters_to_bigwig(counter_d, chromosomes, file, verbose=True):
	if verbose:
		print(file)
	bw = pyBigWig.open(str(file), 'w')
	bw.addHeader(chromosomes)

	for chrom, chrom_length in chromosomes:

		positions = list(counter_d[chrom].keys())

		positions = [p for p in positions if 0 < p <= chrom_length]

		last_p = 0
		last_v = round(float(counter_d[chrom][0]),4)

		starts = [last_p]
		ends   = []
		values = [last_v]

		for p in positions:
			v = round(float(counter_d[chrom][p]), 4)
			# print(p, v)

			if p == 0:
				values[0] = v

			elif p != last_p + 1 and last_p > 0:
				## for sure a zero gap
				ends.append(last_p)

				starts.append(last_p)
				ends.append(p)
				values.append(0.0)

				starts.append(p)
				values.append(v)


			else:
				if v != last_v:
					ends.append(p)

					starts.append(p)
					values.append(v)

			last_p = p
			last_v = v



		ends.append(chrom_length)
		chroms = [chrom]*len(starts)

		# chroms.pop(0)
		# starts.pop(0)
		# ends.pop(0)
		# values.pop(0)


		# chroms.pop()
		# starts.pop()
		# ends.pop()
		# values.pop()


		# for i in range(len(starts)):
		# 	print(i, chroms[i], starts[i], ends[i], values[i], sep='\t')


		# print('chroms', len(chroms), sep='\t')
		# print('starts', len(starts), sep='\t')
		# print('ends', len(ends), sep='\t')
		# print('values', len(values), sep='\t')

		bw.addEntries(chroms, starts, ends=ends, values=values, validate=False)
		# bw.addEntries(chroms[:100], starts[:100], ends=ends[:100], values=values[:100], validate=False)
		# bw.addEntries(['NC_003070.9','NC_003070.9'],
		# 	[1,2], [2,4], [2.0,3.3])
				
		# bw.addEntries(["chr1", "chr1", "chr1"], [100, 0, 125], ends=[120, 5, 126], values=[0.0, 1.0, 200.0], validate=False)
	bw.close()



def samtools_view(bam, rgs='all', contig=None, start=None, stop=None, threads=4, boundary_rule='loose'): #, read_minmax=(15,30)):




	bamf = pysam.AlignmentFile(bam,'rb', threads=threads)

	if not bamf.has_index():
		print(f'   index not found for {bam}. Indexing with samtools.')

		pysam.index(str(bam))
		bamf.close()
		bamf = pysam.AlignmentFile(bam,'rb', threads=threads)



	if rgs == 'all':
		rgs = [rgd['ID'] for rgd in bamf.header['RG']]

	rgs = set(rgs)



	if boundary_rule == 'tight' and start and stop:

		for read in bamf.fetch(contig=contig, start=start, stop=stop):
			read_length = read.infer_read_length()

			if read.get_tag("RG") in rgs and read.is_mapped and read.overlap == read_length:
				strand = "-" if read.is_reverse else "+"

				# if read_length >= read_minmax[0] and read_length <= read_minmax[-1]:

				yield(strand, 
					read_length, 
					False, 
					read.reference_start, 
					read.reference_name, 
					read.get_tag("RG"), 
					read.get_forward_sequence().replace("T","U"), 
					read.query_name)
	else:

		for read in bamf.fetch(contig=contig, start=start, stop=stop):

			read_length = read.infer_read_length()

			if read.get_tag("RG") in rgs and read.is_mapped:
				strand = "-" if read.is_reverse else "+"

				# if read_length >= read_minmax[0] and read_length <= read_minmax[-1]:
				yield(strand, 
					read_length, 
					False, 
					read.reference_start, 
					read.reference_name, 
					read.get_tag("RG"), 
					read.get_forward_sequence().replace("T","U"), 
					read.query_name)

	bamf.close()

# def samtools_view(bam, dcr_range=False, non_range=False, locus=False, rgs=[], boundary_rule='loose', subsample=None, reverse=False, reverse_chunk = 100):




# 	if boundary_rule == 'tight':
# 		lbound = int(locus.split(":")[-1].split("-")[0])
# 		rbound = int(locus.split(":")[-1].split("-")[1])+1
# 	else:
# 		lbound = False


# 	# call = f"samtools view -@ 4 -F 4 {bam}"
# 	call = ['samtools', 'view', '-F', '4']
	
# 	for rg in rgs:
# 		call += ['-r', rg]

# 	if subsample:
# 		call += ['--subsample', str(subsample)]

# 	call += [str(bam)]

# 	size=False




# 	# print(" ".join(map(str,call)))
		
# 	# print(call)

# 	def parse_line(line):
# 		read_id = line[0]
# 		flag = line[1]
# 		seq = line[9]

# 		if flag == "16":
# 			strand = '-'
# 		elif flag == '0':
# 			strand = "+"
# 		else:
# 			strand = False

# 		length = int(line[5].rstrip("M"))

# 		sam_pos = int(line[3])
# 		sam_chrom = line[2]

# 		seq = seq.replace("T", "U")

# 		rg = line[18].lstrip("RG:Z:")
# 		size=False



# 		return(strand, length, size, sam_pos, sam_chrom, rg, seq, read_id)

# 	if not reverse:
# 		if locus:
# 			call.append(locus)

# 		p = Popen(call, stdout=PIPE, stderr=PIPE, encoding=ENCODING)

# 		for i,line in enumerate(p.stdout):

# 			line = line.strip().split()

# 			# read_id, flag, sam_chrom, sam_pos, _, length, _, _, _, _,_,_,_,_,_,_,_,_,rg= line
# 			strand, length, size, sam_pos, sam_chrom, rg, seq, read_id = parse_line(line)

# 			if lbound and boundary_rule == 'tight':
# 				if sam_pos >= lbound and sam_pos + length + 1 <= rbound:
# 					yield(strand, length, size, sam_pos, sam_chrom, rg, seq, read_id)

# 			else:
# 				yield(strand, length, size, sam_pos, sam_chrom, rg, seq, read_id)

# 		p.wait()

# 	else:
# 		chrom = locus.split(":")[0]

# 		coords = locus.split(":")[1].split("-")

# 		if len(coords) == 2:
# 			left, right = map(int, coords)

# 		elif len(coords) == 1:
# 			left = False
# 			right = int(coords[0])
# 		# start, stop = map(int, locus.split(":")[1].split("-")[1]



# 		read_set=set()
# 		chunk_r = right
# 		while chunk_r > left:

# 			chunk = []
# 			chunk_l = chunk_r - reverse_chunk

# 			locus = f"{chrom}:{chunk_l}-{chunk_r}"

# 			sub_call = call.copy()
# 			sub_call.append(locus)

# 			p = Popen(sub_call, stdout=PIPE, stderr=PIPE, encoding=ENCODING)

# 			for i,line in enumerate(p.stdout):

# 				line = line.strip().split()
# 				if line[0] not in read_set:
# 					chunk.append(line)
# 					read_set.add(line[0])


# 			p.wait()

# 			for line in chunk[::-1]:
# 				strand, length, size, sam_pos, sam_chrom, rg, seq, read_id = parse_line(line)

# 				yield(strand, length, size, sam_pos, sam_chrom, rg, seq, read_id)

# 			chunk_r -= reverse_chunk







def get_chromosomes(file):
	chromosomes = []
	rgs = []
	# call = f"samtools view -@ 4 -H {file}"

	call = ['samtools','view','-H', file]
	# print(" ".join(call))

	p = Popen(call, stdout=PIPE, stderr=PIPE, encoding=ENCODING)
	out, err = p.communicate()


	for o in out.strip().split("\n"):
		o = o.strip().split('\t')
		if o[0] == "@SQ":
			name = o[1].split(":")[-1]
			length = int(o[2].split(":")[-1])
			chromosomes.append((name,length))
		if o[0] == "@RG":
			rgs.append(o[1].split(":")[-1])


	# if output_directory:
	# 	with open(f"./{output_directory}/chrom_sizes.txt", 'w') as outf:
	# 		for chrom, size in chromosomes:
	# 			print(chrom, size, sep='\t', file=outf)

	return(chromosomes, rgs)

# def get_chromosome_depths(bam_file):
# 	print("reading chromosome depths from idxstats...")
# 	print()

# 	p = Popen(['samtools', 'idxstats', bam_file], stdout=PIPE, stderr=PIPE, encoding=ENCODING)

# 	chrom_depths = Counter()
# 	for line in p.stdout:
# 		line = line.strip().split("\t")
# 		chrom_depths[line[0]] += int(line[2])
# 		# print(line)

# 	p.wait()

# 	return(chrom_depths)



# def check_rgs(annotation_readgroups, bam_rgs):

# 	if annotation_readgroups[0].lower() == 'all':
# 		annotation_readgroups = bam_rgs
# 	else:
# 		for rg in annotation_readgroups:
# 			if rg not in bam_rgs:
# 				sys.exit(f"Error: provided readgroup '{rg}' not found within bamfile header:\n{bam_rgs}")

# 	annotation_readgroups = set(annotation_readgroups)
# 	return(annotation_readgroups)
	# def get(self):


class color:
	PURPLE = '\033[95m'
	CYAN = '\033[96m'
	DARKCYAN = '\033[36m'
	BLUE = '\033[94m'
	GREEN = '\033[92m'
	YELLOW = '\033[93m'
	RED = '\033[91m'
	BOLD = '\033[1m'
	UNDERLINE = '\033[4m'
	END = '\033[0m'

class color:
	PURPLE = ''
	CYAN = ''
	DARKCYAN = ''
	BLUE = ''
	GREEN = ''
	YELLOW = ''
	RED = ''
	BOLD = ''
	UNDERLINE = ''
	END = ''


class bigwigClass():
	'''A class to handle producing rpm bigwig files from a counter object c[pos] = depth'''

	def __init__(self, file, rpm_threshold, total_reads, chromosomes, strand= "+"):
		self.file = file
		# self.file = f"./{output_directory}/Coverages/{file}.wig"
		self.bw = pyBigWig.open(self.file, 'w')

		self.rpm_threshold = rpm_threshold
		self.total_reads   = total_reads

		self.bw.addHeader(chromosomes)

		if strand == "+":
			self.strand = 1
		elif strand == "-":
			self.strand = -1 


	def add_chrom(self, c, chrom, chrom_length, peaks_only=False):
	
		depths = [0.0000]

		for p in range(chrom_length):

			depth = round(c[p] / self.total_reads * 1000000, 4) * self.strand

			if peaks_only and depth < self.rpm_threshold:
				depth = 0

			if depth != depths[-1]:
				# chr19 49302000 49302300 -1.0
				# print(f"variableStep  chrom={chrom}  span={len(depths)}")
				# print(f"{p-len(depths)+2} {depths[-1]}")

				# self.bw.addEntries(chrom, [p-len(depths)+2], values=[depths[-1]], span=len(depths))
				# print(chrom, p-len(depths)+2, p+2, depths[-1], sep='\t')


				# self.bw.addEntries([chrom], [p-len(depths)+2], [p+2], values=[depths[-1]])
				self.bw.addEntries(chrom, [p-len(depths)+2], values=[depths[-1]], span=len(depths))

				depths = []

			depths.append(depth)


class wiggleClass():
	'''A class to handle producing rpm wig/bigwig files from a counter object c[pos] = depth'''

	def __init__(self, file, rpm_threshold, total_reads):
		self.file = file
		# self.file = f"./{output_directory}/Coverages/{file}.wig"
		self.outf = open(self.file, 'w')

		self.rpm_threshold = rpm_threshold
		self.total_reads   = total_reads



	def add_chrom(self, c, chrom, chrom_length, peaks_only=False):
	
		depths = [0.0000]

		for p in range(chrom_length):

			depth = round(c[p] / self.total_reads * 1000000, 4)

			if peaks_only and depth < self.rpm_threshold:
				depth = 0

			if depth != depths[-1]:
				# chr19 49302000 49302300 -1.0
				print(f"variableStep  chrom={chrom}  span={len(depths)}", file=self.outf)
				print(f"{p-len(depths)+2} {depths[-1]}", file=self.outf)

				depths = []

			depths.append(depth)

	# def add(self, val, pos, chrom):

	# 	if val != self.val:
	# 		span = pos - self.start_pos

	# 		if span > 0:

	# 			print(f"variableStep chrom={chrom} span={span}", file=self.outf)
	# 			print(f"{self.start_pos} {self.val}", file=self.outf)

	# 			self.val = val
	# 			self.start_pos = pos

	def convert(self, output_directory, cleanup=False):

		self.outf.close()

		wig = self.file

		bigwig = wig.replace(".wig", ".bigwig")

		print(f"  {wig} -> {bigwig}", flush=True)

		call = f"wigToBigWig {wig} ./{output_directory}/ChromSizes.txt {bigwig}"

		p = Popen(call.split(), stdout=PIPE, stderr=PIPE, encoding=ENCODING)

		out, err= p.communicate()

		if out.strip() + err.strip() != "":

			print(out)
			print(err)

		if cleanup:
			os.remove(wig)

# def wig_to_bigwig(wig, output_directory):#, cleanup=True):

# 	bigwig = wig.replace(".wig", ".bigwig")

# 	call = f"wigToBigWig {wig} ./{output_directory}/ChromSizes.txt {bigwig}"

# 	print(call)

# 	p = Popen(call.split(), stdout=PIPE, stderr=PIPE, encoding=ENCODING)

# 	out, err= p.communicate()

# 	if out.strip() + err.strip() != "":

# 		print(out)
# 		print(err)

	# if cleanup:
	# 	os.remove(wig)

class Logger(object):
	def __init__(self, file_name):
		self.terminal = sys.stdout
		self.file_name = file_name
		self.log = open(file_name, "w")

		git_dir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
		git_repo = Repo(git_dir)
		git_commit = git_repo.head.commit.tree

		self.terminal.write(f"git_dir:    {git_dir}\n")
		self.terminal.write(f"git_commit: {git_commit}\n\n")
		self.log.write(f"git_dir:    {git_dir}\n")
		self.log.write(f"git_commit: {git_commit}\n\n")

	def clear_ansi(self, message):
		return(message.replace("\033[1m", "").replace("\033[0m",""))

	def write(self, message):
		self.terminal.write(message)
		self.log.write(self.clear_ansi(message))

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
				encoding=ENCODING)
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
# 		p = Popen(call, stdout=PIPE, stderr=PIPE, encoding=ENCODING)
# 		out,err=p.communicate()
# 		# print(out)
# 		# print(err)


# 	# call = f"samtools view -@ 4 -F 4 {bam}"
# 	call = ['samtools', 'view', '-F', '4', bam]


# 	if locus:
# 		call.append(locus)
		
# 	# print(call)
# 	p = Popen(call, stdout=PIPE, stderr=PIPE, encoding=ENCODING)

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

# 	p = Popen(call, stdout=PIPE, stderr=PIPE, encoding=ENCODING)
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


	# print(bam_rgs)


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


def subsample(total_aligned_reads, base_alignment_file, params):

	from random import sample, seed, shuffle

	seed(params['subsample_seed'])

	target_depth = params['subsample']
	seed         = params['subsample_seed']

	seed_string  = f"s{seed}"

	def parse_target_depth():
		multipliers = {'M' : 1000000, "K" : 1000, "G" : 1000000000}

		number = float(re.sub('[A-Za-z]', '', target_depth))
		prefix = re.sub(r'\d', '', target_depth).upper()

		if prefix != '':
			number = round(number * multipliers[prefix])



		if number < 1000:
			subsample_string = str(round(number))
		elif number < 1000000:
			subsample_string = str(round(number/1000)) + 'k'
		else:
			subsample_string = str(round(number/1000000)) + 'M'

		return(number, subsample_string)

	td_number, td_string = parse_target_depth()

	if td_number > total_aligned_reads:
		print("subsample larger than to alignment. Ignoring command....")
		print(f'{total_aligned_reads:,} < {td_number:,} subsample')
		sys.exit()
		return(alignment_file, False, total_aligned_reads)

	nn = floor(total_aligned_reads / td_number)

	subsample_files = [Path(base_alignment_file.parent, f"{base_alignment_file.stem}_{td_string}_n{n}_{seed_string}.bam") for n in range(nn)]


	if params['subsample_n'] >= len(subsample_files):
		sys.exit(f"ERROR:	specified 'subsample_n' ({params['subsample_n']}) >= expected number of subsamples ({nn})")

	if params['subsample_keep_max'] == 'all':
		max_n = nn
	elif params['subsample_keep_max'] < nn:
		max_n = params['subsample_keep_max']
	else:
		max_n = nn

	print(f"{total_aligned_reads:,} -> {td_number:,}")
	print()


	print('looking for subset file(s)')
	all_found = True
	for n in range(max_n):
		file = subsample_files[n]
		print(f"  {isfile(file)}\t{file}")
		if not isfile(file):
			all_found = False



	if all_found and not params['force']:
		print("subset alignments found! skipping...")
		for i,f in enumerate(subsample_files):
			print(f"  n{i}  {isfile(f)}  {f}")
		return subsample_files[params['subsample_n']]


	print()
	print('  sampling...')

	sample_i = sample(range(total_aligned_reads), td_number * max_n)

	print('  sorting...')

	sample_i = sorted(sample_i, reverse=True)

	print(f"  splitting into ({nn}) discrete sub-alignments... (keeping {max_n})")

	sample_n = []
	for n in range(max_n):
		sample_n += [n] * td_number
	shuffle(sample_n)




	this_i = sample_i.pop()
	this_n = sample_n.pop()

	temp_files = []
	for n in range(max_n):
		temp_files.append(Path(base_alignment_file.parent, f"temp_n{n}_s{seed}_{time()}.sam"))

	open_files = []
	for f in temp_files:
		outf = open(f, 'w') 
		open_files.append(outf)

		call = ['samtools', 'view', "-H", base_alignment_file]
		p = Popen(call, stdout=outf, encoding=ENCODING)
		p.wait()
	# fraction = str(number / sum(chrom_depth_c.values()))

	# print(call)




	call = ['samtools', 'view', "-F", '4', base_alignment_file]
	p = Popen(call, stdout=PIPE, encoding=ENCODING)

	for i,line in enumerate(p.stdout):


		if i == this_i:
			open_files[this_n].write(line)

			if not sample_n:
				p.terminate()
				break


			this_i = sample_i.pop()
			this_n = sample_n.pop()

			# print(len(sample_i), '        ', end = '\r')


	p.wait()

	for n in range(max_n):
		open_files[n].close()


	# shutil.move(temp_file, subsample_alignment_file)
	# sys.exit()

	for n in range(max_n):
		file = subsample_files[n]

		with open(file, 'wb') as outf:

			call = ['samtools', 'view', '-h', '--bam', temp_files[n]]
			p = Popen(call, stdout=outf)
			p.wait()

		# print("removing", temp_file)
		os.remove(temp_files[n])

	# print(f'Using {ssamp.files} as alignment for annotation')
	print()

	return subsample_files[params['subsample_n']]



def module_title(module, version):
	print()
	print()
	print(f"{color.BOLD}Module:{color.END} {module}")
	print(f"{color.BOLD}Version:{color.END} {version}")
	print()


