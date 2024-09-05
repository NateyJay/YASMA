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




ENCODING='cp850'



class trackClass():
	def __init__(self, bw_file, chromosomes):
		self.bw = pyBigWig.open(str(bw_file), 'w')
		self.bw.addHeader(chromosomes)

		self.last_start      = 0
		self.interval_length = 1
		self.last_chrom = chromosomes[0][0]
		self.last_val   = 0

		self.chromosomes = chromosomes

	def write(self):
		stop = self.last_start + self.interval_length
		self.bw.addEntries(
						[self.last_chrom], 
						[self.last_start], 
						ends= [stop], 
						values= [float(self.last_val)]
						)


	def add(self, chrom, pos, val):

		if chrom != self.last_chrom:
			self.write()
			self.last_start      = 0
			self.interval_length = 1

		elif pos > self.last_start:

			if val != self.last_val:
				self.write()
				self.last_start = pos
				self.interval_length = 1

			else:
				self.interval_length += 1


		self.last_val   = val
		self.last_chrom = chrom

	def close(self):
		self.write()
		self.bw.close()

	def process_bam(self, bam_file, aligned_depth):
		print('processing .bam to form .bw')
		rpm_factor = 1000000 / aligned_depth

		bamf = pysam.AlignmentFile(str(bam_file),'rb')

		for c,l in self.chromosomes:
			print(" ", c, end = " ")
			depths = bamf.count_coverage(c, quality_threshold=0)


			for i in range(l):

				rpm = round(sum([depths[r][i] for r in range(4)]) * rpm_factor,3)
				self.add(c, i, rpm)

				if i % 1000000 == 0:
					print(".", end = '')

				# print(i, depths[1][i], rpm)
			print()

		# for line in pysam.depth(str(bam_file)):
		# 	print(line)
		# 	sys.exit()


class bigwigClass():
	'''A class to handle producing rpm bigwig files from a counter object c[pos] = depth'''

	def __init__(self, file, total_reads, chromosomes, strand= "+", name=''):
		self.file = str(file)
		# self.file = f"./{output_directory}/Coverages/{file}.wig"
		self.bw = pyBigWig.open(self.file, 'w')

		self.total_reads   = total_reads

		self.bw.addHeader(chromosomes)

		if strand == "+":
			self.strand = 1
		elif strand == "-":
			self.strand = -1 

		self.name= name


	def reset(self, chrom_length):
		# self.last_depth_pos = 1
		self.last_pos = 1
		# self.last_depth = 0
		# self.span = 0
		self.depths = [0] * chrom_length


	def add(self, pos, length):

		for r in range(pos, pos+length+1):
			try:
				self.depths[r] += 1
			except IndexError:
				# print(f"WARNING: position {r:,} is longer than the total chromosome length")
				pass

		self.last_pos = r

	def rle(self, chrom):
		last_depth = 0
		span = 1
		last_pos = 0

		starts = []
		ends   = []
		values = []


		for pos, depth in enumerate(self.depths):

			if depth == last_depth:
				span += 1

			else:


				ends.append(pos)
				starts.append(last_pos)
				values.append(round(last_depth / self.total_reads * 1000000, 8) * self.strand)

				# print(self.name, chrom, last_pos, "->", pos, depth, sep='\t')

				last_depth = depth
				last_pos   = pos
				span       = 1

		if last_pos < pos:
			ends.append(pos)
			starts.append(last_pos)
			values.append(round(last_depth / self.total_reads * 1000000, 8) * self.strand)


		if starts[0] == ends[0]:
			starts.remove(0)
			ends.remove(0)
			values.remove(0)

		self.bw.addEntries([chrom] * len(values), starts, ends=ends, values=values)

	def close(self):
		self.bw.close()



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

	def add_sratools(self):

		try:
			p = Popen(['prefetch', '--version'], stdout=PIPE, stderr=PIPE, encoding=ENCODING)
			out,err = p.communicate()
			out = out.split("\n")
			found=True
			version = out[1].split()[-1]

		except FileNotFoundError:
			found=False
			version=''

		self.reqs.append(('sratools', found, version))

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


		self.inputs = {'project_name' : None}


		self.input_list = [
			"srrs",
			"untrimmed_libraries",
			"trimmed_libraries",
			"adapter",
			"alignment_file",
			# 'annotation_readgroups',
			"conditions",
			"annotation_conditions",
			'genome_file',
			'jbrowse_directory',
			'gene_annotation_file'
			]

		self.rep_groups = {}

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

		self.inputs['project_name'] = project_name

		self.parse(params)
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

		def encode_path(p):

			if p.is_relative_to(od):
				return(str(p.relative_to(od)))
			else:
				return(str(p.absolute()))


		for p in ["untrimmed_libraries", "trimmed_libraries"]:
			if self.inputs[p]:
				for i in range(len(self.inputs[p])):
					self.inputs[p][i] = encode_path(self.inputs[p][i])

		for p in ["alignment_file", 'genome_file', 'jbrowse_directory', 'gene_annotation_file']:
			if self.inputs[p]:
				self.inputs[p] = encode_path(self.inputs[p])




	def decode_inputs(self):
		od = self.output_directory

		def decode_path(p):
			return(Path(od, p))


		for p in ["untrimmed_libraries", "trimmed_libraries"]:
			if self.inputs[p]:
				for i in range(len(self.inputs[p])):
					self.inputs[p][i] = decode_path(self.inputs[p][i])


		for p in ["alignment_file", 'genome_file', 'jbrowse_directory', 'gene_annotation_file']:
			if self.inputs[p]:
				self.inputs[p] = decode_path(self.inputs[p])





	def parse(self, params):

		## list type parameters
		for p in ['trimmed_libraries','untrimmed_libraries', 'srrs', 'annotation_conditions']:
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

			# print(option, "->", value)

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

		if not value:
			return

		# elif saved_value != value and self.override:

		if saved_value == value:
			return

		# print((self.inputs[option], value))
		# print(f"  Override!")
		print(f"  Input value '{option}'")
		print(f"      changed: '{saved_value}' -> ")
		print(f"               '{value}'")
		print()
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


			genf = pysam.FastaFile(self.inputs['genome_file'])

			genome_chromosomes = genf.references

			genf.close()



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

			with pysam.AlignmentFile(self.inputs['alignment_file'], 'rb') as bamf:
				header = bamf.header.to_dict()

			for entry in header['SQ']:
				alignment_chromosomes.add(entry['SN'])


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

		print()
		print("Checking chromosome/scaffold/contig overlaps between inputs...")
		print()
		print('genome', 'align', 'gene', 'chrom/scaffold', sep='\t')
		print("--------------------------------------")
		for chrom in all_chroms:
			print(is_found(chrom, genome_chromosomes), is_found(chrom, alignment_chromosomes), is_found(chrom, gene_annotation_chromosomes), chrom, sep='\t')


	def check_paired_end(self):
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
				else:
					first_pairs.append(input_lib)


			if len(libs) != len(possible_basenames):
				# print("Likely paired ends!!!")
				# print("Use these libraries for the first pairs")

				self.inputs['untrimmed_libraries'] = first_pairs




def reverse_conditions(c):

	d = {}

	for cond, libs in c.items():
		for lib in libs:
			d[lib] = cond

	return(d)






			
def validate_glob_path(ctx, param, value):


	if len(value) == 0:
		# print(param)
		# raise click.UsageError("Error: Missing or empty option '-l'")
		return(None)

	paths = list(value)

	full_paths = []
	for path in paths:
		full_paths.append(str(Path(path).absolute()))

		if not isfile(path) and not isdir(path):
			raise click.BadParameter(f"path not found: {path}")

	full_paths = tuple(full_paths)
	return(full_paths)

def validate_library_paths(ctx, param, value):

	paths = validate_glob_path(ctx, param, value)
	if not paths:
		return

	for path in paths:
		path = Path(path)


		if not set(['.fa','.fna','.fasta']).isdisjoint(set(path.suffixes)):
			print("fasta")
			first_char = ">"
		elif not set(['.fq','.fastq']).isdisjoint(set(path.suffixes)):
			print("fastq")
			first_char = "@"
		else:
			sys.exit(f"FileTypeError: '{str(path)}' does not look like a library (allowed: .fa, .fna, .fasta, .fq, .fastq)")


		if ".gz" in path.suffixes:
			import gzip
			f = gzip.open(str(path), 'rb')
			line = f.readline().decode(ENCODING)

		else:
			f = open(str(path), 'r')
			line = f.readline()

		if not line.startswith(first_char):
			sys.exit(f"FileTypeError: '{str(path)}' does not start with the expected character '{first_char}'. Are you sure this file is OK?")

		f.close()

	return(paths)


def validate_path(ctx, param, value):


	if not value:
		return(None)

	path = value.strip()

	if not isfile(path) and not isdir(path):
		raise click.BadParameter(f"path not found: {path}")

	full_path = str(Path(path).absolute())

	return(full_path)

def validate_condition(ctx, param, input_tuple):

	if not input_tuple:
		return(None)

	d = {}
	rep_groups = set()

	for entry in input_tuple:

		if Path(entry).is_file():

			with open(entry, 'r') as f:
				for line in f:
					val, key = line.strip().split('\t')
					try:
						d[key].append(val)
					except KeyError:
						d[key] = [val]

			continue

		if entry.count(":") != 1:
			raise click.BadParameter(f"conditions groups must be a valid path or contain a single colon: {entry}\nexample: SRR123456:WT")


		val, key = entry.split(":")

		if val in rep_groups:
			raise click.BadParameter(f"library base-names can only be mentioned once: {entry}")

		rep_groups.add(val)

		try:
			d[key].append(val)
		except KeyError:
			d[key] = [val]



	return(d)		


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
		self.last_percent = 0

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
		perc = self.get_percent(self.running)

		if perc:
			self.last_percent = perc

		return(perc)








def read_loci(params):


	results_file = f"{params['output_directory']}/tradeoff/loci.txt"
	with open(results_file, 'r') as f:
		header = f.readline().strip().split("\t")
		header = [h.lower() for h in header]

		for line in f:
			line = line.strip().split("\t")
			d = dict(zip(header, line))

			yield d

def make_depth_file(alignment_file, verbose=True):

	header = ['rg','chrom','length','abundance']
	depth_file = alignment_file.with_suffix(".depth.txt")
	try:
		os.remove(str(alignment_file) + ".bai")
	except FileNotFoundError:
		pass

	c = Counter()

	if verbose:
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

	# print(c.keys())
	with open(depth_file, 'w') as outf:
		print("\t".join(header), file=outf)
		for rg in rgs:
			for chrom in chroms:
				for length in lengths:
					print(rg, chrom, length, c[(rg,length, chrom)], sep='\t', file=outf)

def get_global_depth(alignment_file, force=False, aggregate_by=['rg','chrom','length']):

	header = ['rg','chrom','length','abundance']
	depth_file = alignment_file.with_suffix(".depth.txt")



	if not isfile(depth_file) or stat(depth_file).st_size < 50 or force:
		make_depth_file(alignment_file)
		

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
			if not read.get_tag("RG") in rgs:
				continue

			if read.is_unmapped:
				continue

			read_length = read.infer_read_length()

			if read.overlap != read_length:
				continue

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

			if not read.get_tag("RG") in rgs:
				continue

			if read.is_unmapped:
				continue

			read_length = read.infer_read_length()

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



def get_chromosomes(file):
	chromosomes = []
	rgs = []
	# call = f"samtools view -@ 4 -H {file}"

	with pysam.AlignmentFile(file, 'rb') as bamf:
		header = bamf.header.to_dict()

	# pprint(header['SQ'])
	# sys.exit()

	chromosomes = []
	for entry in header['SQ']:
		chromosomes.append((entry['SN'], entry['LN']))

	rgs = []
	for entry in header['RG']:
		rgs.append(entry['ID'])

	return(chromosomes, rgs)



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

	def write(self, message, terminal_only=False):
		self.terminal.write(message)
		if not terminal_only:
			self.log.write(self.clear_ansi(message))


	def flush(self):
		self.terminal.flush()
		self.log.flush()


	def overwrite_lines(self, n=None, text=None):
		if not n:
			n = text.count("\n")
		## used to be called write_over_terminal_lines()
		for r in range(n):
			self.write("\x1b[1A\x1b[2K", terminal_only = True)



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



def inf_counter():
	i = 1
	while True:
		yield(i)
		i += 1


def subsample(total_aligned_reads, base_alignment_file, params, inputs):

	from random import sample, seed, shuffle

	seed(params['subsample_seed'])

	target_depth          = params['subsample']
	seed                  = params['subsample_seed']

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


