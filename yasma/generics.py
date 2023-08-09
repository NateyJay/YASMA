# generic functions

import sys
import click
from os import stat
from pprint import pprint

# import click

from subprocess import PIPE, Popen, call
# from pathlib import Path

from pathlib import Path
from os.path import isfile, isdir

from time import time, sleep
from collections import Counter, deque
from math import floor
# from itertools import count, chain

from statistics import median, mean

# from Levenshtein import distance

# from timeit import timeit

from math import log10

import pyBigWig

import json



ENCODING='cp850'
# ENCODING=ENCODING



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
			p = Popen(['samtools', 'version'], stdout=PIPE, stderr=PIPE, encoding=ENCODING)
			out,err = p.communicate()
			out = out.split("\n")
			found=True
			version = out[0].split()[-1]

		except FileNotFoundError:
			found=False
			version=''

		self.reqs.append(('samtools', found, version))

	def add_shortstack(self):
		try:
			p = Popen(['ShortStack', '-v'], stdout=PIPE, stderr=PIPE, encoding=ENCODING)
			out,err = p.communicate()
			out = out.split("\n")
			print(out)
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

		output_directory = params['output_directory']
		output_directory = output_directory.rstrip("/")

		self.output_directory = Path(output_directory)
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


		for i in self.input_list:
			self.inputs[i] = None

		if isfile(self.file):
			try:
				self.read()
			except json.decoder.JSONDecodeError:
				pass

		pprint(self.inputs)
		self.parse(params)
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
		od = str(self.output_directory.absolute())

		for p in ["untrimmed_libraries", "trimmed_libraries"]:
			if self.inputs[p]:
				for i in range(len(self.inputs[p])):
					self.inputs[p][i] = self.inputs[p][i].replace(od, "<OUTPUT_DIRECTORY>")

		for p in ["alignment_file", 'genome_file', 'jbrowse_directory', 'gene_annotation_file']:
			if self.inputs[p]:
				self.inputs[p] = self.inputs[p].replace(od, "<OUTPUT_DIRECTORY>")





	def decode_inputs(self):
		od = str(self.output_directory)

		for p in ["untrimmed_libraries", "trimmed_libraries"]:
			if self.inputs[p]:
				for i in range(len(self.inputs[p])):
					self.inputs[p][i] = self.inputs[p][i].replace("<OUTPUT_DIRECTORY>", od)

		for p in ["alignment_file", 'genome_file', 'jbrowse_directory', 'gene_annotation_file']:
			if self.inputs[p]:
				self.inputs[p] = self.inputs[p].replace("<OUTPUT_DIRECTORY>", od)



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

		elif saved_value != value:
			print(f"  Warning: input for option '{color.BOLD}{option}{color.END}' does not match logged value")

			print(f"  Replace: ... '{self.inputs[option]}'")
			print(f"  with: ...... '{value}'")

			res = input("   (y)es or (n)o?\n")

			if res == 'y':
				self.inputs[option] = value



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

		if self.inputs['genome_file'] and not isfile(self.inputs['genome_file']+".fai"):
			print(" ".join(['samtools', 'faidx', self.inputs['genome_file']]))
			p = Popen(['samtools', 'faidx', self.inputs['genome_file']], stdout=PIPE, stderr=PIPE)
			p.wait()

		genome_chromosomes = set()
		if self.inputs['genome_file']:

			with open(self.inputs['genome_file'] + ".fai", 'r') as f:

				for line in f:
					genome_chromosomes.add(line.split()[0])

			# print(genome_chromosomes)


		gene_annotation_chromosomes = set()
		if self.inputs['gene_annotation_file']:

			with open(self.inputs['gene_annotation_file']) as f:

				for line in f:
					if not line.startswith("#"):
						gene_annotation_chromosomes.add(line.split()[0])

			gene_annotation_chromosomes = list(gene_annotation_chromosomes)
			# print(gene_annotation_chromosomes)


		alignment_chromosomes = set()
		if self.inputs['alignment_file']:
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
		# print(path)
		if not isfile(path):
			raise click.BadParameter(f"path not found: {path}")

	full_paths = tuple(full_paths)
	return(full_paths)		


def validate_path(ctx, param, value):


	if not value:
		return(None)

	path = value.strip()

	if not isfile(path):
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

	p = Popen(call, stdout=PIPE, stderr=PIPE, encoding=ENCODING)

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



def get_window_depths(alignment_file, chromosome, window_length):
	i = 0
	window = deque()
	d_out = []

	sam_iter = samtools_view(alignment_file, locus=chromosome)

	read = next(sam_iter)

	while True:

		_, length, _, pos, chromosome, _, _, _ = read

		if pos > i:
			i += 1

			try:
				d = window.popleft()
			except IndexError:
				d = 0

			d_out.append(d)


		elif pos == i:

			for r in range(length + window_length):

				try:
					window[r] += 1
				except IndexError:
					window.append(1)

			try:
				read = next(sam_iter)
			except StopIteration:
				return(d_out)
			# if i > 1000000:
			# 	return(d_out)
		elif pos < i:
			print(pos, i)
			sys.exit("this shouldn't happen")


def get_lambda(alignment_file, chromosomes, window_length, output_directory):

	lambda_file = f"./{output_directory}/Lambdas.txt"
	lambda_d = {}


	if isfile(lambda_file):
		with open(lambda_file, 'r') as f:
			for line in f:
				c,m = line.strip().split('\t')
				m = float(m)
				lambda_d[c] = m

	else:
		window_d = {}


		for c, l in chromosomes:

			window_d[c] = get_window_depths(alignment_file, c, window_length)


			with open(f"./{output_directory}/{c}.dist.txt", 'w') as outf:
				for d in window_d[c]:
					print(d, file=outf)


		print('\n')
		print('chrom','median','mean', sep='\t')
		for c,l in chromosomes:

			window_d[c] = sample(window_d[c],10000)

			print(c, median(window_d[c]), round(mean(window_d[c]),4), sep='\t')
			
			lambda_d[c] = mean(window_d[c])


		with open(lambda_file, 'w') as outf:
			for c,m in lambda_d.items():
				print(c,m, sep='\t', file=outf)


	return(lambda_d)


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





def get_global_depth(output_directory, alignment_file, force=False, aggregate_by=['rg','chrom','length']):
	depth_file = f"./{output_directory}/global_depth.txt"

	header = ['rg','chrom','length','abundance']

	if not isfile(depth_file) or stat(depth_file).st_size < 50 or force:
		call = ['samtools', 'view', '-F', '4', alignment_file]
		print(" ".join(call))
		c = Counter()

		p = Popen(call, stdout=PIPE, encoding=ENCODING)

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





def samtools_view(bam, dcr_range=False, non_range=False, locus=False, rgs=[], boundary_rule='loose'):

	if bam.endswith('.bam'):
		index_file = f"{bam}.bai"
	elif bam.endswith('.cram'):
		index_file = f"{bam}.crai"

	if not isfile(index_file):
		# call = f"samtools index {bam}"
		call= ['samtools','index',bam]
		p = Popen(call, stdout=PIPE, stderr=PIPE, encoding=ENCODING)
		out,err=p.communicate()
		# print(out)
		# print(err)

		# print("WHY AM I INDEXING???")

	if boundary_rule == 'tight':
		lbound = int(locus.split(":")[-1].split("-")[0])
		rbound = int(locus.split(":")[-1].split("-")[1])+1
	else:
		lbound = False


	# call = f"samtools view -@ 4 -F 4 {bam}"
	call = ['samtools', 'view', '-F', '4']
	
	for rg in rgs:
		call += ['-r', rg]

	call += [bam]


	if locus:
		call.append(locus)
		
	# print(call)
	p = Popen(call, stdout=PIPE, stderr=PIPE, encoding=ENCODING)

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


		if lbound and boundary_rule == 'tight':
			if sam_pos >= lbound and sam_pos + length + 1 <= rbound:
				yield(strand, length, size, sam_pos, sam_chrom, rg, seq, read_id)

		else:
			yield(strand, length, size, sam_pos, sam_chrom, rg, seq, read_id)

	p.wait()

def get_chromosomes(file, output_directory=False):
	chromosomes = []
	rgs = []
	# call = f"samtools view -@ 4 -H {file}"

	call = ['samtools','view','-H', file]
	# print(call)

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








