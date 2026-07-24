# generic functions

import sys
import os

import click
from .cli import cli
from click_option_group import optgroup


from subprocess import PIPE, Popen, call, DEVNULL
from pathlib import Path, PurePath
from time import time
from collections import Counter
from os import stat
from math import floor


import pysam
import json
# from pprint import pprint
from glob import glob



ENCODING='cp850'




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

	def add_hmmer(self):
		try:
			p = Popen(['nhmmer', '-h'], stdout=PIPE, stderr=PIPE, encoding=ENCODING)
			out,err = p.communicate()
			out = out.split("\n")
			# print(err)
			found=True
			version = out[1].split()[2]

		except FileNotFoundError:
			found=False
			version=''

		self.reqs.append(('hmmsearch', found, version))

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

	# def add_rnafold(self):
	# 	try:
	# 		p = Popen(['RNAfold', '--version'], stdout=PIPE, stderr=PIPE, encoding=ENCODING)
	# 		out,err = p.communicate()
	# 		out = out.split("\n")
	# 		found=True
	# 		version = out[0].split()[-1]

	# 	except FileNotFoundError:
	# 		found=False
	# 		version=''

	# 	self.reqs.append(('RNAfold', found, version))

	# def add_RNAfold(self):

	# 	try:
	# 		p = Popen(['RNAfold', '--version'], stdout=PIPE, stderr=PIPE, encoding=ENCODING)
	# 		out,err = p.communicate()
	# 		out = out.split("\n")
	# 		found=True
	# 		version = out[0].split()[-1]

	# 	except FileNotFoundError:
	# 		found=False
	# 		version=''

	# 	self.reqs.append(('RNAfold', found, version))

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
			version = out[1].split()[2]

		except FileNotFoundError:
			found=False
			version=''

		self.reqs.append(('prefetch', found, version))

		try:
			p = Popen(['fasterq-dump', '-V'], stdout=PIPE, stderr=PIPE, encoding=ENCODING)
			out,err = p.communicate()
			out = out.split("\n")
			found=True
			version = out[1].split()[2]

		except FileNotFoundError:
			found=False
			version=''

		self.reqs.append(('fasterq-dump', found, version))

	def check(self):

		fail=False

		print("Requirements:")

		for tool, found, version in self.reqs:
			if not found:
				fail = True

			if found:
				print('[x]', tool, "->", version)
			else:
				print('[ ]', tool)

		
		if fail:
			sys.exit("Error: requirements not met!")

		print()

class inputClass():

	def __init__(self, params):
		self.params = params

		try:
			if params['override']:
				self.override = True
			else:
				self.override = False
		except KeyError:
			self.override = False

		self.output_directory = params['output_directory']
		self.output_directory.mkdir(parents=True, exist_ok=True)


		project_name = self.output_directory.name

		self.file = Path(self.output_directory, "inputs.json")

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
			'gene_annotation_file',
			'annotation_files',
			'min_length',
			'max_length'
			]

		self.rep_groups = {}

		self.paths = [
			"untrimmed_libraries",
			"trimmed_libraries",
			"alignment_file",
			'genome_file',
			'jbrowse_directory',
			'gene_annotation_file',
			'annotation_files'
			]


		for i in self.input_list:
			self.inputs[i] = None


		if self.file.is_file():
			try:
				self.read()
			except json.decoder.JSONDecodeError:
				print("inputs.json - DECODER ERROR!")
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

		self.remove_duplicates()

		with open(self.file, 'w') as outf:
			outf.write(json.dumps(self.inputs, indent=2))
		self.decode_inputs()


	def remove_duplicates(self):

		keys = ['srrs','untrimmed_libraries','trimmed_libraries']

		for key in keys:
			if self.inputs[key]:
				self.inputs[key] = list(dict.fromkeys(self.inputs[key]))




	def encode_inputs(self):

		od = self.output_directory

		def encode_path(p):

			## is_relative_to() was added in python 3.9.
			try:
				if p.is_relative_to(od):
					return(str(p.relative_to(od)))
				else:
					return(str(p.absolute()))
			except AttributeError:
				print("Warning: python 3.9+ required for Path.is_relative_to(). Using compatibility function.")
				return(str(p.absolute()))


		for p in ["untrimmed_libraries", "trimmed_libraries", 'annotation_files']:

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


		for p in ["untrimmed_libraries", "trimmed_libraries", 'annotation_files']:
			if p in self.inputs and self.inputs[p]:
				for i in range(len(self.inputs[p])):
					self.inputs[p][i] = decode_path(self.inputs[p][i])


		for p in ["alignment_file", 'genome_file', 'jbrowse_directory', 'gene_annotation_file']:
			if p in self.inputs:
				if self.inputs[p]:
					self.inputs[p] = decode_path(self.inputs[p])
			else:
				self.inputs[p] = None





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

				if isinstance(value, tuple):
					value = list(value)

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

		print("Required options:")
		pass_check = True
		offset = 25

		for option in required_options:
			value = self.inputs[option]

			warnings=''

			if not value:
				pass_check=False
				check_str = " "

			else:

				if isinstance(value, Path):

					if not value.is_file() and not value.is_dir():
						pass_check=False
						check_str = "?"
						warnings=f"\n  Warning:  Path() not found -> {value}"

					else:
						check_str = "x"

				else:
					check_str = "x"


			print(f"[{check_str}] {option}:", "." * (offset-len(option)), value, warnings)


			# 	pass_check = False
			# 	print(f'Error: required option {option} not supplied or in logged inputs...')

		if not pass_check:
			sys.exit("Error: one or more essential options are not provided in 'inputs.json' or program call")

		headers = ["Other Options:", "Other Params:"]

		other_options = list(self.inputs.keys())
		other_options += [o for o in self.inputs.keys() if o not in other_options]
		other_options = [o for o in other_options if o not in required_options]

		other_params = list(self.params.keys())
		other_params = [o for o in other_params if o not in other_options]


		for options in [other_options, other_params]:
			print()
			print(headers.pop(0))
			for option in options:

				if option in self.inputs:
					value = self.inputs[option]
				else:
					value = self.params[option]

				if isinstance(value, Path):
					value = value.relative_to(self.output_directory)




				if type(value) in [list, tuple, set]:
					if not value:
						print(f"    {option}:", "." * (offset-len(option)), str(value))
						continue

					if isinstance(value[0], Path):
						value[0] = value[0].relative_to(self.output_directory)
					print(f"    {option}:", "." * (offset-len(option)), str(value[0]))
					for v in value[1:]:

						if isinstance(v, Path):
							v = v.relative_to(self.output_directory)

						print(f"     ", " " * (offset), str(v))

				elif type(value) == dict:

					first = True
					for k,v in value.items():

						if isinstance(v, Path):
							v = v.relative_to(self.output_directory)

						if first:
							prefix = f"    {option}: " + "." * (offset-len(option)) + " "
							first=False
						else:
							prefix = " " * len(prefix)

						print(f"{prefix}{k} : {v}")

				else:
					print(f"    {option}:", "." * (offset-len(option)), str(value))

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

			if not self.inputs['gene_annotation_file'].is_file():
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

	if not c:
		return(None)

	d = {}

	for cond, libs in c.items():
		for lib in libs:
			d[lib] = cond

	return(d)




def validate_outdir(ctx, param, od):

	if od is None:
		od = Path().cwd()
		if not Path(od, 'inputs.json').is_file():
			sys.exit("Error: cannot run yasma in an uninitialized directiory (doesn't contain inputs.json) without specifying -o/--output_directory.\n\nTo initialize as a yasma directory, run the same command including the directory you want to use as the home directory of the analysis (-o ./path/to/my_directory) or the current directory with (-o .). This directory name will be used as the project name.")



	else:
		od = Path(od)
		if not od.is_dir():
			sys.exit(f"InputError: --output_directory '{od}' not found!")

	if od == Path("."):
		od = Path().cwd()

	return(od)


def validate_glob_path(ctx, param, value):


	if len(value) == 0:
		# print(param)
		# raise click.UsageError("Error: Missing or empty option '-l'")
		return(None)

	input_paths = list(value)


	full_paths = []
	for path in input_paths:
		paths = glob(path)

		for path in paths:

			path = Path(path)
			full_paths.append(path.absolute())

			if not path.is_file() and not path.is_dir():
				print(f"Warning: bad path ({path}) removed!")
				# raise click.BadParameter(f"path not found: {path}")

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
			print(line)
			sys.exit(f"FileTypeError: '{str(path)}' does not start with the expected character '{first_char}'. Are you sure this file is OK?")

		f.close()

	return(paths)


def validate_path(ctx, param, value):

	if not value:
		return(None)

	path = Path(value.strip())

	if not path.is_file() and not path.is_dir():
		raise click.BadParameter(f"path not found: {path}")
		

	full_path = Path(path).absolute()

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


def process_range(unprocessed_peaks):
	peaks = set()
	for peak in unprocessed_peaks:
		if peak.count("-") == 1:
			for r in range(int(peak.split('-')[0]),int(peak.split('-')[1])+1):
				peaks.add(r)

		else:
			try:
				peaks.add(int(peak))
			except ValueError:
				sys.exit(f"Error: '{peak}' if incorrectly formated. Must be a number or range.")


	peaks = list(peaks)
	peaks.sort()
	return(peaks)

def parse_locus(locus):
	locus = locus.replace("..", "-").strip()

	chrom = locus.split(":")[0]
	start = int(locus.split(":")[1].split("-")[0])
	stop  = int(locus.split(":")[1].split("-")[1])

	return(chrom, start, stop)

def samtools_faidx(locus=None, strand=None, genome_file=None):

	if not genome_file:
		sys.exit("genome_file not defined - why is this?")


	if not locus:
		call = ['samtools', 'faidx', genome_file]
		p = Popen(call, stdout=PIPE, stderr=PIPE, encoding=ENCODING)
		p.wait()
		return


	genf = pysam.FastaFile(genome_file)

	out = genf.fetch(region=locus)
	out = out.upper()


	# call = ['samtools', 'faidx', genome_file, locus]
	# # print(" ".join(call))
	# p = Popen(call, stdout=PIPE, stderr=PIPE, encoding=ENCODING)

	# out, err = p.communicate()

	# if err != "":
	# 	print(f"WARNING: {err}")
	# 	# sys.exit(err)

	# out = "".join(out.split("\n")[1:]).strip().upper()

	out = out.replace("T","U")
	# print(out)

	if strand == '-':
		out = out[::-1]
		out = complement(out)
	return(out)



class sizeClass():
	def __init__(self, 
		sizes=[],
		minmax=(15,30)):


		self.min_size=minmax[0]
		self.max_size=minmax[1]

		self.depth = 0

		self.size_c = Counter()

		if len(sizes) > 0:
			self.update(sizes)

	def get_keys(self, size):
		keys = set()

		keys.add((size,))

		keys.add((size-1, size+0,))
		keys.add((size-0, size+1,))

		keys.add((size-2, size-1, size+0,))
		keys.add((size-1, size+0, size+1,))
		keys.add((size-0, size+1, size+2,))

		keys = [k for k in keys if min(k) >= self.min_size or max(k) <= self.max_size]
		return(keys)



	def update(self, sizes):

		# if type(sizes) == int:
		# 	sizes = [sizes]

		if not sizes:
			return

		if type(sizes) == int:
			sizes = [sizes]

		for size in sizes:


			if 15 <= size <= 30:

				self.depth += 1
				self.size_c.update(self.get_keys(size))

				# for mer in [1,2,3]:
				# 	self.size_d[mer].update(self.size_key_d[mer][size])


	def get(self):

		mc = self.size_c.most_common()

		self.size_1_depth, self.size_2_depth, self.size_3_depth = 0,0,0
		self.size_1_key, self.size_2_key, self.size_3_key = None, None, None

		for key, depth in mc:
			if self.size_1_depth == 0 and len(key) == 1:
				self.size_1_key, self.size_1_depth = key, depth

			if self.size_2_depth == 0 and len(key) == 2:
				self.size_2_key, self.size_2_depth = key, depth

			if self.size_3_depth == 0 and len(key) == 3:
				self.size_3_key, self.size_3_depth = key, depth

			if self.size_1_depth * self.size_2_depth * self.size_3_depth > 0:
				break

		# pprint(self.size_c.most_common(10))
		# print(self.depth, "->", round(self.depth/2), "min")
		# print(self.size_1_key, self.size_1_depth, sep="\t")
		# print(self.size_2_key, self.size_2_depth, sep="\t")
		# print(self.size_3_key, self.size_3_depth, sep="\t")

		####################################
		## revision Jul 12 2024
		### these used to just require majority for all... (> 0.5), but that is a really weak standard. For example, to have contiguous sizes make up just a bare majority?
		### I think it is more reasonable to say that it is -> freq(n.sizes) > n.sizes / (n.sizes+1)
		### depth safe guards are still important... otherwise we will get some weird loci
		### size_1 did not have a depth threshold... i have added it to 15 to make sure we're not calling loci with virtually no reads to be selective. (8/15 reads must be one size)
		####################################
		## revision Jul 18 2024
		## On second thought, i have opted for just a bare majority to consider a locus size specific.
		## This increasing threshold looks to hold many loci ~just outside~ of consideration. This makes some sense, where there is probably a single predominant size, and adding in peripheral off-sized reads is unlikely to add 1/6 (1-size) or 1/4 (2-sizes) of total locus abundance. 
		## I also upped the minimums abundances a bit. Seems like high-duplication loci (loci skewed towards one or a few sequences) are a problem and maybe this can help it.	
		####################################

		if self.size_1_depth > self.depth * 0.5 and self.depth > 30:
			sizecall = self.size_1_key

		elif self.size_2_depth > self.depth * 0.5 and self.depth > 45:
			sizecall = self.size_2_key

		elif self.size_3_depth > self.depth * 0.5 and self.depth > 60:
			sizecall = self.size_3_key

		else:
			sizecall = tuple("N")

		self.sizecall = sizecall
		# print(self.depth)
		# print(self.sizecall)
		# sys.exit()

		return(sizecall)


	def __str__(self):
		out = self.get()
		# print(out)
		out = "_".join(map(str, out))

		return(out)


	def __eq__(self, other):

		self.get()
		other.get()

		# print(self.sizecall)
		# print(other.sizecall)


		scall = set(self.sizecall)
		ocall = set(other.sizecall)


		def expand_call(call):
			if len(call) == 3:
				call.add("N")
				return(call)

			if "N" in call:
				return(call)

			call.add(min(call)-1)
			call.add(max(call)+1)

			return(call)


		scall = expand_call(scall)
		ocall = expand_call(ocall)


		common = scall.intersection(ocall)

		if len(common) > 1:
			return True
		elif "N" in common:
			return True
		else:
			return False

	def __add__(self, other):
		self.size_c += other.size_c
		self.depth  += other.depth
		return(self)



class assessClass():
	'''produces a line assessment of a locus, similar to ShortStack3'''

	def __init__(self):

		self.header =  ['Locus','Name','Length','Reads','RPM']
		self.header += ['UniqueReads','FracTop','Strand','MajorRNA','MajorRNAReads','Complexity']
		self.header += ['Gap', 'skew', 'size_1n','size_1n_depth', 'size_2n','size_2n_depth', 'size_3n','size_3n_depth', 'sizecall']



	def format(self, locus, seq_c, strand_c, sizecall, aligned_depth, last_stop, project, annotation):

		name, chrom, start, stop = locus


		### Basic information

		depth = sum(seq_c.values())
		rpm = depth / aligned_depth * 1000000

		gap = start - last_stop

		if depth == 0:
			gff_line = [
			chrom, 'yasma_locus', 'empty', start, stop, '.', '.', '.',
			f'ID={name};project={project};annotation={annotation};depth={depth}']
			result_line = [f"{chrom}:{start}-{stop}", name, stop-start, depth, rpm, 0, 'NA', "NA", 'NA', 'NA', 'NA',
						gap, 'NA', 
						'NA', 'NA',
						'NA', 'NA',
						'NA', 'NA',
						'NA'
			]

			return(result_line, gff_line)

		### ShortStack standard metrics

		unique_reads = len(seq_c.keys())
		frac_top = strand_c["+"] / sum(strand_c.values())

	
		if frac_top > 0.8:
			strand = "+"
		elif frac_top < 0.2:
			strand = "-"
		else:
			strand = "."

		major_rna = seq_c.most_common()[0][0]
		major_rna_depth = seq_c.most_common()[0][1]


		# complexity = unique_reads / depth


		### More derived metrics


		complexity = unique_reads / (stop - start)

		skew = major_rna_depth / depth



		



		frac_top   = round(frac_top,3)
		complexity = round(complexity,3)
		rpm        = round(rpm,3)
		skew       = round(skew, 3)

		sizecall.get()

		result_line = [f"{chrom}:{start}-{stop}", name, stop-start, depth, rpm]
		result_line += [unique_reads, frac_top, strand, major_rna, major_rna_depth, complexity]
		result_line += [
			gap, skew, 
			sizecall.size_1_key, sizecall.size_1_depth,
			sizecall.size_2_key, sizecall.size_2_depth,
			sizecall.size_3_key, sizecall.size_3_depth,
			sizecall
		]


		if 'N' in sizecall.sizecall:
			feature_type = "OtherRNA"
		else:
			feature_type = f"RNA_{sizecall}"

		if start < 1:
			start = 1
		gff_line = [
			chrom, 'yasma_locus',feature_type, start, stop, '.', strand, '.',
			f'ID={name};project={project};annotation={annotation};sizecall={sizecall};depth={depth};rpm={rpm};fracTop={frac_top};complexity={complexity};skew={skew};majorRNA={major_rna}'
		]


		return(result_line, gff_line)


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








# def read_loci(results_file):


# 	# results_file = f"{params['output_directory']}/tradeoff/loci.txt"
# 	with open(results_file, 'r') as f:
# 		header = f.readline().strip().split("\t")
# 		header = [h.lower() for h in header]

# 		for line in f:
# 			line = line.strip().split("\t")
# 			d = dict(zip(header, line))

# 			yield d

def make_depth_file(alignment_file, verbose=True):

	if not alignment_file.is_file():
		sys.exit(f"Error: Cannot make depth file, alignment_file not found\n{str(alignment_file)}")

	header = ['rg','chrom','length','abundance']
	depth_file = alignment_file.with_suffix(".depth.txt")

	try:
		alignment_file.with_suffix(alignment_file.suffix + ".bai").unlink()
	except FileNotFoundError:
		pass

	try:
		alignment_file.with_suffix(alignment_file.suffix + ".csi").unlink()
	except FileNotFoundError:
		pass


	try:
		pysam.index(str(alignment_file))
	except AttributeError as err:
		print(err)
		print()
		sys.exit(f"Error: the alignment file failed basic indexing checks. Are you sure this is a complete bamfile?\n{alignment_file}\n")

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
	print()

def get_global_depth(alignment_file, force=False, aggregate_by=['rg','chrom','length']):

	header = ['rg','chrom','length','abundance']
	depth_file = alignment_file.with_suffix(".depth.txt")



	if not depth_file.is_file() or stat(depth_file).st_size < 50 or force:
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





def samtools_view(bam, rgs='all', contig=None, start=None, stop=None, threads=4, boundary_rule='loose'): #, read_minmax=(15,30)):

	bamf = pysam.AlignmentFile(bam,'rb', threads=threads)

	if contig == '*':
		until_eof = True
	else:
		until_eof = False

	if not bamf.has_index():
		print(f'   index not found for {bam}. Indexing with samtools.')

		pysam.index(str(bam))
		bamf.close()
		bamf = pysam.AlignmentFile(bam,'rb', threads=threads)



	if rgs == 'all':
		rgs = [rgd['ID'] for rgd in bamf.header['RG']]

	rgs = set(rgs)



	if boundary_rule == 'tight' and start and stop:

		for read in bamf.fetch(contig=contig, start=start, stop=stop, until_eof=until_eof):
			if not read.get_tag("RG") in rgs:
				continue

			if read.is_unmapped and contig != '*':
				continue

			read_length = read.infer_read_length()

			if read.overlap != read_length:
				continue

			strand = "-" if read.is_reverse else "+"

			# if read_length >= read_minmax[0] and read_length <= read_minmax[-1]:

			yield(strand, 
				read_length, 
				read.get_tag("MD"), 
				read.reference_start, 
				read.reference_name, 
				read.get_tag("RG"), 
				read.get_forward_sequence().upper().replace("T","U"), 
				read.query_name)
	else:

		for read in bamf.fetch(contig=contig, start=start, stop=stop, until_eof=until_eof):

			if not read.get_tag("RG") in rgs:
				continue

			if read.is_unmapped and contig != '*':
				continue

			seq = read.get_forward_sequence().replace("T","U")
			read_length = len(seq)
			# read_length = read.infer_read_length()

			strand = "-" if read.is_reverse else "+"

			# if read_length >= read_minmax[0] and read_length <= read_minmax[-1]:
			yield(strand, 
				read_length, 
				read.get_tag("MD"), 
				read.reference_start, 
				read.reference_name, 
				read.get_tag("RG"), 
				seq, 
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


class Logger(object):
	def __init__(self, file_name):
		self.terminal = sys.stdout
		self.file_name = file_name
		self.log = open(file_name, "w")


	def write(self, message):

		self.terminal.write(message)
		self.log.write(message)


	def flush(self):
		self.terminal.flush()
		self.log.flush()


	def overwrite_lines(self, n=None, text=None):
		if not n:
			n = text.count("\n")+1
		## used to be called write_over_terminal_lines()
		for r in range(n):
			self.terminal.write("\x1b[1A\x1b[2K")




def RNAfold(seq):
	import RNA

	# assert len(seq) > 0, f"hairpin is length 0\n{seq}"

	# call = ['RNAfold', '--noPS']

	# p = Popen(call,
	# 			  stdout=PIPE,
	# 			stderr=PIPE,
	# 			stdin=PIPE,
	# 			encoding=ENCODING)
	# out, err = p.communicate(f">{time()}\n{seq}")


	# mfe = float(out.strip().split()[-1].strip(")").strip("("))
	# fold = out.strip().split("\n")[2].split()[0]

	fc  = RNA.fold_compound(seq)
	fold, mfe = fc.mfe()


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



def complement(s, dna=False):
	d = {"U":"A", 
	"A":"U", "G":"C", "C":"G", "N":"N"}

	if dna:
		d['T'] = 'A'
		d['A'] = 'T'

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



def get_rg(lib):
	lib = Path(lib)
	extensions = "".join(lib.suffixes)

	while lib.suffix in {'.gz', '.zip', '.t', '.fastq', '.fq', '.fasta', '.fa', '.fna'}:
		lib = lib.with_suffix("")

	return str(lib.name)


def get_library_format(file):

	if file.suffix == ".gz":
		import gzip
		with gzip.open(file, 'rb') as f:
			first_line = f.readline().decode(ENCODING)


	else:
		with open(file, 'r') as f:
			first_line = f.readline()

	if first_line.startswith(">"):
		return ".fa"
	elif first_line.startswith("@"):
		return ".fq"
	elif len(first_line) == 0:
		print(f"Warning: file ({str(file)}) appears do contain zero reads...")
	else:
		print(f"Error: file ({str(file)}) does not look like a fasta or fastq (or gzipped version) file...")
		sys.exit()



# def subsample(total_aligned_reads, base_alignment_file, params, inputs):

# 	from random import sample, seed, shuffle

# 	seed(params['subsample_seed'])

# 	target_depth          = params['subsample']
# 	seed                  = params['subsample_seed']

# 	seed_string  = f"s{seed}"

# 	def parse_target_depth():
# 		multipliers = {'M' : 1000000, "K" : 1000, "G" : 1000000000}

# 		number = float(re.sub('[A-Za-z]', '', target_depth))
# 		prefix = re.sub(r'\d', '', target_depth).upper()

# 		if prefix != '':
# 			number = round(number * multipliers[prefix])



# 		if number < 1000:
# 			subsample_string = str(round(number))
# 		elif number < 1000000:
# 			subsample_string = str(round(number/1000)) + 'k'
# 		else:
# 			subsample_string = str(round(number/1000000)) + 'M'

# 		return(number, subsample_string)

# 	td_number, td_string = parse_target_depth()

# 	if td_number > total_aligned_reads:
# 		print("subsample larger than to alignment. Ignoring command....")
# 		print(f'{total_aligned_reads:,} < {td_number:,} subsample')
# 		sys.exit()
# 		return(alignment_file, False, total_aligned_reads)

# 	nn = floor(total_aligned_reads / td_number)

# 	subsample_files = [Path(base_alignment_file.parent, f"{base_alignment_file.stem}_{td_string}_n{n}_{seed_string}.bam") for n in range(nn)]


# 	if params['subsample_n'] >= len(subsample_files):
# 		sys.exit(f"ERROR:	specified 'subsample_n' ({params['subsample_n']}) >= expected number of subsamples ({nn})")

# 	if params['subsample_keep_max'] == 'all':
# 		max_n = nn
# 	elif params['subsample_keep_max'] < nn:
# 		max_n = params['subsample_keep_max']
# 	else:
# 		max_n = nn

# 	print(f"{total_aligned_reads:,} -> {td_number:,}")
# 	print()


# 	print('looking for subset file(s)')
# 	all_found = True
# 	for n in range(max_n):
# 		file = subsample_files[n]
# 		print(f"  {isfile(file)}\t{file}")
# 		if not isfile(file):
# 			all_found = False



# 	if all_found and not params['force']:
# 		print("subset alignments found! skipping...")
# 		for i,f in enumerate(subsample_files):
# 			print(f"  n{i}  {isfile(f)}  {f}")
# 		return subsample_files[params['subsample_n']]


# 	print()
# 	print('  sampling...')

# 	sample_i = sample(range(total_aligned_reads), td_number * max_n)

# 	print('  sorting...')

# 	sample_i = sorted(sample_i, reverse=True)

# 	print(f"  splitting into ({nn}) discrete sub-alignments... (keeping {max_n})")

# 	sample_n = []
# 	for n in range(max_n):
# 		sample_n += [n] * td_number
# 	shuffle(sample_n)




# 	this_i = sample_i.pop()
# 	this_n = sample_n.pop()

# 	temp_files = []
# 	for n in range(max_n):
# 		temp_files.append(Path(base_alignment_file.parent, f"temp_n{n}_s{seed}_{time()}.sam"))

# 	open_files = []
# 	for f in temp_files:
# 		outf = open(f, 'w') 
# 		open_files.append(outf)

# 		call = ['samtools', 'view', "-H", base_alignment_file]
# 		p = Popen(call, stdout=outf, encoding=ENCODING)
# 		p.wait()
# 	# fraction = str(number / sum(chrom_depth_c.values()))

# 	# print(call)




# 	call = ['samtools', 'view', "-F", '4', base_alignment_file]
# 	p = Popen(call, stdout=PIPE, encoding=ENCODING)

# 	for i,line in enumerate(p.stdout):


# 		if i == this_i:
# 			open_files[this_n].write(line)

# 			if not sample_n:
# 				p.terminate()
# 				break


# 			this_i = sample_i.pop()
# 			this_n = sample_n.pop()

# 			# print(len(sample_i), '        ', end = '\r')


# 	p.wait()

# 	for n in range(max_n):
# 		open_files[n].close()


# 	# shutil.move(temp_file, subsample_alignment_file)
# 	# sys.exit()

# 	for n in range(max_n):
# 		file = subsample_files[n]

# 		with open(file, 'wb') as outf:

# 			call = ['samtools', 'view', '-h', '--bam', temp_files[n]]
# 			p = Popen(call, stdout=outf)
# 			p.wait()

# 		# print("removing", temp_file)
# 		os.remove(temp_files[n])

# 	# print(f'Using {ssamp.files} as alignment for annotation')
# 	print()

# 	return subsample_files[params['subsample_n']]





def module_title(module, version):
	print()
	print()
	print(f"Module:  {module}")
	print(f"Version: {version}")
	print()


