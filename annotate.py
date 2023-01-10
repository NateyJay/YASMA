#!/usr/bin/env python3 

import sys
import os
from pprint import pprint

import click

from subprocess import PIPE, Popen, call
from pathlib import Path

from os.path import isfile, isdir

from time import time, sleep
from collections import Counter, deque
from itertools import count, chain

from statistics import median, mean

from Levenshtein import distance

from timeit import timeit

from math import log10


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


class foldClass():
	def __init__(self, name, seq, alignment_file, locus, strand, mas, output_directory):

		self.name    = name
		self.seq     = seq
		self.locus   = locus
		self.strand  = strand

		self.mas = mas

		self.start = int(locus.split(":")[-1].split("-")[0])
		self.stop  = int(locus.split(":")[-1].split("-")[1])

		self.alignment_file   = alignment_file
		self.output_directory = output_directory


		self.RNAfold()

		self.coor     = []
		self.sequence = []
		self.pairs    = []

		self.lines    = []


		self.read()

		xs = [c[0] for c in self.coor]
		ys = [c[1] for c in self.coor]

		# print(round(max(xs) - min(xs),1), "x range")
		# print(round(max(ys) - min(ys),1), 'y range')
		# print("bounding box:", self.bounding_box)

		self.get_depth()

		self.write()
		
	def get_depth(self):

		if self.strand == "+":
			flag = "0"
		elif self.strand == "-":
			flag = "16"
		else:
			sys.exit("what???")

		p1 = Popen(['samtools', 'view', "-h", self.alignment_file, self.locus], 
			stdout=PIPE, encoding='utf-8')

		p2 = Popen(['samtools', 'depth', '-'], 
			stdin=PIPE, stdout=PIPE, encoding='utf-8')
		# out, err = p1.communicate()

		for line in p1.stdout:
			if line.startswith("@") or line.split("\t")[1] == flag:
				p2.stdin.write(line)

		out, err = p2.communicate()




		# print(out)

		depth_d = {}
		for o in out.strip().split("\n"):
			o = o.strip().split()

			key, val = [int(val) for val in o[1:]]

			depth_d[key] = val

		depths = []
		for r in range(self.start, self.stop + 1):

			try:
				depths.append(depth_d[r])
			except KeyError:
				depths.append(0)


		if self.strand == '-':
			depths = depths[::-1]

		self.depths = depths



	def RNAfold(self):

		temp_name = time()

		call = ['RNAfold']

		p = Popen(call,
					  stdout=PIPE,
					stderr=PIPE,
					stdin=PIPE,
					encoding='utf-8')
		out, err = p.communicate(f">{temp_name}\n{self.seq}")


		self.fold_file = f"{self.output_directory}/Folds/{self.name}.eps"
		os.rename(f"{temp_name}_ss.ps", self.fold_file)

		# sys.exit()


	def read(self):


		with open(self.fold_file, 'r') as f:
			for line in f:
				line = line.strip()

				self.lines.append(line)

				if line.startswith("/sequence"):
					line = f.readline().strip()
					self.sequence = line.rstrip("\\")
					self.lines.append(line)

				if line.startswith("/coor"):
					while True:
						line = f.readline().strip()
						self.lines.append(line)

						if 'def' in line:
							break

						line = line.strip().lstrip("[").rstrip("]").split()
						line = [float(l) for l in line]

						self.coor.append(line)

				if line.startswith("%%BoundingBox"):
					bound = line.split()[1:]
					bound = [int(b) for b in bound]

					self.bounding_box = bound





	def write(self):

		outf = open(f"{self.output_directory}/Folds/{self.name}.eps", 'w')

		for line in self.lines:

			if line == '%%EndComments':

				# leg_x = 72
				# leg_y = 720
				# text_x = 150
				# text_y = 134

				leg_x  = 72
				leg_y  = self.bounding_box[3] + 120
				text_x = 200
				text_y = self.bounding_box[3] + 120

				print(f'''
0 0 1 setrgbcolor
{leg_x} {leg_y} 4 0 360 arc closepath fill stroke
0 0.5 1 setrgbcolor
{leg_x} {leg_y - 10 * 1} 4 0 360 arc closepath fill stroke
0 1 1 setrgbcolor
{leg_x} {leg_y - 10 * 2} 4 0 360 arc closepath fill stroke
0.5 1 0.5 setrgbcolor
{leg_x} {leg_y - 10 * 3} 4 0 360 arc closepath fill stroke
1 1 0 setrgbcolor
{leg_x} {leg_y - 10 * 4} 4 0 360 arc closepath fill stroke
1 0.5 0 setrgbcolor
{leg_x} {leg_y - 10 * 5} 4 0 360 arc closepath fill stroke
1 0 0 setrgbcolor
{leg_x} {leg_y - 10 * 6} 4 0 360 arc closepath fill stroke
1 0 0.5 setrgbcolor
{leg_x} {leg_y - 10 * 7} 4 0 360 arc closepath fill stroke
1 0 1 setrgbcolor
{leg_x} {leg_y - 10 * 8} 4 0 360 arc closepath fill stroke

0.2 0.2 0.2 setrgbcolor
{leg_x + 80} {leg_y - 10 * 8} 4 0 360 arc 1.3 setlinewidth stroke


0 0 0 setrgbcolor
/Helvetica findfont
8 scalefont
setfont
{leg_x + 8} {leg_y - 2} moveto
(10) show
/Helvetica findfont
4 scalefont
setfont
{leg_x + 18} {leg_y + 2} moveto
(0) show

/Helvetica findfont
8 scalefont
setfont
{leg_x + 8} {leg_y - 22} moveto
(10) show
/Helvetica findfont
4 scalefont
setfont
{leg_x + 18} {leg_y - 18} moveto
(1) show

/Helvetica findfont
8 scalefont
setfont
{leg_x + 8} {leg_y - 42} moveto
(10) show
/Helvetica findfont
4 scalefont
setfont
{leg_x + 18} {leg_y - 38} moveto
(2) show

/Helvetica findfont
8 scalefont
setfont
{leg_x + 8} {leg_y - 62} moveto
(10) show
/Helvetica findfont
4 scalefont
setfont
{leg_x + 18} {leg_y - 58} moveto
(3) show

/Helvetica findfont
8 scalefont
setfont
{leg_x + 8} {leg_y - 82} moveto
(>=10) show
/Helvetica findfont
4 scalefont
setfont
{leg_x + 27} {leg_y - 78} moveto
(4) show



/Helvetica findfont
8 scalefont
setfont
{leg_x + 8 + 80} {leg_y - 83} moveto
(Most Abundant Sequence \(MAS\)) show

/Helvetica findfont
8 scalefont
setfont
{leg_x - 4} {leg_y + 10} moveto
(Depth of Coverage) show


% Information at bottom page.


/Helvetica findfont
8 scalefont setfont
{text_x + 20} {text_y - 10 * 0} moveto
(Name:  {self.name}) show


/Helvetica findfont
8 scalefont setfont
{text_x + 11} {text_y - 10 * 1} moveto
(Location:  {self.locus}) show

/Helvetica findfont
8 scalefont setfont
{text_x + 18} {text_y - 10 * 2} moveto
(Strand:  {self.strand}) show

/Helvetica findfont
8 scalefont setfont
{text_x + 0} {text_y - 10 * 3} moveto
(MAS length:  {len(self.mas)} nt) show


/Helvetica findfont
8 scalefont setfont
{text_x + 2} {text_y - 10 * 4} moveto
(Alignments:  {self.alignment_file}) show



''', file=outf)
			elif line == '% switch off outline pairs or bases by removing these lines':
				
				print('''/maplemark { % i r g b maplemark  draw filled circle around base i
  setrgbcolor
  newpath 1 sub coor exch get aload pop
  fsize 2 div 
  0 360 arc closepath fill stroke
} bind def

/borderdraw { % i borderdraw  draw filled circle around base i
  0.2 0.2 0.2 setrgbcolor
  newpath 1 sub coor exch get aload pop
  fsize 2 div 0 360 arc 
  1.3 setlinewidth
  stroke
} bind def

/show5 { % i mark 5-prime end at base i
  0 0 0 setrgbcolor
  newpath 1 sub coor exch get aload pop moveto
  -5 0 rmoveto
  -15 10 rlineto
  -8 0 rmoveto (5') show stroke
} bind def
''', file=outf)

				print("drawoutline", file=outf)


				for i,d in enumerate(self.depths): 

					if d == 0:
						r, g, b = (0.9, 0.9, 0.9)
					else:
						r, g, b = abundance_to_rgb(d)

					print(f"{i+1} {r} {g} {b} maplemark", file=outf)


				print(line, file=outf)

				print("1 show5", file=outf)


				for i in range(self.seq.index(self.mas),self.seq.index(self.mas) + len(self.mas)):
					print(f"{i} borderdraw", file=outf)

			elif line.startswith("%%BoundingBox"):

				line = line.split()
				# line[3] = int(line[3]) + 200
				line[4] = int(line[4]) + 200

				# bound = line.split()[1:]
				# bound = [int(b) for b in bound]

				# bound[-1] += 500
				# bound[-2] += 500

				line = " ".join(map(str, line))


				# line = f""
				# print(line)




			elif line == "drawoutline":
				line = "% " + line

			elif line == "drawpairs":
				line = "% " + line

			print(line, file=outf)



		outf.close()
		os.remove(self.fold_file)
		


def get_rg_depth(output_directory, alignment_file):
	depth_file = f"./{output_directory}/RGdepth.txt"

	if isfile(depth_file): #and not force:
		with open(depth_file, 'r') as f:

			line = f.readline()
			c = Counter()
			
			for line in f:
				line = line.strip().split("\t")


				c[(line[0], int(line[1]))] += int(line[3])



		return(c)



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


def complement(s):
	d = {"U" : "A", 
	"A":"U", "G":"C", "C":"G", "N":"N"}

	s = "".join([d[letter] for letter in s])
	return(s)


def samtools_view(bam, dcr_range=False, non_range=False, locus=False):

	if not isfile(f"{bam}.bai"):
		# call = f"samtools index {bam}"
		call= ['samtools','index',bam]
		p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')
		out,err=p.communicate()
		# print(out)
		# print(err)


	# call = f"samtools view -@ 4 -F 4 {bam}"
	call = ['samtools', 'view', '-F', '4', bam]


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



		yield(strand, length, size, sam_pos, sam_chrom, rg, seq)

	p.wait()

def get_chromosomes(file,output_directory):
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
		self.log.write(self.clear_ansi(message))

	def flush(self):
		self.terminal.flush()
		self.log.flush()



def inf_counter():
	i = 1
	while True:
		yield(i)
		i += 1





@click.group()
def cli():
	pass


# @click.command()
# @click.option('--count', default=1, help='Number of greetings.')
# @click.option('--name', prompt='Your name',
#               help='The person to greet.')
# def hello(count, name):
#     """Simple program that greets NAME for a total of COUNT times."""
#     for x in range(count):
#         click.echo(f"Hello {name}!")

@click.command()

@click.option("-a", "--alignment_file", 
	required=True, 
	type=click.Path(exists=True),
	help='Alignment file input (bam or cram).')

@click.option('-r', '--annotation_readgroups', 
	required=True,
	multiple=True,
	help="List of read groups (RGs, libraries) to be considered for the annotation. 'ALL' uses all readgroups for annotation, but often pertainent RGs will need to be specified individually.")

@click.option("-o", "--output_directory", 
	default=f"Annotation_{round(time())}", 
	type=click.Path(),
	help="Directory name for annotation output")

@click.option("-f", "--force",
	is_flag=True,
	help='Force remake of supporting files')

def precheck(alignment_file, annotation_readgroups, output_directory, force):
	"""Runs precheck to identify likely dicer sizes."""

	Path(output_directory).mkdir(parents=True, exist_ok=True)


	log_file = f"{output_directory}/Log_precheck.txt"

	message = f"log_file is already exists ({log_file}). The annotator will not over-write by default (use --force to override). Be warned: this will trigger the overwrite of some files in this folder!"
	assert not isfile(log_file) or force, message

	sys.stdout = Logger(log_file)


	chromosomes, bam_rgs = get_chromosomes(alignment_file, output_directory)
	annotation_readgroups = check_rgs(annotation_readgroups, bam_rgs)


	# def get_rg_depth():
	# 	depth_file = f"./{output_directory}/readgroup_depth.txt"
	# 	# if isfile(depth_file): #and not force:
	# 	# 	with open(depth_file, 'r') as f:
	# 	# 		line = f.readline()
	# 	# 		line = line.split("\t")

	# 	# 		if line[0] == file:
	# 	# 			# print(f'\nread depth from {depth_file}...')
	# 	# 			return(int(line[1]))



	# 	print('reading annotation RG depth...')

	# 	c = Counter()

	# 	call = ['samtools', 'view', '-F', '4']

	# 	for r in annotation_readgroups:
	# 		call += ['-r', r]

	# 	call += [alignment_file]
	# 	# print(call)

	# 	p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')

	# 	depth = 0
	# 	for line in p.stdout:
	# 		line = line.strip().split("\t")

	# 		length = int(line[5][:-1])

	# 		c.update([length])

	# 		depth += 1

	# 		# if depth > 1000000:
	# 		# 	p.kill()
	# 		# 	break

	# 	p.wait()

	# 	highest = max(c.values())

	# 	print("length\tprop\tprop_highest\tabundance")
	# 	for r in range(15,31):
	# 		prop = round(c[r] / depth, 4)
	# 		prop_highest = round(c[r] / highest, 4)
	# 		print(r, prop, prop_highest, f"{c[r]:,}", sep='\t')

	# 	with open(depth_file, 'w') as outf:
	# 		print(alignment_file, depth, sep='\t', file=outf)

	depth_c = get_rg_depth(output_directory, alignment_file)


	len_c = Counter()
	for key in depth_c.keys():
		if key[0] in annotation_readgroups:
			len_c[key[1]] += depth_c[key]

	print(f"annotation_readgroups: {annotation_readgroups}")

	print()
	print("length\tp_depth\tp_highest\tabundance")

	dep = sum(len_c.values())
	highest = len_c.most_common(1)[0][1]

	for r in range(15, 31):
		print(r, 
			round(len_c[r] / dep, 4),
			round(len_c[r] / highest, 4),
			len_c[r],
			sep='\t')

	# pprint(depth_c)






@click.command()
@click.option("-a", "--alignment_file", 
	required=True, 
	type=click.Path(exists=True),
	help='Alignment file input (bam or cram).')

@click.option('-r', '--annotation_readgroups', 
	required=True,
	multiple=True,
	help="List of read groups (RGs, libraries) to be considered for the annotation. 'ALL' uses all readgroups for annotation, but often pertainent RGs will need to be specified individually.")

@click.option("-d", "--dicercall", 
	default=[20,21,22,23,24],
	multiple=True, 
	help='List of sRNA lengths that derive from dicer.')

@click.option("-o", "--output_directory", 
	default=f"Annotation_{round(time())}", 
	type=click.Path(),
	help="Directory name for annotation output")

@click.option("-f", "--force",
	is_flag=True,
	help='Force remake of supporting files')

@click.option("--partial_wigs",
	is_flag=True,
	help='Only make wiggle files associated with essential functions (ignoring size and strand specific coverages. (May improve speed)')

@click.option("--window",
	default=100,
	help="Window size (centered on position) for determining DCR vs non-DCR read ratio (counting overlapping reads).")

@click.option("--merge_dist",
	default=150,
	help="Maximum gap size between valid regions to merge to a single locus.")

@click.option('--pad',
	default=10,
	help='Number of bases arbitrarily added to either end of a defined locus.')

@click.option('--rpm_cutoff',
	default=1.0,
	help='RPM depth threshold for DCR-sized reads to be considered as a valid region.')

@click.option('--extension_ratio',
	default=0.5,
	help='Fraction of RPM threshold to be considered for extending a locus boundaries')

@click.option('--dicer_ratio',
	default=3.0,
	help='Ratio of dicer to non-dicer reads to be considered for a valid region')

def annotate(alignment_file, annotation_readgroups, dicercall, output_directory, force, partial_wigs, window, merge_dist, pad, rpm_cutoff, extension_ratio, dicer_ratio):
	"""Main annotation suite."""
	print('run annotation')




	assert isfile(alignment_file), f"{alignment_file} does not exist"


	Path(output_directory).mkdir(parents=True, exist_ok=True)
	Path(f'./{output_directory}/Coverages').mkdir(parents=True, exist_ok=True)




	log_file = f"{output_directory}/Log_annotation.txt"
	message = f"log_file is already exists ({log_file}). The annotator will not over-write by default (use --force to override). Be warned: this will trigger the overwrite of some files in this folder!"
	assert not isfile(log_file) or force, message

	if force and isfile(log_file):
		print("force flag included, overwrite possible", end='', flush=True)
		counter = 0
		while counter < 5:
			print('.', end='', flush=True)
			sleep(1)
			counter += 1
		

	sys.stdout = Logger(log_file)



	dicercall = [int(d) for d in dicercall]
	dcr_range = set([r for r in range(min(dicercall), max(dicercall) + 1)])
	non_range = set([r for r in range(15,30) if r not in dcr_range])


	assert window % 2 == 0, "Window must be an even number!"
	half_window = int(window / 2)



	# possible names:
	# DicerLocus
	# smallDicer


	print()
	print()
	print("\033[1m-- annotator v0.3x --\033[0m")

	print()
	print()
	print(f"\033[1m[Prerequisites]\033[0m")


	def check_reqs():
		tool_responses = {
		'samtools version' : 'Samtools compilation details:',
		# 'gt --version' : 'gt (GenomeTools)',
		# 'bgzip --version' : 'bgzip (htslib)',
		# 'tabix --version' : 'tabix (htslib)',
		'wigToBigWig' : 'wigToBigWig v 2.8',
		}


		fails = []

		for tool, response in tool_responses.items():
			p = Popen(tool.split(), stdout=PIPE, stderr=PIPE, encoding='utf-8')

			out, err = p.communicate()

			merged = out + err

			tool = tool.split()[0]


			# print(out)
			# print(err)
			if response in merged:
				pass_str = "[x]"
			else:
				pass_str = "[ ]"
				fails.append(tool)




			print(" ", pass_str, tool)
			# sys.exit()

		# do_not_prepare_gff = False
		do_not_make_bigwig = False

		if 'samtools' in fails:
			sys.exit("Error: samtools not found in PATH (required)")

		# for tool in ['gt','bgzip','tabix']:
		# 	if tool in fails:
		# 		do_not_prepare_gff = True
		# 		break

		if 'wigToBigWig' in fails:
			do_not_make_bigwig = True

		# if do_not_prepare_gff:
		# 	print("Warning: will not prepare indexed gff for jbrowse due to missing reqs")
		if do_not_make_bigwig:
			print("Warning: will not prepare bigwig files due to missing reqs")

		return(do_not_make_bigwig)


	do_not_make_bigwig = check_reqs()

	# rpm_cutoff = round(rpm_cutoff / window, 6)

	chromosomes, bam_rgs = get_chromosomes(alignment_file, output_directory)
	annotation_readgroups = check_rgs(annotation_readgroups, bam_rgs)





	## initiating output files
	gff_file   = f"{output_directory}/Annotation.gff3"

	with open(gff_file, 'w') as outf:
		print("##gff-version 3", file=outf)

		for c, l in chromosomes:
			print(f"##sequence-region   {c} 1 {l}", file=outf)


	count_file = f'{output_directory}/Counts.txt'
	with open(count_file, 'w') as outf:
		print("cluster", 'ann_depth', 'tot_depth', "\t".join(bam_rgs), sep='\t', file=outf)



	results_file = f"{output_directory}/Results.txt"
	with open(results_file, 'w') as outf:
		print("#name\tlocus\tlength\tgap\tdepth\trpm\tdepth:length\tfrac_top\tstrand\tfrac_dicer\tdcr_reads\tnon_reads\tdicercall\tfrac_dicercall\t" + "\t".join(map(str, dcr_range)), file=outf)


	reads_file = f"{output_directory}/TopReads.txt"
	with open(reads_file, 'w') as outf:
		print("cluster\tseq\trank\tdepth\trpm\tlocus_prop", file=outf)





	print()
	print(f"\033[1m[General settings]\033[0m")
	print(f"             alignment_file: {alignment_file}")
	print(f"      output_directory: {output_directory}")
	print(f" annotation_readgroups: {list(annotation_readgroups)}")
	print(f"           dicer_sizes: {list(dcr_range)}")
	print(f"                 force: {force}")
	print(f"          partial_wigs: {partial_wigs}")
	print(f"              log_file: {log_file}")
	print()


	print(f"\033[1m[Annotation settings]\033[0m")
	print(f"     window: {window}")
	print(f" merge_dist: {merge_dist}")
	print(f"        pad: {pad}")



	# def get_library_depth(output_directory, file):
	# 	depth_file = f"./{output_directory}/library_depth.txt"
	# 	if isfile(depth_file): #and not force:
	# 		with open(depth_file, 'r') as f:
	# 			line = f.readline()
	# 			line = line.split("\t")

	# 			if line[0] == file:
	# 				# print(f'\nread depth from {depth_file}...')
	# 				return(int(line[1]))



	# 	print('reading annotation RG depth...')

	# 	call = ['samtools', 'view', '-F', '4']

	# 	for rg in annotation_readgroups:
	# 		call += ['-r', rg]

	# 	call += [file]

	# 	p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')

	# 	depth = 0
	# 	for line in p.stdout:
	# 		depth += 1

	# 	p.wait()


	# 	with open(depth_file, 'w') as outf:
	# 		print(file, depth, sep='\t', file=outf)

	# 	return(depth)

	# library_depth = get_library_depth(output_directory, alignment_file)

	depth_c = get_rg_depth(output_directory, alignment_file)

	library_depth = 0
	for key in depth_c.keys():
		if key[0] in annotation_readgroups:
			library_depth += depth_c[key]


	read_equivalent = 1 / library_depth * 1000000
	depth_cutoff = library_depth / rpm_cutoff / 1000000
	ext_cutoff = rpm_cutoff * extension_ratio


	print()
	print('\033[1m[Depth settings]\033[0m')
	print(f'    ann_rg_depth: {library_depth:,} reads')
	print(f'          1 read: {round(read_equivalent,5)} rpm')
	print(f"      rpm_cutoff: {rpm_cutoff} rpm -> {round(depth_cutoff,2)} reads")
	print(f"      ext_cutoff: {ext_cutoff} rpm > {round(depth_cutoff*extension_ratio,2)} reads")
	print(f" extension_ratio: {extension_ratio}")
	print(f"     dicer_ratio: {dicer_ratio}")





	class wiggleClass():
		def __init__(self, file):
			self.file = f"./{output_directory}/Coverages/{file}.wig"
			self.outf = open(self.file, 'w')
			self.reset()


		def reset(self):
			self.val = 0
			self.start_pos = 1


		def add(self, val, pos, chrom):


			if val != self.val:
				span = pos - self.start_pos

				if span > 0:

					print(f"variableStep chrom={chrom} span={span}", file=self.outf)
					print(f"{self.start_pos} {self.val}", file=self.outf)

					self.val = val
					self.start_pos = pos

		def convert(self, cleanup=False):

			self.outf.close()

			wig = self.file

			bigwig = wig.replace(".wig", ".bigwig")

			print(f"  {wig} -> {bigwig}", flush=True)

			call = f"wigToBigWig {wig} ./{output_directory}/ChromSizes.txt {bigwig}"

			p = Popen(call.split(), stdout=PIPE, stderr=PIPE, encoding='utf-8')

			out, err= p.communicate()

			if out.strip() + err.strip() != "":

				print(out)
				print(err)

			if cleanup:
				os.remove(wig)

	coverage_names = ['dcr','non']

	if not partial_wigs:
		for s in ["+","-"]:
			for l in list(dcr_range) + ['non']:
				coverage_names.append(f"{l}{s}")



	wig_d = {c : wiggleClass(c) for c in coverage_names + ['rpm_passing', 'ratio_passing', 'passing_all']}


	cluster_counter = inf_counter()



	sam_iter = samtools_view(alignment_file, dcr_range, non_range)
	read = next(sam_iter)
	sam_strand, sam_length, sam_size, sam_pos, sam_chrom, sam_rg, sam_read = read


	total_locus_count = 0
	chrom_count = 0

	class locusClass():
		def __init__(self):
			self.reads = deque([])


			self.in_locus = False
			self.nucleated = False
			self.last_hit_pos = 0
			self.last_end = 0
			self.start = False
			self.stop  = False

		def hit(self, pos, hit_type):
			self.last_hit_pos = pos

			if not self.in_locus:
				self.start = pos

			self.in_locus = True


			if hit_type == 'nuc':
				self.nucleated = True


		def add(self, read):
			self.reads[-1].append(read)

		def check(self, pos):
			# print(self.reads)

			# if pos > 1392000:
			# 	print(pos, self.last_hit_pos)

			if self.in_locus:

				if pos - self.last_hit_pos > merge_dist:
					self.in_locus = False

					self.stop = self.last_hit_pos

					if self.nucleated:
						self.nucleated = False
						return(True)

			else:

				while True:

					if len(self.reads) == 0:
						break


					if len(self.reads[0]) == 0:
						self.reads.popleft()
					elif self.reads[0][0][0] + 35 < pos - merge_dist:
						self.reads.popleft()

					else:
						break
						
			self.reads.append([])






		def summarize(self, chrom):

			# start = self.reads[0][1]
			start = self.start# - self.pad
			stop  = self.stop# + self.pad


			# pprint(self.reads)
			# print(self.reads.values())

			reads = chain.from_iterable(self.reads)


			reads = [r for r in reads if r[0] + r[1] >= start and r[0] <= stop]

			len_c    = Counter()
			size_c   = Counter()
			strand_c = Counter()
			rg_c     = Counter()
			read_c   = Counter()

			read_starts = []
			read_stops  = []


			for r in reads:
				sam_pos, sam_length, sam_size, sam_strand, sam_rg, sam_read = r

				read_starts.append(sam_pos)
				read_stops.append(sam_pos + sam_length)

				len_c.update([sam_length])
				size_c.update([sam_size])
				strand_c.update([sam_strand])
				rg_c.update([sam_rg])
				if sam_strand == "-":
					sam_read = complement(sam_read[::-1])
				read_c.update([sam_read])


			# print(start, stop)
			start = min(read_starts) - pad
			stop  = max(read_stops) + pad


			name = f"Cl_{next(cluster_counter)}"


			if len(reads) == 0:
				print(f"WARNING: {name} detected no reads. Likely an error. Skipping.")
				return(0,0)

			# pprint(len_c)
			# pprint(read_c)


			dist_to_last = start - self.last_end
			self.last_end = stop


			if start < 0:
				start = 0

			coords = f"{chrom}:{start}..{stop}"
			length = stop - start
			n_reads = len(reads)
			frac_top = round(strand_c["+"] / n_reads,3)
			frac_dicercall = round(size_c['dcr'] / n_reads, 3)
			rpm = round(n_reads / library_depth * 1000000, 2)


			cum_count = 0
			top_reads = read_c.most_common(100)
			with open(reads_file, 'a') as outf:
				for rank, read in enumerate(top_reads):

					seq, dep = read
					rpm = round(dep * read_equivalent, 4)

					cum_count += dep

					loc_prop = round(cum_count / n_reads, 4)

					print(name, seq, rank, dep, rpm, loc_prop, file=outf, sep='\t')

					if loc_prop >= 0.3:
						break




			predominant_length, predominant_length_depth = len_c.most_common(1)[0]
			predominant_length_depth = round(predominant_length_depth/n_reads,3)
			# print(predominant_length, predominant_length_depth)

			if frac_top >= 0.8:
				strand = '+'
			elif frac_top <= 0.2:
				strand = "-"
			else:
				strand = '.'

			depth_by_length = round(n_reads / length, 3)


			to_print = [name, coords, length, dist_to_last, n_reads, rpm, depth_by_length, frac_top, strand]
			to_print += [frac_dicercall, size_c['dcr'],  size_c['non']]
			to_print += [predominant_length, predominant_length_depth]
			to_print += [len_c[d] for d in dcr_range]
			to_print = "\t".join(map(str, to_print))


			with open(results_file, 'a') as outf:
				print(to_print, sep='\t', file=outf)

			# print(name, coords, length, dist_to_last, n_reads, rpm, depth_by_length, frac_top, strand, 
			# 		frac_dicercall, size_c['dcr'],  size_c['non'], 
			# 		predominant_length, round(predominant_length_depth/n_reads,3), "\t".join([str(len_c[d]) for d in dcr_range]), sep='\t')


			# sys.exit()

			with open(gff_file, 'a') as outf:
				print(f"{chrom}\tsmoothLocus\tnc_RNA\t{start}\t{stop}\t.\t.\t.\tID={name};dicercall={predominant_length};frac_dicercall={predominant_length_depth}", file=outf)


			to_print = [name]

			to_print.append(sum([rg_c[rg] for rg in annotation_readgroups]))
			to_print.append(sum([rg_c[rg] for rg in bam_rgs]))
			to_print += [rg_c[rg] for rg in bam_rgs]
			to_print = "\t".join(map(str, to_print))

			with open(count_file, 'a') as outf:
				print(to_print, file=outf)


			# self.reads = {}

			return(length, dist_to_last)

	class coverageClass():
		def __init__(self, bandwidth=0):
			self.ls = deque()
			self.bandwidth = bandwidth

		def get(self):
			try:
				d = self.ls.popleft()
			except IndexError:
				d = 0

			return(d)

		def add(self, length):
			for r in range(length + self.bandwidth + 1):
				try:
					self.ls[r] += 1
				except IndexError:
					self.ls.append(1)



	def test1(rpm):
		if rpm == None:
			return('-')

		if rpm >= rpm_cutoff:
			return('n')
		elif rpm >= ext_cutoff:
			return('e')
		else:
			return('-')

	def test2(dcr, non):
		if dcr == None or non == None:
			return('-')

		if dcr >= non * dicer_ratio and dcr > 0:
			return('x')
		else:
			return('-')


	for chrom, chrom_length in chromosomes:
		chrom_count += 1


		print()
		print()
		print(f"{chrom_count} / {len(chromosomes)}")
		print(f"chrom: {chrom}")
		print(f"       {chrom_length:,} bp")
		pos = 0

		print("  ", end='')

		locus_lengths = []
		locus_gaps = []



		coverages = {c : coverageClass() for c in coverage_names}
		coverage_buffer = {c : deque([0]*half_window) for c in coverage_names}

		window_coverages = {'dcr' : coverageClass(window), 'non' : coverageClass(window)}
		locus = locusClass()#merge_dist, chrom, pad, cluster_counter, library_depth, read_equivalent, dcr_range, annotation_readgroups, bam_rgs)



		while sam_chrom == chrom:

			# if pos == 1393300:
			# 	sys.exit()

			# if pos == 1000000:
			# 	sys.exit("timeup!")


			corrected_pos = pos - half_window


			if pos == sam_pos:

				if sam_size:

					if sam_rg in annotation_readgroups:

						coverages[sam_size].add(sam_length)
						window_coverages[sam_size].add(sam_length)

						if not partial_wigs:
							if sam_size == 'non':
								coverages[f'{sam_size}{sam_strand}'].add(sam_length)
							else:
								coverages[f'{sam_length}{sam_strand}'].add(sam_length)



					locus.add((sam_pos, sam_length, sam_size, sam_strand, sam_rg, sam_read))

				try:
					read = next(sam_iter)
					sam_strand, sam_length, sam_size, sam_pos, sam_chrom, sam_rg, sam_read = read
					# print(pos, sam_id, sam_size, sep='\t')
				except StopIteration:
					break



			elif pos < sam_pos:

				read_count = {}
				# dens_rpm = {}
				# dens = {}
				win_cov = {}
				rpms = {}



				for size in coverage_names:

					cov = coverages[size].get()
					coverage_buffer[size].append(cov)
					coverage_buffer[size].popleft()

					cov = coverage_buffer[size][0]

					if size[-1] == "-":
						cov = cov * -1




					wig_d[size].add(round(cov * read_equivalent,4), corrected_pos, chrom)

					if size in ['dcr','non']:
						win_cov[size] = window_coverages[size].get()
						rpms[size] = round(win_cov[size] * read_equivalent, 4)


				t1 = test1(round(coverage_buffer['dcr'][0] * read_equivalent, 4))
				t2 = test2(win_cov['dcr'], win_cov['non'])

				tests = f"{t1}{t2}"


				# print(tests)

				if t1 == 'n':
					wig_d['rpm_passing'].add(1, corrected_pos, chrom)

				elif t1 == 'e':
					wig_d['rpm_passing'].add(0.3, corrected_pos, chrom)

				else:
					wig_d['rpm_passing'].add(0, corrected_pos, chrom)


				if not win_cov['dcr'] or win_cov['dcr'] == 0:
					ratio = 0
				else:
					try:
						ratio = round(win_cov['dcr'] / win_cov['non'], 2)
					except ZeroDivisionError:
						ratio = dicer_ratio


				wig_d['ratio_passing'].add(ratio, corrected_pos, chrom)


				if tests == "nx":
					locus.hit(corrected_pos, 'nuc')
					wig_d['passing_all'].add(1, corrected_pos, chrom)

				elif tests == "n-" or t1 == 'e':
					locus.hit(corrected_pos, 'ext')
					wig_d['passing_all'].add(0.3, corrected_pos, chrom)

				else:
					wig_d['passing_all'].add(0, corrected_pos, chrom)



				# if tests and "-" not in tests:
				# if rds['dcr'] and rds['dcr'] > 0 and pos-half_window > 0:
				# if pos >

				# if pos > 1392000:
				# 	if "n" in tests or 'e' in tests:
				# 		print(chrom, pos-half_window, 
				# 			"||", coverage_buffer['dcr'][0], round(win_cov['dcr'], 4), 
				# 			"||", coverage_buffer['non'][0], round(win_cov['non'], 4),
				# 			'||', tests,
				# 			sep='\t')




				if locus.check(corrected_pos):
					length, gap = locus.summarize(chrom)

					locus_lengths.append(length)
					locus_gaps.append(gap)
					total_locus_count += 1


				pos += 1


				if pos % 100000 == 0:
					print(".", end='', flush=True)
				if pos % 1000000 == 0:
					print(" ", end='', flush=True)




		for key in wig_d.keys():
			wig_d[key].add(0, pos, chrom)
			wig_d[key].reset()
		# for size in ['dcr','non']:
		# 	wig_densities[size].add(0, pos)
		# wig_rpm_pass.add(0, pos)
		# wig_pass.add(0, pos)

		locus_count = len(locus_gaps)
		med_length = median(locus_lengths)
		med_gap = median(locus_gaps)
		print()
		print(f"  {locus_count:,} loci found")
		print(f"  {med_length:,} median length")
		print(f"  {med_gap:,} median gap")

		# break


	print()
	print(f"{total_locus_count:,} loci found in total")

	print()

	if not do_not_make_bigwig:
		print("converting wigs to bigwigs...")
		for key in wig_d.keys():
			wig_d[key].convert()
	else:
		print("Not making bigwig files due to missing req...")



	def prepare_gff(gff_input):

		sorted_input = gff_input.replace(".gff3", ".sorted.gff3")
		zipped_input = sorted_input.replace(".gff3", ".gff3.gz")


		print("  sorting...")
		c2 = f"gt gff3 -retainids -sortlines -tidy {gff_input}"
		with open(sorted_input, 'w') as f:
			c2 = Popen(c2.split(), stdout=f)

		print("  zipping...")
		c3 = f"bgzip -f {sorted_input}"
		call(c3.split())

		print("  indexing...")
		c4 = f"tabix -f -p gff {zipped_input}"
		call(c4.split())




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


class hairpinClass():
	def __init__(self, stranded, short_enough, name, locus, strand, input_mas, genome_file, alignment_file, output_directory):


		self.valid   = False
		self.status  = []

		self.stranded = stranded
		self.short_enough = short_enough


		self.name = name
		self.locus = locus
		self.strand = strand
		self.input_mas = input_mas
		self.genome_file = genome_file
		self.alignment_file = alignment_file
		self.output_directory = output_directory.rstrip("/")

		self.chrom = self.locus.split(":")[0]
		self.start = int(self.locus.split(":")[-1].split("-")[0])
		self.stop  = int(self.locus.split(":")[-1].split("-")[1])

		self.length = self.stop - self.start

		self.seq = '-'
		self.fold = '-'
		self.mfe = '-'
		self.pairing = '-'
		self.pos_d = '-'
		self.input_mas_coords = '-'


		self.mas = '-'
		self.star = '-'
		self.duplex_mas = '-'
		self.duplex_fold = '-'
		self.duplex_star = '-'


		self.ruling_d = {
		'mfe_per_nt' : '-',
		'mismatches_total' : '-',
		'mismatches_asymm' : '-',
		'no_mas_structures' : '-',
		'no_star_structures' : '-',
		'precision' : '-',
		'star_found' : '-'
		}

		self.ruling = '-'



		if not stranded:
			self.status.append("hairpin not stranded")
			return

		if not short_enough:

			self.status.append("hairpin too long")
			return


		self.seq, self.fold, self.mfe, self.pairing, mas_d, self.pos_d, self.input_mas_coords = self.get_locus(locus, strand, input_mas)


		self.mas_c = mas_d['all']
		self.mas = self.mas_c.most_common(1)[0][0]


		self.star = '-'
		self.duplex_mas, self.duplex_fold, self.duplex_star = '-','-','-'



		## this needs to be salvaged to work on reproducibility.
		# for rg in mas_d.keys():

		# 	if rg != 'all':
		# 		print(rg)

		# 	mas_c = mas_d[rg]
		# 	mas = mas_c.most_common(1)[0][0]



		if self.mas not in self.seq:
			self.status.append("MAS not found in hairpin sequence")
			return


		self.mas_positions = [r + self.seq.index(self.mas) for r in range(len(self.mas) + 1)]

		self.mas_structures = self.find_secondary_structures("".join([self.fold[p] for p in self.mas_positions]))

		if self.mas_structures:
			self.status.append("secondary structure found in MAS")
			return

		self.star_found = self.find_star()

		if self.star_found:

			# print(self.star)
			# print(" " * self.seq.index(self.star) + self.star)



			self.star_structures = self.find_secondary_structures("".join([self.fold[p] for p in self.star_positions]))

			if self.star_structures:
				self.status.append("secondary structure found in MAS")

			else:
				self.valid = True

				Path(f'./{self.output_directory}/Folds').mkdir(parents=True, exist_ok=True)
				fold = foldClass(self.name, self.seq, self.alignment_file, self.locus, self.strand, self.mas, self.output_directory)

				self.assess_miRNA()



	def __str__(self):

		hp_len = len(self.seq)

		def read_string(pos, read, depth):
			s = "." * pos + read + "." *(hp_len - len(read) - pos) + " a=" + str(depth)
			return(s)

		def mismatch_to_lower(pos, read):
			out = ''

			for i,r in enumerate(read):

				if self.seq[i+pos] != r:
					out += r.lower()
				else:
					out += r
			return(out)


		out = []
		out.append("\nPreliminary tests:")
		out.append("stranded:", self.stranded)
		out.append("short_enough:", self.short_enough)



		out.append("\nHairpin sequence:")
		out.append(self.seq)
		out.append(self.fold)
		out.append(read_string(self.seq.index(self.mas), self.mas, self.mas_c[self.mas]))

		if self.star_found:
			out.append(read_string(self.seq.index(self.star), self.star, self.mas_c[self.star]))

		out.append('')

		for i in range(hp_len):

			try:
				reads = self.pos_d[i]
			except KeyError:
				reads = []

			for read in reads:
				print_read = mismatch_to_lower(i, read)
				out.append(read_string(i, print_read, self.mas_c[read]))


		out.append("\nDuplex sequence:")

		out.append(self.duplex_mas)
		out.append(self.duplex_fold)
		out.append(self.duplex_star)


		out.append("\nRuling:")
		out.append(self.ruling)


		out.append("")
		for key, val in self.ruling_d.items():
			out.append(f"{key} : {val}")


		out.append("\nStatus:")
		out += self.status



		return("\n".join(map(str,out)))


	# def strucVis(self):


	# 	Path(f'./{self.output_directory}Folds').mkdir(parents=True, exist_ok=True)

	# 	bam_file = f"./{self.output_directory}Folds/{self.name}.bam"


	# 	with open(bam_file, 'wb') as outf:
	# 		sam_call = ['samtools', 'view', '-b', '-h', self.alignment_file, self.locus]
	# 		p = Popen(sam_call, stdout=outf, stderr=PIPE)
	# 		p.wait()


	# 	index_call = ['samtools', 'index', bam_file]
	# 	p = Popen(index_call, stdout=PIPE, stderr=PIPE)
	# 	p.wait()


	# 	output_file = f"./{self.output_directory}Folds/{self.name}.ps"

	# 	if self.strand == "+":
	# 		strand = 'plus'
	# 	elif self.strand == '-':
	# 		strand = 'minus'
	# 	else:
	# 		sys.exit("unrecognized strand")

	# 	call = ['strucVis', 
	# 	'-b', bam_file, 
	# 	'-g', self.genome_file, 
	# 	'-c', self.locus, 
	# 	'-s', strand,
	# 	'-p', output_file,
	# 	'-n', self.name
	# 	]

	# 	p = Popen(call, stdout=PIPE, stderr=PIPE)
	# 	out, err = p.communicate()

	# 	# print(out, err)


	# 	os.rename(output_file, output_file.replace(".ps", ".eps"))
	# 	os.remove(bam_file)
	# 	os.remove(bam_file + ".bai")


	def get_locus(self, locus, strand, input_mas):

		locus, strand, input_mas = self.locus, self.strand, self.input_mas

		chrom, start, stop = self.chrom, self.start, self.stop

		genome_file, alignment_file = self.genome_file, self.alignment_file

		seq = samtools_faidx(locus, strand, genome_file)
		fold, mfe, pairing = RNAfold(seq)

		mas_d = {}
		pos_d = {}
		input_mas_coords = False
		# print(seq)
		# print(fold)

		for read in samtools_view(alignment_file, locus=locus):

			sam_strand, sam_length, _, sam_pos, sam_chrom, sam_rg, sam_read = read

			# if ignore_replication:
				# sam_rg = 'all'
			sam_rg = 'all'

			# sam_read = sam_read.replace("T",'U')

			if sam_strand == "-":
				sam_read = complement(sam_read[::-1])

			if sam_pos >= start and sam_pos + sam_length <= stop:
				if sam_strand == strand:

					# if sam_read == input_mas:

					# 	sys.exit()

					if not input_mas_coords and sam_read == input_mas:
						input_mas_coords = f"{sam_chrom}:{sam_pos}-{sam_pos + sam_length+1}"

					if strand == '+':
						corrected_pos = sam_pos - start
					else:
						corrected_pos = stop - sam_pos - sam_length + 1

					try:
						pos_d[corrected_pos].add(sam_read)
					except KeyError:
						pos_d[corrected_pos] = set([sam_read])

					try:
						mas_d[sam_rg].update([sam_read])
					except KeyError:
						mas_d[sam_rg] = Counter()
						mas_d[sam_rg].update([sam_read])

		return(seq, fold, mfe, pairing, mas_d, pos_d, input_mas_coords)

	def find_secondary_structures(self, fold):
		# print(fold)

		if "(" in fold and ")" in fold:
			return(True)
		else:
			return(False)
		# sys.exit()


	def find_star(self, offset=2):

		# print('mas')
		# print(" "* offset + "".join([self.seq[p] for p in self.mas_positions]))
		# print(" "* offset + "".join([self.fold[p] for p in self.mas_positions]))

		for left_off, left_pos in enumerate(self.mas_positions):
			if self.pairing[left_pos] != ".":
				break

		for right_off, right_pos in enumerate(self.mas_positions[::-1]):
			if self.pairing[right_pos] != ".":
				break


		# print(left_off, left_pos, right_off, right_pos)
		# print(self.pairing[right_pos])
		# print(self.pairing[left_pos])
		star_right_pos = self.pairing[left_pos] + left_off + offset
		star_left_pos  = self.pairing[right_pos] - right_off + offset


		self.star_positions = [r for r in range(star_left_pos, star_right_pos+1)]


		if self.star_positions == []:
			self.status.append("no star positions found")
			return False

		if max(self.star_positions) >= len(self.seq) or min(self.star_positions) < 0:
			self.status.append("star expands outside of hairpin")
			return False

		if len(set(self.mas_positions).intersection(self.star_positions)) > 0:
			self.status.append("mas and proposed star overlap")
			return False




		# print(self.star_positions)
		star = "".join([self.seq[p] for p in self.star_positions])
		star_fold = "".join([self.fold[p] for p in self.star_positions])
		# print(star_fold[::-1])
		# print(star[::-1])

		# print(self.mas_positions)
		# print(self.star_positions)

		# print(star_left_pos, star_right_pos)


		m_seq  = deque([self.seq[p] for p in self.mas_positions])
		m_fold = deque([self.fold[p] for p in self.mas_positions])
		s_seq  = deque([self.seq[p] for p in self.star_positions[::-1]])
		s_fold = deque([self.fold[p] for p in self.star_positions[::-1]])

		m_out = ' ' * offset
		s_out = s_seq.popleft() + s_seq.popleft()
		f_out = ' ' * offset

		s_fold.popleft()
		s_fold.popleft()

		while len(m_seq) > 0:

			m_f = m_fold.popleft()
			m_s = m_seq.popleft()


			try:
				s_f = s_fold.popleft()
				s_s = s_seq.popleft()
			except IndexError:
				s_f = " "
				s_s = " "

			if s_s == ' ':
				f_out += " "
				s_out += s_s
				m_out += m_s

			elif m_f == "." and s_f == ".":
				f_out += "."
				s_out += s_s
				m_out += m_s

			elif m_f == ".":
				f_out += "."
				s_out += "-"
				m_out += m_s

				s_seq.appendleft(s_s)
				s_fold.appendleft(s_f)


			elif s_f == ".":
				f_out += "."
				s_out += s_s
				m_out += "-"

				m_seq.appendleft(m_s)
				m_fold.appendleft(m_f)


			else:
				f_out += "|"
				s_out += s_s
				m_out += m_s

			# input()



		# print(m_out)
		# print(f_out)
		# print(s_out)

		self.star = star

		self.duplex_mas  = m_out
		self.duplex_fold = f_out
		self.duplex_star = s_out


		return(True)
		# dup_mas = self.mas[:offset]
		# dup_mas_fold = self.
		# dup_star = " " * offset

		# while True:

	def assess_miRNA(self):


		locus_depth = sum(self.mas_c.values())

		def test_mfe():
			# <0.2 kcal/mol/nucleotide

			mfe_per_nt = self.mfe / (self.stop - self.start)
			self.ruling_d['mfe_per_nt'] = mfe_per_nt

			if mfe_per_nt < -0.2:
				return("x")

			return("-")


		def test_duplex_mismatch():

			total_mismatches = 0
			asymetric_mismatches = 0

			m_length = 0
			s_length = 0

			for i in range(len(self.duplex_mas)):

				m = self.duplex_mas[i]
				f = self.duplex_fold[i]
				s = self.duplex_star[i]

				if f == ".":
					if m != "-":
						m_length += 1

					if s != "-":
						s_length += 1

				else:

					if m_length + s_length > 0:

						total_mismatches += max([m_length, s_length])
						asymetric_mismatches += abs(m_length - s_length)



					m_length = 0
					s_length = 0


			out = ''

			self.ruling_d['mismatches_total'] = total_mismatches
			self.ruling_d['mismatches_asymm'] = asymetric_mismatches

			if total_mismatches <= 5:
				out += "x"
			else:
				out += '-'

			if asymetric_mismatches <= 3:
				out += 'x'
			else:
				out += '-'

			return(out)


		def test_secondary_structure():
			# u = duplex_mas_set.intersection(duplex_star_set)
			# print(u)
			out = ''

			out += "-" if self.mas_structures else "x"
			out += "-" if self.star_structures else "x"

			self.ruling_d['no_mas_structures'] =  not self.mas_structures
			self.ruling_d['no_star_structures'] =  not self.star_structures

			return(out)

		def test_precision():

			single_variants = 0
			for key, val in self.mas_c.items():

				# if strand == "-":
				# 	key = complement(key[::-1])
				# print(key, val)
				# print(mas, distance(key, mas))
				# print()
				if distance(key, self.mas) <= 1 or distance(key, self.star) <= 1:
					single_variants += val

			# print(single_variants)

			precision = single_variants / locus_depth 
			self.ruling_d['precision'] =  precision

			if precision > 0.75:
				return("x")
			else:
				return("-")

		def test_star_found():


			self.ruling_d['star_found'] = False


			for key, val in self.mas_c.items():
				if distance(key, self.star) <= 1:

					self.ruling_d['star_found'] = True
					return('x')

			# if self.mas_c[self.star] > 0:
			# 	return('x')
			else:
				return("-")



		test_str = ''
		test_str += test_mfe()
		test_str += " " + test_duplex_mismatch()
		test_str += " " + test_secondary_structure()
		test_str += " " + test_precision()
		test_str += " " + test_star_found()


		# print(test_str)
		# return(test_str)
		self.ruling = test_str





	def table(self):

		line = [self.name, self.locus, self.strand]
		line += [self.stranded, self.length, self.short_enough]
		line += [self.seq, self.fold, self.mfe, self.mas, self.star] 
		line += [self.duplex_mas, self.duplex_fold, self.duplex_star]
		line += [self.valid]

		line += [self.ruling] + list(self.ruling_d.values())

		return("\t".join(map(str,line)))

	# sys.exit()
	# def check_fold(start, stop):



@click.command()
@click.option("-a", "--alignment_file", 
	required=True, 
	type=click.Path(exists=True),
	help='Alignment file input (bam or cram).')

# @click.option('-r', '--annotation_readgroups', 
# 	required=False,
# 	default = False,
# 	multiple=True,
# 	help="List of read groups (RGs, libraries) to be considered for the annotation. 'ALL' uses all readgroups for annotation, but often pertainent RGs will need to be specified individually.")

@click.option("-o", "--output_directory", 
	default=f"Annotation_{round(time())}", 
	type=click.Path(),
	help="Directory name for annotation output")

@click.option('-i', "--ignore_replication",
	is_flag=True,
	help='Evaluate all readgroups together, ignoring if a miRNA is replicated')

@click.option("-m", "--max_length",
	default=300,
	help='Maximum hairpin size (default 300). Longer loci will not be considered for miRNA analysis.')

def hairpin(alignment_file, output_directory, ignore_replication, max_length):


	def get_genome_file():
		call = ['samtools', 'view', '-H', alignment_file]

		p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')

		for line in p.stdout:
			if line.startswith("@SQ"):
				# print(line)
				break
		p.wait()
		line = line.strip().split()[4]
		# print(line)
		genome = line.lstrip("UR:")

		return(genome)



	def trim_hairpin(hpc, offset=2, wiggle = 5):
		chrom, start, stop, strand = hpc.chrom, hpc.start, hpc.stop, hpc.strand

		d2d = hpc.mas_positions + hpc.star_positions
		# print(d2d)
		
		left  = min(d2d) - offset - wiggle
		right = max(d2d)          + wiggle

		# print(left, right)
		if strand == "-":
			left, right = stop-right-1, stop-left-1

		elif strand == "+":
			left, right = start+left, start+right

		else:
			sys.exit("ONLY STRANDED EXPECTED")


		# print()

		# print(chrom, start, stop, strand)

		trimmed_locus = f"{chrom}:{left}-{right}"

		# print(trimmed_locus)
		return(trimmed_locus)




	genome_file = get_genome_file()

	results_file = f"{output_directory}/Results.txt"
	assert isfile(results_file), f"results_file {results_file} not found... (Have you run annotation with this directory?)"

	input_mas_d = {}
	tops_file = f"{output_directory}/TopReads.txt"
	with open(tops_file, 'r') as f:
		header = f.readline()
		for line in f:
			line = line.strip().split('\t')

			name = line[0]
			mas  = line[1].upper().replace("T","U")

			if name not in input_mas_d.keys():
				input_mas_d[name] = mas
				# input_mas_d[line[0]] = line[1]



	header_line = "name\tlocus\tstrand\tstranded\tlength\tshort_enough\tseq\tfold\tmfe\tmas\tstar\tduplex_mas\tduplex_fold\tduplex_star\tvalid_fold\truling\tmfe_per_nt\tmismatches_asymm\tmismatches_total\tno_mas_structures\tno_star_structures\tprecision\tstar_found"

	hairpin_file = f"{output_directory}/Hairpins.txt"
	with open(hairpin_file, 'w') as outf:
		print(header_line, file=outf)



	with open(results_file, 'r') as f:
		header = f.readline()
		for line in f:
			line = line.strip().split('\t')
			name, locus = line[:2]
			strand = line[8]
			length = int(line[2])


			locus = locus.replace("..", "-")

			chrom = locus.split(":")[0]
			start = int(locus.split(":")[1].split("-")[0])
			stop  = int(locus.split(":")[1].split("-")[1])




			# print(seq, fold, mfe, sep='\n')

			cluster_selected = True
			# cluster_selected = name == 'Cl_125'

			# status = []

			stranded = strand in ["-", "+"]
			# if stranded:
			# 	status.append(f"{strand} stranded")
			# else:
			# 	status.append("not stranded")

			short_enough = length <= max_length
			# if short_enough:
			# 	status.append(f"length {length} <= {max_length}")
			# else:
			# 	status.append(f"too long {length}")





			input_mas = input_mas_d[name]

			hpc = hairpinClass(stranded, short_enough, name,locus, strand, input_mas, genome_file, alignment_file, output_directory)
			hpc.table()

			if hpc.valid:
				print()
				print(f"{hpc.ruling}\t\033[1m{name}\033[0m", length, sep='\t')
				# hpc.strucVis()

			with open(hairpin_file, 'a') as outf:
				print(hpc.table(), file=outf)

			

			if hpc.valid:
				trimmed_locus = trim_hairpin(hpc)


				trimmed_hpc = hairpinClass(stranded, short_enough, name+"-t", trimmed_locus, strand, input_mas, genome_file, alignment_file, output_directory)
				if trimmed_hpc.valid:
					print(f"{trimmed_hpc.ruling}\t\033[1m{name}\033[0m", len(trimmed_hpc.seq), 'trimmed', sep='\t')

					with open(hairpin_file, 'a') as outf:
						print(trimmed_hpc.table(), file=outf)
				# print()



			# input_mas_coords = hpc.input_mas_coords



			# imas_start = int(input_mas_coords.split(":")[-1].split("-")[0])
			# imas_stop  = int(input_mas_coords.split(":")[-1].split("-")[1])
			# # print()

			# gemini_size = 120
			# gemini_offset = 20

			# gemini = [
			# f"{chrom}:{imas_start - gemini_size + gemini_offset}-{imas_stop + gemini_offset}",
			# f"{chrom}:{imas_start - gemini_offset}-{imas_stop + gemini_size - gemini_offset}"
			# ]

			# twins = ['castor', 'pollux']

			# for gem_i, gem in enumerate(gemini):

			# 	twin = twins[gem_i]

			# 	gem_hpc = hairpinClass(gem, strand, input_mas, genome_file, alignment_file)
			# 	if gem_hpc.valid:
			# 		print(f"{gem_hpc.ruling}\t\033[1m{name}\033[0m", len(gem_hpc.seq), twin, sep='\t')







cli.add_command(precheck)
cli.add_command(annotate)
cli.add_command(hairpin)

if __name__ == '__main__':
	cli()








