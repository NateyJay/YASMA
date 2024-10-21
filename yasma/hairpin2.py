# Hairpin detection module

import os
import click
from subprocess import PIPE, Popen, call
from pathlib import Path
from os.path import isfile, isdir
from time import time, sleep
from math import log10, sqrt
import pysam
from collections import Counter


from .generics import *
from .cli import cli
from Levenshtein import distance




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

		# self.find_5p_angle()

		self.write()

		# sys.exit()

		
	def get_depth(self):

		if self.strand == "+":
			flag = "0"
		elif self.strand == "-":
			flag = "16"
		else:
			sys.exit("what???")

		p1 = Popen(['samtools', 'view', "-h", self.alignment_file, self.locus], 
			stdout=PIPE, encoding=ENCODING)

		p2 = Popen(['samtools', 'depth', '-'], 
			stdin=PIPE, stdout=PIPE, encoding=ENCODING)
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
					encoding=ENCODING)
		out, err = p.communicate(f">{temp_name}\n{self.seq}")


		self.fold_file = f"{self.output_directory}/hairpin/folds/{self.name}_unannotated.eps"

		# print(self.fold_file)
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


	def find_5p_angle(self):
		x, y = self.coor[0]

		rotation = [
		[ 15, 0  ],
		[ 10, 5  ],
		[  5, 10 ],
		[  0, 15 ],
		[ -5, 10 ],
		[-10, 5  ],
		[-15, 0  ],
		[-10, -5 ],
		[ -5, -10],
		[  0, -15],
		[  5, -10],
		[ 10, -5 ]
		]



		for r in rotation:
			rx = x + r[0]
			ry = y + r[1]

			for cx, cy in self.coor:

				print(15.0^2)
				print((rx - cx)^2)

				print((rx - cx)^2 + (ry - cy)^2)

				err = sqrt((rx - cx)^2 + (ry - cy)^2)
				print(err)
				sys.exit()




	def write(self):

		outf = open(f"{self.output_directory}/hairpin/folds/{self.name}.eps", 'w')

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
0.4 0.4 1 setrgbcolor
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

					print(f"{i + 1} {r} {g} {b} maplemark", file=outf)


				print(line, file=outf)

				print("1 show5", file=outf)


				for i in range(self.seq.index(self.mas),self.seq.index(self.mas) + len(self.mas)):
					print(f"{i + 1} borderdraw", file=outf)

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

		try:
			mas_d['all']
		except KeyError:
			return

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


		self.mas_positions = [r + self.seq.index(self.mas) for r in range(len(self.mas))]

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

				Path(f'./{self.output_directory}/hairpin/folds').mkdir(parents=True, exist_ok=True)
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


		contig = locus.split(":")[0]
		start = int(locus.split(":")[-1].split("-")[0])
		stop  = int(locus.split(":")[-1].split("-")[1])

		for read in samtools_view(alignment_file, contig=contig, start=start, stop=stop):

			sam_strand, sam_length, _, sam_pos, sam_chrom, sam_rg, sam_read, sam_read_id = read
			# strand, length, size, sam_pos, sam_chrom, rg, seq, read_id

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
		# print(self.pairing)
		# print(self.pairing[left_pos])
		# print(left_off)
		# print(offset)
		# # sys.exit()
		# print()

		try:
			star_right_pos = self.pairing[left_pos] + left_off + offset
			star_left_pos  = self.pairing[right_pos] - right_off + offset
		except TypeError:
			self.status.append("star positioning error")
			return False


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
				f_out += ":"
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



@cli.command(group="Calculation", help_priority=3, name='hairpin2')

@click.option("-a", "--alignment_file", 
	required=False, 
	type=click.Path(exists=True),
	help='Alignment file input (bam or cram).')

# @click.option('-r', '--annotation_readgroups', 
# 	required=False,
# 	default = False,
# 	multiple=True,
# 	help="List of read groups (RGs, libraries) to be considered for the annotation. 'ALL' uses all readgroups for annotation, but often pertainent RGs will need to be specified individually.")

@click.option("-o", "--output_directory",
	required=True, 
	type=click.Path(),
	help="Directory name for annotation output")

@click.option('-i', "--ignore_replication",
	is_flag=True,
	help='Evaluate all readgroups together, ignoring if a miRNA is replicated')

@click.option("-m", "--max_length",
	default=300,
	help='Maximum hairpin size (default 300). Longer loci will not be considered for miRNA analysis.')

# @click.option("--method", 
# 	default="Poisson", 
# 	help="Annotator algorithm used (Poisson or Dicer)")

def hairpin2(**params):
	"""Evaluates annotated loci for hairpin or miRNA structures."""

	

	rc = requirementClass()
	rc.add_samtools()
	rc.add_RNAfold()
	rc.check()

	ic = inputClass(params)
	ic.check(['alignment_file'])

	output_directory     = str(ic.output_directory)
	alignment_file       = ic.inputs['alignment_file']
	genome_file          = ic.inputs['genome_file']

	ignore_replication   = params['ignore_replication']
	max_length           = params['max_length']







	annotation_file = Path(output_directory, 'tradeoff', 'loci.txt')

	genf = pysam.FastaFile(genome_file)


	def reverse_complement(seq):
		comp = {"A":"T","T":"A","G":"C","C":"G","N":"N"}

		seq = [comp[s.upper()] for s in seq[::-1]]

		seq = "".join(seq)
		return seq


	def fold_from_coords(contig, start, stop, strand):

		seq = genf.fetch(reference=contig, start=start, end=stop)

		if strand == '-':
			seq = reverse_complement(seq)


		call = ['RNAfold', '--noPS']

		p = Popen(call,
					  stdout=PIPE,
					stderr=PIPE,
					stdin=PIPE,
					encoding=ENCODING)
		out, err = p.communicate(f">.\n{seq}")

		out = out.strip().split()

		return(seq, out[2])



	def get_seqs(contig, start, stop):

		seqs = Counter()

		bamf = pysam.AlignmentFile(alignment_file)

		for read in bamf.fetch(contig=contig, start=start, stop=stop):

			if read.is_forward:
				strand = "+"
			else:
				strand = "-"


			key = (read.reference_start, read.query_length, strand)
			
			seqs[key] += 1


		bamf.close()

		return(seqs)



	def get_folds(start, stop):
		fold_sizes = [100, 200, 300]

		out = []

		outer_start = start - max(fold_sizes)
		outer_stop  = stop  + max(fold_sizes)

		if outer_start < 0:
			outer_start = 0


		for size in fold_sizes:

			increment = round(size / 20)

			f_start = outer_start
			while True:
				f_stop = f_start + size

				if f_stop > outer_stop:
					out.append((size,outer_stop - size, outer_stop))
					break

				out.append((size, f_start, f_stop))

				f_start += increment

		return out


	def get_unpaired_positions(fold):

		out = []
		for i,f in enumerate(fold):
			if f == ".":
				out.append(i)

		return(out)






	def get_pins(fold):
		# print(fold)

		inflections = []

		last_f = ")"
		for i,f in enumerate(fold):

			if f == ".":
				continue

			if f != last_f:
				if f == "(":
					inflections.append(i)
				# else:
				# 	inflections.append(last_i)

			last_f = f
			last_i = i


		pin_buffer = []
		pins       = [[]]
		for i,f in enumerate(fold):

			if f == ".":
				continue

			if f == "(":
				pin_buffer.append(i)

			else:
				vals = [pin_buffer.pop(), i]
				pins[-1] += vals

				for p in vals:
					if p in inflections:
						pins.append([])

		pins = [p for p in pins if len(p) > 0]
	
		for i, pin in enumerate(pins):
			pin.sort()

			unpaired_positions = []
			for r in range(pin[0], pin[-1]):
				if r not in pin:
					unpaired_positions.append(r)

			pin = [pin[0], pin[-1]]

			length = pin[1]-pin[0]
			pin.append(length)

			score = Counter(fold[pin[0]: pin[1]+1])
			pin.append(score['.'])
			pin.append(score['('] + score[')'])

			score = round(score['.']/score['('], 4)
			pin.append(score)

			pin.append(unpaired_positions)

			pins[i] = pin


		# for i,f in enumerate(fold):
		# 	to_print = [i, f]

		# 	for pi,p in enumerate(pins):
		# 		if p[0] <= i <= p[1]:
		# 			to_print.append(pi)
		# 		else:
		# 			to_print.append("")	

		# 	print("\t".join(map(str, to_print)))

		return(pins)







	with open(annotation_file, 'r') as f:
		header = f.readline()
		header = header.strip().split('\t')
		# print(header)

		for line in f:
			line = line.strip().split('\t')

			strand       = line[header.index("Strand")]
			region       = line[header.index("Locus")]
			contig       = region.split(":")[0]
			start, stop  = [int(r) for r in region.split(":")[1].split("-")]


			## MIR166A
			contig = 'NC_003071.7'
			start  = 19176108
			stop   = 19176278

			## a nearby negative control
			# start  = 19375108
			# stop   = 19375278
			strand = '+'


			seq_c = get_seqs(contig, start, stop)


			# fold = fold_from_coords(contig, start, stop, strand)
			# process_fold(fold)

			with open('pins.txt', 'w') as outf:
				print('contig', 'fold_start', 'fold_stop','fold_i', 'fold_length', 'pin_start', 'pin_stop', 'pin_i', 'pin_length', 'pin_paired','pin_unpaired','pin_score', file=outf, sep='\t')
				fold_windows = get_folds(start, stop)

				for fold_i,fold_window in enumerate(fold_windows):

					fold_length, fold_start, fold_stop = fold_window


					seq, fold = fold_from_coords(contig, fold_window[1], fold_window[2], strand)
					# fold = fold_from_coords(contig, 19375108, 19375808, strand)
					# pins = get_pins(fold)

					# for pin_i,pin in enumerate(pins):

					# 	pin_start, pin_stop, pin_length, paired, unpaired, pin_score, unpaired_positions = pin

					# 	print()
					# 	print(f"{contig} : {pin_start + start} - {pin_stop + start}")

					# print(pin)
					# unpaired_positions = [u + fold_start for u in get_unpaired_positions(fold)]

					seq_keys = [k for k in seq_c.keys() if fold_start <= k[0] <= fold_stop-k[1]]



					if len(seq_keys) == 0:
						continue

					strand_c = Counter()

					for seq_key in seq_keys:
						pos, length, strand = seq_key
						strand_c[strand] += seq_c[seq_key]

					frac_top = strand_c['+'] / sum(strand_c.values())

					if frac_top < 0.2:
						fold_strand = "-"
					elif frac_top > 0.8:
						fold_strand = "+"
					else:
						continue

					seq_keys = [k for k in seq_keys if k[2] == fold_strand]

					# print()
					# print(f"{contig}:{fold_start}-{fold_stop}")
					# print(f"window_n: {fold_i}")
					# print(f"window_length: {fold_length}")

					# print(seq_keys)

					pairing_c = Counter()
					for seq_key in seq_keys:
						pos, length, strand = seq_key
						depth = seq_c[seq_key]


						rel_start = pos - fold_start
						rel_stop  = rel_start + length

						# if rel_stop > len(fold):
						# 	continue

						unp_c = Counter(fold[rel_start:rel_stop])

						for k in unp_c.keys():
							if k == "(":
								pair_score = unp_c['('] / sum(unp_c.values())
								break

							elif k == ")":
								pair_score = unp_c[')'] / sum(unp_c.values())
								break

						pair_score = round(pair_score, 4)
						pairing_c[pair_score] += depth


					pairing_scores = [k for k in pairing_c.keys()]
					pairing_scores.sort()

					running = 0
					for score in pairing_scores:
						flag = ''

						if running <= sum(pairing_c.values())/2 <= running + pairing_c[score]:
							flag = ' <-'

							print(f"{contig}:{fold_start}-{fold_stop}", fold_i, fold_length, score, pairing_c[score], sep='\t')



					# pprint(pairing_c)
					# print(seq_keys)
					# input()

					# print(contig, fold_start, fold_stop, fold_i, fold_length, pin_start, pin_stop, pin_i, pin_length,  paired, unpaired, pin_score, file=outf, sep='\t')


				# structured_positions = [i+fold_window[1] for i,p in enumerate(fold) if p != "."]
				# c.update(structured_positions)

			# min_window = min([w[1] for w in fold_windows])
			# max_window = max([w[2] for w in fold_windows])

			# with open('struc.txt', 'w') as outf:
			# 	print("position", "structure", file=outf, sep='\t')

			# 	for p in range(min_window, max_window + 1):

			# 		print(p , c[p], sep='\t', file=outf)

			# all_positions = [i for sub in dist.keys() for i in sub]
			# with open("dist.txt", 'w') as outf:
			# 	out_line = []
			# 	for j in range(min(all_positions), max(all_positions) + 1):
			# 		out_line.append(f"j{j}")
			# 	print("\t".join(map(str, out_line)), file=outf)

			# 	for i in range(min(all_positions), max(all_positions) + 1):
			# 		out_line = [f"i{i}"]

			# 		for j in range(min(all_positions), max(all_positions) + 1):
			# 			try:
			# 				out_line.append(dist[(i,j)])
			# 			except KeyError:
			# 				out_line.append(0)

			# 		print("\t".join(map(str, out_line)), file=outf)



			sys.exit()




















	# def get_genome_file():
	# 	call = ['samtools', 'view', '-H', alignment_file]

	# 	p = Popen(call, stdout=PIPE, stderr=PIPE, encoding=ENCODING)

	# 	for line in p.stdout:
	# 		if line.startswith("@SQ"):
	# 			# print(line)
	# 			break
	# 	p.wait()
	# 	line = line.strip().split()[4]
	# 	# print(line)
	# 	genome = line.lstrip("UR:")

	# 	return(genome)



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




	# genome_file = get_genome_file()

	results_file = f"{output_directory}/tradeoff/loci.txt"
	assert isfile(results_file), f"results_file {results_file} not found... (Have you run annotation with this directory?)"

	input_mas_d = {}
	tops_file = f"{output_directory}/tradeoff/reads.txt"
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

	Path(output_directory+ "/hairpin/folds").mkdir(parents=True, exist_ok=True)
	hairpin_file = f"{output_directory}/hairpin/hairpins.txt"
	with open(hairpin_file, 'w') as outf:
		print(header_line, file=outf)



	print("""
mfe_per_nt
┋
┋ mismatches_total
┋ ┋
┋ ┋mismatches_asymm
┋ ┋┋
┋ ┋┋ no_mas_structures
┋ ┋┋ ┋
┋ ┋┋ ┋no_star_structures
┋ ┋┋ ┋┋
┋ ┋┋ ┋┋ precision
┋ ┋┋ ┋┋ ┋
┋ ┋┋ ┋┋ ┋ star_found
┋ ┋┋ ┋┋ ┋ ┋
v vv vv v v""")



	for entry in read_loci(params):

		name     = entry['name']
		locus    = entry['locus']
		strand   = entry['strand']
		length   = int(entry['length'])
		sizecall = entry['sizecall']


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
				print(f"{trimmed_hpc.ruling}\t\033[1m{name}\033[0m", len(trimmed_hpc.seq), sizecall, 'trimmed', sep='\t')

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








# if __name__ == '__main__':
# 	cli()








