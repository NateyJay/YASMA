# Hairpin detection module


from .generics import *

from time import time, sleep
from math import log10, sqrt
from Levenshtein import distance
from collections import deque






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
	def __init__(self, name, seq, alignment_file, locus, strand, mas, output_directory, hairpin_dir):

		self.name    = name
		self.seq     = seq
		self.locus   = locus
		self.strand  = strand

		self.hairpin_dir = hairpin_dir

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


		self.fold_file = Path(self.output_directory, self.hairpin_dir, "folds", f"{self.name}_unannotated.eps")

		# print(self.fold_file)
		Path(f"{temp_name}_ss.ps").rename(self.fold_file)

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

		outf = open(Path(self.output_directory, self.hairpin_dir, "folds", f"{self.name}.eps"), 'w')

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
		self.fold_file.unlink()
		






class hairpinClass():
	def __init__(self, params, inputs, #stranded, short_enough, 
		name, locus, strand, length, input_mas, full_pos_d):#, genome_file, alignment_file, output_directory, hairpin_dir):


		self.valid   = False
		self.status  = []

		stranded = strand in ["-", "+"]

		self.stranded         = stranded

		self.name             = name
		self.locus            = locus
		self.strand           = strand
		self.input_mas        = input_mas
		self.genome_file      = inputs['genome_file']
		self.alignment_file   = inputs['alignment_file']
		self.output_directory = params['output_directory']
		self.hairpin_dir      = params['hairpin_dir']

		self.chrom = self.locus.split(":")[0]
		self.contig = self.chrom
		self.start = int(self.locus.split(":")[-1].split("-")[0])
		self.stop  = int(self.locus.split(":")[-1].split("-")[1])

		self.length = self.stop - self.start

		self.seq = '-'
		self.fold = '-'
		self.mfe = '-'
		self.pairing = '-'
		self.pos_d = '-'
		self.input_mas_coords = '-'


		self.mas = input_mas
		self.star = '-'
		self.duplex_mas = '-'
		self.duplex_fold = '-'
		self.duplex_star = '-'


		self.ruling_d = {
		'mfe_per_nt' : '.',
		'mismatches_total' : '.',
		'mismatches_asymm' : '.',
		'no_mas_structures' : '.',
		'no_star_structures' : '.',
		'precision' : '.',
		'star_found' : '.'
		}

		self.ruling = '- . .. .. . .'



		if not stranded:
			self.status.append("hairpin not stranded")
			return




		self.seq, self.fold, self.mfe, self.pairing, self.read_c = self.get_locus(locus, strand, input_mas, full_pos_d)


		if sum(self.read_c.values()) == 0:
			self.status.append("locus depth is 0 reads")
			return


		# print(self.seq, self.fold, self.mfe, self.pairing)
		# print(self.pos_d)
		# input()

		# try:
		# 	mas_d['all']
		# except KeyError:
		# 	return

		# self.read_c = mas_d['all']
		# self.mas = self.read_c.most_common(1)[0][0]



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
			self.ruling = self.ruling[:7] + "-" + self.ruling[8:]
			return

		self.star_found = self.find_star()

		if self.star_found:

			# print(self.star)
			# print(" " * self.seq.index(self.star) + self.star)



			self.star_structures = self.find_secondary_structures("".join([self.fold[p] for p in self.star_positions]))

			if self.star_structures:
				self.status.append("secondary structure found in STAR")
				self.ruling = self.ruling[:8] + "-" + self.ruling[9:]


			else:
				self.valid = True

				# Path(self.output_directory, self.hairpin_dir, 'folds').mkdir(parents=True, exist_ok=True)
				fold = foldClass(self.name, self.seq, self.alignment_file, self.locus, self.strand, self.mas, self.output_directory, self.hairpin_dir)

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
		out.append(f"stranded: {self.stranded}")
		out.append(f"short_enough: {self.short_enough}")



		out.append("\nHairpin sequence:")
		out.append(self.seq)
		out.append(self.fold)
		out.append(read_string(self.seq.index(self.mas), self.mas, self.read_c[self.mas]))

		if self.star_found:
			out.append(read_string(self.seq.index(self.star), self.star, self.read_c[self.star]))

		out.append('')

		for i in range(hp_len):

			try:
				reads = self.pos_d[i]
			except KeyError:
				reads = []

			for read in reads:
				print_read = mismatch_to_lower(i, read)
				out.append(read_string(i, print_read, self.read_c[read]))


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


	def get_locus(self, locus, strand, input_mas, pos_d):

		locus, strand, input_mas    = self.locus, self.strand, self.input_mas
		chrom, start, stop          = self.chrom, self.start, self.stop
		genome_file, alignment_file = self.genome_file, self.alignment_file

		seq = samtools_faidx(locus, strand, genome_file)
		fold, mfe, pairing = RNAfold(seq)

		# print(seq)
		# print(input_mas)
		# print(input_mas in seq)
		# sys.exit()


		contig = locus.split(":")[0]
		start = int(locus.split(":")[-1].split("-")[0])
		stop  = int(locus.split(":")[-1].split("-")[1])

		read_c = Counter()


		for pos in range(start, stop):

			try:
				reads = pos_d[pos]
			except KeyError:
				continue

			for read in reads:
				if self.strand == "-":
					read = complement(read[::-1])

				if len(read) + pos <= stop:
					read_c[read] += 1




		# def sub_process_positions(pos_d):
		# 	'''Processes locus for the sub or trimmed locus.

		# 	Returns
		# 	'''

		# 	out_d = dict()

		# 	for pos in pos_d.keys():

		# 		reads = pos_d[pos]

		# 		for read in reads:

		# 			# if strand == '+':
		# 			# 	corrected_pos = pos - start
		# 			# else:
		# 			# 	corrected_pos = stop - pos - len(read) + 1

		# 			read_c[read] += 1


		# 		# out_d[corrected_pos] = reads

		# 		# if pos  > stop:
		# 		# 	break

		# 	return read_c

		# read_c = sub_process_positions(pos_d)



		return(seq, fold, mfe, pairing, read_c)

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


		locus_depth = sum(self.read_c.values())

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
			for key, val in self.read_c.items():

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


			for key, val in self.read_c.items():
				if distance(key, self.star) <= 1:

					self.ruling_d['star_found'] = True
					return('x')

			# if self.read_c[self.star] > 0:
			# 	return('x')
			else:
				return("-")

		test_str = 'x'
		test_str += " " + test_mfe()
		test_str += " " + test_duplex_mismatch()
		test_str += " " + test_secondary_structure()
		test_str += " " + test_precision()
		test_str += " " + test_star_found()

		self.ruling = test_str

	def table(self):

		line = [self.name, self.locus, self.strand]
		line += [self.stranded, self.length]
		line += [self.seq, self.fold, self.mfe, self.mas, self.star] 
		line += [self.duplex_mas, self.duplex_fold, self.duplex_star]
		line += [self.valid]

		line += [self.ruling] + list(self.ruling_d.values())

		return("\t".join(map(str,line)))

	# sys.exit()
	# def check_fold(start, stop):

	def status_line(self, sizecall, mature_d):

		miRNA = ''
		if mature_d:
			for p in range(self.start, self.stop + 1):
				try:
					miRNA = mature_d[(self.contig, p, self.strand)]
					break
				except KeyError:
					pass		

		if not self.valid:
			status = self.status
		else:
			status = ''

		print(self.ruling, len(self.seq), sizecall, self.name, self.locus, miRNA, status, sep='\t', flush=True)




def read_locus(alignment_file, contig, start, stop, strand):
	pos_d = {}
	# read_c = Counter()

	for read in samtools_view(alignment_file, contig=contig, start=start, stop=stop):

		sam_strand, sam_length, _, sam_pos, sam_chrom, sam_rg, sam_read, sam_read_id = read

		if sam_strand == "-":
			sam_read = complement(sam_read[::-1])

		if sam_pos >= start and sam_pos + sam_length <= stop:
			if sam_strand == strand:

				if strand == '+':
					corrected_pos = sam_pos 
				else:
					corrected_pos = sam_pos - sam_length + 1

				try:
					pos_d[corrected_pos].append(sam_read)
				except KeyError:
					pos_d[corrected_pos] = [sam_read]


				# read_c[sam_read] += 1

	return(pos_d)#, read_c)


def bowtie_matures(ma_file, genome_file):

	out_d = dict()

	temp_file = "temp.fa"
	with open(temp_file, 'w') as outf:
		with open(ma_file, 'r') as f:
			for line in f:
				if not line.startswith(">"):
					line = line.upper()
					line = line.replace("U","T")

				outf.write(line)


	bowtie_call = ['bowtie', '-f']

	# if ".fa" in ma_file.suffixes or ".fasta" in ma_file.suffixes:
	# 	pass
	# else:
	# 	sys.exit("matures only accepts .fa/.fasta format")


	bowtie_call += ['-v', '0', '-p', '1', '-S', '-a']

	# if bowtie_version >= 1.3:
	# 	bowtie_call.append("-x")

	bowtie_call.append(str(genome_file.with_suffix('')))

	bowtie_call.append(str(temp_file))
	p = Popen(bowtie_call, encoding=ENCODING, stdout=PIPE, stderr=PIPE)

	for line in p.stdout:
		line = line.strip().split('\t')
		# print(line)

		if line[0].startswith("@"):
			continue

		elif line[1] == '4':
			continue

		name   = line[0]
		contig = line[2]

		strand = "+" if line[1] == '0' else "-"
		start  = int(line[3])
		stop   = start + len(line[9])


		for pos in range(start-5, stop+5+1):

			out_d[(contig, pos, strand)] = name





	p.wait()

	return out_d






@cli.command(group="Calculation", help_priority=3)

@click.option("-a", "--alignment_file", 
	required=False, 
	type=click.UNPROCESSED, callback=validate_path,
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

@click.option("-g", "--genome_file", 
	# default=f"Annotation_{round(time())}", 
	required=False,
	# type=click.Path(exists=True),
	type=click.UNPROCESSED, callback=validate_path,
	help='Genome or assembly which was used for the original alignment.')

@click.option('-i', "--ignore_replication",
	is_flag=True,
	help='Evaluate all readgroups together, ignoring if a miRNA is replicated')

@click.option("-m", "--max_length",
	default=300,
	help='Maximum hairpin size (default 300). Longer loci will not be considered for miRNA analysis.')


@click.option("--matures",
	type=click.Path(),
	help='location for a fasta of mature miRNAs which will be used to spot orthologs.')

@click.option("--annotation_folder")

@click.option("-n", "--name")


@click.option('--ignore_subhairpins', is_flag=True, default=False, help='This prevents folding of sub-hairpins in long stranded loci')

# @click.option("--method", 
# 	default="Poisson", 
# 	help="Annotator algorithm used (Poisson or Dicer)")

def hairpin(**params):
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

	matures              = params['matures']
	ignore_replication   = params['ignore_replication']
	max_length           = params['max_length']

	name = params['name']

	if name:
		hairpin_dir = Path(output_directory, f'hairpin_{name}')
	else:
		hairpin_dir = Path(output_directory, f'hairpin')

	hairpin_dir.mkdir(parents=True, exist_ok=True)

	params['hairpin_dir'] = hairpin_dir
	params['output_directory'] = output_directory

	if params['matures']:
		mature_d = bowtie_matures(matures, genome_file)

	else:
		mature_d = None


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


	def locus_steps(locus, length, step):
		chrom = locus.split(":")[0]
		full_start = int(locus.split(":")[1].split("-")[0])
		full_stop  = int(locus.split(":")[1].split("-")[1])

		start = full_start
		stop  = full_start + length


		while True:
			yield (f"{chrom}:{start}-{stop}", length)

			start += step
			stop  += step

			if stop > full_stop:
				return (f"{chrom}:{full_stop-length}-{full_stop}", length)


	results_file = Path(output_directory, 'tradeoff', "loci.txt")


	assert results_file.is_file(), f"results_file {results_file} not found... (Have you run annotation with this directory?)"

	input_mas_d = {}
	# tops_file = f"{output_directory}/tradeoff/reads.txt"
	tops_file = Path(output_directory, 'tradeoff', "reads.txt")
	with open(tops_file, 'r') as f:
		header = f.readline()
		for line in f:
			line = line.strip().split('\t')

			name = line[0]
			mas  = line[1].upper().replace("T","U")

			if name not in input_mas_d.keys():
				input_mas_d[name] = mas
				# input_mas_d[line[0]] = line[1]




	other_mas_d = {}

	other_mas_file = Path(output_directory, 'tradeoff', 'reads.txt')
	with open(other_mas_file, 'r') as f:
		f.readline()

		for line in f:
			line = line.strip().split("\t")

			cluster, seq = line[:2]

			try:
				other_mas_d[cluster].append(seq)
			except:
				other_mas_d[cluster] = [seq]
		




	header_line = "name\tlocus\tstrand\tstranded\tlength\tseq\tfold\tmfe\tmas\tstar\tduplex_mas\tduplex_fold\tduplex_star\tvalid_fold\truling\tmfe_per_nt\tmismatches_asymm\tmismatches_total\tno_mas_structures\tno_star_structures\tprecision\tstar_found"

	Path(output_directory, hairpin_dir, "folds").mkdir(parents=True, exist_ok=True)
	hairpin_file = Path(output_directory, hairpin_dir, "hairpins.txt")
	with open(hairpin_file, 'w') as outf:
		print(header_line, file=outf)



	print("""
stranded
┋
┋ mfe_per_nt
┋ ┋
┋ ┋ mismatches_total
┋ ┋ ┋
┋ ┋ ┋mismatches_asymm
┋ ┋ ┋┋
┋ ┋ ┋┋ mas_duplex_structure
┋ ┋ ┋┋ ┋
┋ ┋ ┋┋ ┋star_duplex_structure
┋ ┋ ┋┋ ┋┋
┋ ┋ ┋┋ ┋┋ precision
┋ ┋ ┋┋ ┋┋ ┋
┋ ┋ ┋┋ ┋┋ ┋ star_found
┋ ┋ ┋┋ ┋┋ ┋ ┋
v v vv vv v v""")



	for entry in read_loci(results_file):

		name     = entry['name']
		locus    = entry['locus']
		strand   = entry['strand']
		length   = int(entry['length'])
		sizecall = entry['sizecall']


		locus = locus.replace("..", "-")

		chrom = locus.split(":")[0]
		start = int(locus.split(":")[1].split("-")[0])
		stop  = int(locus.split(":")[1].split("-")[1])





		other_mas = other_mas_d[name]


		# print(seq, fold, mfe, sep='\n')

		cluster_selected = True




		stranded = strand in ["-", "+"]
		short_enough = length <= params['max_length']


		if not stranded:
			hpc = hairpinClass(params, ic.inputs, name, locus, strand, length, 'None', {})
			hpc.status_line(sizecall, mature_d)
			continue

		full_pos_d = read_locus(alignment_file, chrom, start, stop, strand)

		if short_enough:

			input_mas = input_mas_d[name]

			hpc = hairpinClass(params, ic.inputs, name, locus, strand, length, input_mas, full_pos_d)
			hpc.table()

			hpc.status_line(sizecall, mature_d)

			with open(hairpin_file, 'a') as outf:
				print(hpc.table(), file=outf)

			

			if hpc.valid:

				trimmed_locus = trim_hairpin(hpc)
				trimmed_hpc = hairpinClass(params, ic.inputs, name+'.trim', trimmed_locus, strand, length, input_mas, full_pos_d)

				if trimmed_hpc.valid:
					# print(f"{trimmed_hpc.ruling}\t\033[1m{name}\033[0m", len(trimmed_hpc.seq), sizecall, 'trimmed', sep='\t')
					trimmed_hpc.status_line(sizecall, mature_d)

					with open(hairpin_file, 'a') as outf:
						print(trimmed_hpc.table(), file=outf)


		elif not params['ignore_subhairpins']:
			# continue

			found_valid_hairpin = False

			for i,sub in enumerate(locus_steps(locus, length=200, step=50)):
				sub_locus, sub_length = sub
				# print(sub_locus)
				seq = samtools_faidx(sub_locus, strand, genome_file)
				# print(seq)
				seq = seq.upper()
				seq = seq.replace("T", "U")

				found = False
				for mas in other_mas:
					# print(mas)
					if mas in seq:
						found = True
						break

				# input(found)
				if found:

					hpc = hairpinClass(params, ic.inputs, f"{name}.sub{i}", sub_locus, strand, sub_length, mas, full_pos_d)
					# print(hpc.valid, hpc.ruling)
					# print(hpc.status)
					# input()
					if hpc.valid:
						trimmed_locus = trim_hairpin(hpc)
						trimmed_hpc = hairpinClass(params, ic.inputs, f"{name}.sub{i}", trimmed_locus, strand, sub_length, mas, full_pos_d)

						if trimmed_hpc.valid:
							found_valid_hairpin = True
							# print()
							# print(f"{hpc.ruling}\t\033[1m{name}\033[0m", length, sep='\t')
							# print(f"{hpc.ruling}\t\033[1m{name}\033[0m", len(trimmed_hpc.seq), sizecall, f'sub{i}', sep='\t')
							hpc.status_line(sizecall, mature_d)

							with open(hairpin_file, 'a') as outf:
								print(trimmed_hpc.table(), file=outf)


			if not found_valid_hairpin:
				# print(f"\t\t\t\t{name} <- no valid sub-hairpins")
				hpc.status_line(sizecall, mature_d)







