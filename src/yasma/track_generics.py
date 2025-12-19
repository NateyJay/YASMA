
import pyBigWig


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

	# def process_bam(self, bam_file, aligned_depth):

	# 	import pysam

	# 	print('processing .bam to form .bw')
	# 	rpm_factor = 1000000 / aligned_depth

	# 	bamf = pysam.AlignmentFile(str(bam_file),'rb')

	# 	for c,l in self.chromosomes:
	# 		print(" ", c, end = " ")
	# 		depths = bamf.count_coverage(c, quality_threshold=0)


	# 		for i in range(l):

	# 			rpm = round(sum([depths[r][i] for r in range(4)]) * rpm_factor,3)
	# 			self.add(c, i, rpm)

	# 			if i % 1000000 == 0:
	# 				print(".", end = '')

	# 		print()



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


	def add(self, pos, length, val=1):
		'''add depth for positons along an alignment. Val may be used to select normalized values (default is 1 read) '''

		for r in range(pos, pos+length+1):
			try:
				self.depths[r] += val
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
				if self.total_reads == None:
					values.append(float(last_depth) * self.strand)
				else:
					values.append(round(last_depth / self.total_reads * 1000000, 8) * self.strand)

				# print(self.name, chrom, last_pos, "->", pos, depth, sep='\t')

				last_depth = depth
				last_pos   = pos
				span       = 1

		if last_pos < pos:
			ends.append(pos)
			starts.append(last_pos)

			if self.total_reads == None:
				values.append(float(last_depth) * self.strand)
			else:
				values.append(round(last_depth / self.total_reads * 1000000, 8) * self.strand)


		if starts[0] == ends[0]:
			starts.remove(0)
			ends.remove(0)
			values.remove(0)

		# print(values)

		self.bw.addEntries([chrom] * len(values), starts, ends=ends, values=values)

	def close(self):
		self.bw.close()

