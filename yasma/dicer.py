## sRNA annotation module

import sys
import click
from pathlib import Path
from os.path import isfile, isdir
from time import time, sleep
from collections import Counter, deque
from itertools import count, chain
from statistics import median, mean
from pprint import pprint


from .generics import *
from .cli import cli



@cli.command(group='Annotation', help_priority=2)
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

def dicer(alignment_file, annotation_readgroups, dicercall, output_directory, force, partial_wigs, window, merge_dist, pad, rpm_cutoff, extension_ratio, dicer_ratio):
	"""Annotator based on expected dicer-derived sRNA profiles."""
	# print('run annotation')




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
	gff_file   = f"{output_directory}/Dicer.annotation.gff3"

	with open(gff_file, 'w') as outf:
		print("##gff-version 3", file=outf)

		for c, l in chromosomes:
			print(f"##sequence-region   {c} 1 {l}", file=outf)


	# count_file = f'{output_directory}/Counts.txt'
	# with open(count_file, 'w') as outf:
	# 	print("cluster", 'ann_depth', 'tot_depth', "\t".join(bam_rgs), sep='\t', file=outf)



	results_file = f"{output_directory}/Dicer.results.txt"
	with open(results_file, 'w') as outf:
		print("locus\tname\tlocus_peak\tlength\tgap\tdepth\trpm\tdepth:length\tfrac_top\tstrand\tfrac_dicer\tdcr_reads\tnon_reads\tdicercall\tfrac_dicercall\t" + "\t".join(map(str, dcr_range)), file=outf)


	reads_file = f"{output_directory}/Dicer.reads.txt"
	with open(reads_file, 'w') as outf:
		print(TOP_READS_HEADER, file=outf)





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



	# depth_c = get_rg_depth(output_directory, alignment_file)


	# library_depth = 0
	# for key in depth_c.keys():
	# 	if key[0] in annotation_readgroups:
	# 		library_depth += depth_c[key]

	# print(library_depth)


	depth_c = get_global_depth(output_directory, alignment_file, aggregate_by=['rg'])

	library_depth = sum([depth_c[rg] for rg in annotation_readgroups])
	# print(library_depth)
	# sys.exit()


	read_equivalent = 1 / library_depth * 1000000
	depth_cutoff = library_depth / rpm_cutoff / 1000000
	ext_cutoff = rpm_cutoff * extension_ratio


	print()
	print('\033[1m[Depth settings]\033[0m')
	print(f'    ann_rg_depth: {library_depth:,} reads')
	print(f'          1 read: {round(read_equivalent,5)} rpm')
	print(f"      rpm_cutoff: {rpm_cutoff} rpm -> {round(depth_cutoff,2)} reads")
	print(f"      ext_cutoff: {ext_cutoff} rpm -> {round(depth_cutoff*extension_ratio,2)} reads")
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
	sam_strand, sam_length, sam_size, sam_pos, sam_chrom, sam_rg, sam_seq, sam_read = read
	# strand,    length,     size,     sam_pos,  sam_chrom, rg,    seq, read_id)


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
			depth_c  = Counter()

			read_starts = []
			read_stops  = []


			for r in reads:
				sam_pos, sam_length, sam_size, sam_strand, sam_rg, sam_seq = r



				for r in range(sam_length):
					depth_c[r] += 1


				read_starts.append(sam_pos)
				read_stops.append(sam_pos + sam_length)

				len_c.update([sam_length])
				size_c.update([sam_size])
				strand_c.update([sam_strand])
				rg_c.update([sam_rg])
				if sam_strand == "-":
					sam_seq = complement(sam_seq[::-1])
				read_c.update([sam_seq])




			peak_positions = [depth_c.most_common()[0][0]]
			deepest = depth_c.most_common()[0][1]
			for i,d in depth_c.most_common():
				if d < deepest:
					peak_positions.append(last_i)
					break

				last_i = i




			# print(max(depth_c))
			# pprint(depth_c)
			# print(peak_positions)
			# sys.exit()
			# pprint(depth_buffer)

			# depth_buffer = [d for d in samtools_depth(alignment_file, annotation_readgroups, coords)]
			# peak_positions = [i for r in if d == max(depth_buffer)]
			# print(peak_positions)
			locus_peak = f'{chrom}:{start + peak_positions[0] + pad}-{start + peak_positions[-1] + pad}'
			# locus_peak = 'not calculated'
			# print(locus_peak)

			# sys.exit()



			# print(start, stop)
			start = min(read_starts) - pad
			stop  = max(read_stops) + pad


			name = f"Cl_{next(cluster_counter)}"


			if len(reads) == 0:
				print(f"WARNING: {name} detected no reads. Likely an error. Skipping.")
				return(0,0)



			dist_to_last = start - self.last_end
			self.last_end = stop


			if start < 0:
				start = 0

			coords = f"{chrom}:{start}-{stop}"
			length = stop - start
			n_reads = len(reads)
			frac_top = round(strand_c["+"] / n_reads,3)
			frac_dicercall = round(size_c['dcr'] / n_reads, 3)
			rpm = round(n_reads / library_depth * 1000000, 2)


			top_reads_save(read_c, reads_file, read_equivalent, name)




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




			to_print = [coords, name, locus_peak, length, dist_to_last, n_reads, rpm, depth_by_length, frac_top, strand]
			to_print += [frac_dicercall, size_c['dcr'],  size_c['non']]
			to_print += [predominant_length, predominant_length_depth]
			to_print += [len_c[d] for d in dcr_range]
			to_print = "\t".join(map(str, to_print))




			with open(results_file, 'a') as outf:
				print(to_print, sep='\t', file=outf)



			with open(gff_file, 'a') as outf:
				if start == 0:
					start += 1
				print(f"{chrom}\tsmoothLocus\tnc_RNA\t{start}\t{stop}\t.\t{strand}\t.\tID={name};dicercall={predominant_length};frac_dicercall={predominant_length_depth}", file=outf)


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
		locus = locusClass() #merge_dist, chrom, pad, cluster_counter, library_depth, read_equivalent, dcr_range, annotation_readgroups, bam_rgs)



		while sam_chrom == chrom:



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



					locus.add((sam_pos, sam_length, sam_size, sam_strand, sam_rg, sam_seq))

				try:
					read = next(sam_iter)
					sam_strand, sam_length, sam_size, sam_pos, sam_chrom, sam_rg, sam_seq, sam_read = read


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

		if locus_count > 0:
			med_length = median(locus_lengths)
			med_gap = median(locus_gaps)
		else:
			med_length = -1
			med_gap = -1

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










