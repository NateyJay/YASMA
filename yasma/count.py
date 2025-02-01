# Simple quantification of features

from .generics import *

from shutil import rmtree

# def call_count(job):

# 	name, locus = job

# 	contig = locus.split(":")[0]
# 	start = int(locus.split(":")[-1].split("-")[0])
# 	stop  = int(locus.split(":")[-1].split("-")[1])

# 	c = Counter()

# 	sam_iter = samtools_view(alignment_file, contig=contig, start=start, stop=stop)

# 	for read in sam_iter:
# 		sam_strand, sam_length, _, sam_pos, _, sam_rg, sam_seq = read


# 		if sam_pos >= start and sam_pos + sam_length <= stop:

# 			c.update([sam_rg])


# 	line = [name, locus, sum(c.values())]
# 	for rg in rgs:
# 		line.append(c[rg])

# 	# pprint(c)

# 	lock.acquire()
# 	# print(".", end='', flush=True)

# 	with open(output_file, 'a') as outf:
# 		print("\t".join(map(str, line)), file=outf)


# 	# pbar.update(1)
# 	lock.release()



# def init(l, r, a, o, ):
# 	global lock
# 	global rgs
# 	global alignment_file
# 	global output_file
# 	lock = l
# 	rgs = r
# 	alignment_file = a
# 	output_file = o

@cli.command(group='Calculation', help_priority=3)


@optgroup.group('\n  Basic options',
				help='')

@optgroup.option("-a", "--alignment_file", 
	required=False, 
	type=click.UNPROCESSED, callback=validate_path,
	help='Alignment file input (bam or cram).')

@optgroup.option("-o", "--output_directory", 
	# default=f"Annotation_{round(time())}", 
	required=False,
	type=click.UNPROCESSED, callback=validate_outdir,
	help="Directory name for annotation output. Defaults to the current directory, with this directory name as the project name.")

@optgroup.option("-c", "--conditions", 
	required=False, 
	multiple=True,
	type=click.UNPROCESSED, callback=validate_condition,
	help='Values denoting condition groups (sets of replicate libraries) for projects with multiple tissues/treatments/genotypes. Can be entered here as space sparated duplexes, with the library base_name and condition groups delimited by a colon. E.g. SRR1111111:WT SRR1111112:WT SRR1111113:mut SRR1111114:mut')

@optgroup.option("-an", "--annotation_file", 
	required=False, 
	type=click.UNPROCESSED, callback=validate_path,
	multiple=False,
	help='Locus annotation in gff, gff3, gtf, bed, or tabular format. Tabular requires contig:start-stop and locus_namein the first two columns (tab-delimited, "#" escape char).')

@optgroup.option("-n", "--name", 
	required=False,
	type=str,
	help="Name for resulting counts and analysis. Required if not using the default annotation file.")

@optgroup.option("-f", "--features", 
	required=False, 
	multiple=True,
	help='A filter to only analyze certain features of an input gff/gtf annotations (field 3). Permits all features by default. Multiple features can be included as space delimited entries. Does not apply for other annotation file types.')

@optgroup.option("--reanalyze", 
	default=False,
	is_flag=True,
	help="Flag calling for the reanalysis of annotated features. Produces a file similar to the 'loci.txt' output of tradeoff, which includes major sRNA locus dimensions.")

@optgroup.option("--include_zeros",
	is_flag=False,
	help="Include to save 0-depth rows in the deep counts. By default, these are excluded to save space (except for one entry to make sure downstream analyses will include un-found entries)")


# @optgroup.option("--ignore_unaligned",
# 	is_flag=False,
# 	help="Include to skip counting unaligned reads in deepcounts.txt. These are useful for some analyses, but it can be faster to ignore.")



def count(** params):
	"""Gets counts for all readgroups, loci, strand, and sizes."""

	rc = requirementClass()
	rc.check()

	ic = inputClass(params)
	ic.check(['alignment_file','annotation_file'])

	output_directory     = ic.output_directory
	alignment_file       = ic.inputs['alignment_file']
	conditions           = ic.inputs['conditions']
	annotation_file      = ic.inputs['annotation_file']
	project_name         = ic.inputs['project_name']

	include_zeros        = params['include_zeros']
	name                 = params['name']
	feature_filter       = params['features']

	Path(output_directory, "counts").mkdir(parents=True, exist_ok=True)
	temp_dir = Path(output_directory, 'counts', 'temp')

	try:
		rmtree(temp_dir)
	except FileNotFoundError:
		pass
	temp_dir.mkdir(parents=True, exist_ok=True)

	if annotation_file.absolute() != Path(output_directory, "tradeoff", "loci.gff3").absolute():
		if not name:

			sys.exit("Error: if annotation_file != tradeoff/loci.gff3 a name is required")


	if feature_filter:
		feature_filter = set(feature_filter)
		print("filtering annotation to include only:", feature_filter)

		if annotation_file.suffix not in ['.gff3','gff2','.gff','.gtf']:
			sys.exit(f"Error: feature filtering is not compatible with '{annotation_file.suffix}' files")


	chromosomes, libraries = get_chromosomes(alignment_file)

	if not conditions:
		conditions = {'all' : libraries}


	if name:
		name_str = name + "_"
	else:
		name_str = "tradeoff_"

	counts_file      = Path(temp_dir, f'{name_str}counts.txt')
	deep_counts_file = Path(temp_dir, f'{name_str}deepcounts.txt') 
	analysis_file    = Path(temp_dir, f'{name_str}loci.txt') 




	chrom_depth_c = get_global_depth(alignment_file, aggregate_by=['rg','chrom'])

	# keys = list(chrom_depth_c.keys())
	# for key in keys:
	# 	if key[0] in libraries:
	# 		chrom_depth_c[key[1]] += chrom_depth_c[key]

	# 	del chrom_depth_c[key]

	# for key in list(chrom_depth_c.keys()):
	# 	if key not in [c for c,l in chromosomes]:
	# 		del chrom_depth_c[key]

	aligned_read_count = sum(chrom_depth_c.values())


	print()
	print(f'counting annotation: {annotation_file}')
	print()

	c = Counter()




	feature_c = Counter()
	loci = []
	coord_d = {}

	with open(annotation_file, 'r') as f:
		if annotation_file.suffix == '.txt':
			f.readline()

		for i,line in enumerate(f):

			if not line.startswith("#"):
				line = line.strip().split("\t")

				if annotation_file.suffix == ".txt":
					coords, name = line[:2]
					coords = coords.replace("..", '-')

				elif annotation_file.suffix in ['.gff', '.gff2', '.gff3']:
					coords  = f"{line[0]}:{line[3]}-{line[4]}"
					name    = line[8].split(";")[0].split("=")[-1].strip('"')
					feature = line[2]
					feature_c[feature] += 1

					if feature_filter and feature in feature_filter:
						continue

				elif annotation_file.suffix == '.gtf':
					coords  = f"{line[0]}:{line[3]}-{line[4]}"
					name    = line[8].split(";")[0].split()[-1].strip('"')
					feature = line[2]
					feature_c[feature] += 1

					if feature_filter and feature in feature_filter:
						continue

				elif annotation_file.suffix == '.bed':
					name   = f"bed_{i}"
					coords = f"{line[0]}:{line[1]}-{line[2]}"

				else:
					print(f'file.suffix "{annotation_file.suffix}" not expected in annotation file...')
					sys.exit()


				loci.append((name, coords))
				coord_d[name] = coords

	if feature_filter:
		print("included features found:")
		found = 0
		for feature in feature_filter:
			found += feature_c[feature]
			print(f"  {feature_c[feature]}\t{feature}")

		print()
		print(f"  {sum(feature_c.values()) - found} <- feature(s) not included")



	print('')
	print('processing annotation...')


	chroms, rgs = get_chromosomes(alignment_file)

	# pprint(conditions)

	rev_conditions = {}
	for k,vs in conditions.items():
		for v in vs:
			rev_conditions[v] = k

	try:
		cond_rgs = [f"{rev_conditions[r]}.{r}" for r in rgs]
	except Exception as err:
		print(f"Unexpected {err=}, {type(err)=}")
		print("rgs:", rgs)
		print("rev_conditions:", rev_conditions)
		raise

	with open(counts_file, 'w') as outf:
		print('name','locus', "\t".join(cond_rgs), sep='\t', file=outf)

	with open(deep_counts_file, 'w') as outf:
		print('name', 'condition', 'rg','length','strand','count', sep='\t', file=outf)

	if params['reanalyze']:
		with open(analysis_file, 'w') as outf:
			print("\t".join(assessClass().header), file=outf)


	### new

	print("finding annotated positions...")
	chrom_d = {}
	for chrom, chrom_length in chroms:
		chrom_d[chrom] = [None] * chrom_length



	for i, locus in enumerate(loci):

		name, locus = locus
		# print(name, locus)

		contig = locus.split(":")[0]
		start  = int(locus.split(":")[1].split("-")[0])
		stop   = int(locus.split(":")[1].split("-")[1])

		for r in range(start, stop):
			chrom_d[contig][r] = name





	deep_c = Counter()


	bamf = pysam.AlignmentFile(alignment_file,'rb')

	if not bamf.has_index():
		print(f'   index not found for {bam}. Indexing with samtools.')

		pysam.index(str(bam))
		bamf.close()
		bamf = pysam.AlignmentFile(alignment_file,'rb')




	seq_c    = Counter()
	strand_c = Counter()
	size_c   = sizeClass()

	locus_stop   = 0
	locus_contig = ''

	read_i = 0

	in_locus = False
	
	total_aligned_depth = sum(chrom_depth_c.values())

	class openLocusClass():
		def __init__(self, name):

			self.name = name
			self.coords = coord_d[name]
			self.contig = self.coords.split(":")[0]
			self.start = int(self.coords.split(":")[1].split("-")[0])
			self.stop  = int(self.coords.split(":")[1].split("-")[1])

			self.seq_c    = Counter()
			self.strand_c = Counter()
			self.size_c   = sizeClass()

		def check_passed(self, pos):
			if pos > self.stop + 1:
				return True
			else:
				return False

		def add(self, read):
			strand = "-" if read.is_reverse else "+"
			self.seq_c[read.get_forward_sequence()] += 1
			self.strand_c[strand] += 1
			self.size_c.update([read.infer_read_length()])

		def write_to_file(self):
			locus_tuple = (self.name, self.contig, self.start, self.stop)
			results_line, gff_line = assessClass().format(locus_tuple, self.seq_c, self.strand_c, self.size_c, total_aligned_depth, 0, project_name)
					
			with open(analysis_file, 'a') as outf:
				print("\t".join(map(str,results_line)), file=outf)

	locus_d = {}
	max_info_length = 15


	for read_i, read in enumerate(bamf.fetch(until_eof=True)):



		if read.is_mapped:
			l_locus = chrom_d[read.reference_name][read.reference_start-1]
			
			strand = "-" if read.is_reverse else "+"

			if l_locus is None:
				locus_name = 'unannotated'

			elif chrom_d[read.reference_name][read.reference_end-1] is None:
				locus_name = 'unannotated'

			else:
				locus_name = l_locus

				if params['reanalyze']:
					if locus_name not in locus_d:
						locus_d[locus_name] = openLocusClass(locus_name)

					locus_d[locus_name].add(read)

		else:
			strand = "*"
			locus_name  = 'unaligned'






		deep_c[(locus_name, read.get_tag("RG"), read.infer_read_length(), strand)] += 1

		if read_i % 1000000 == 0:

			if read.reference_name is not None:
				info = read.reference_name
			else:
				info = locus_name

			if len(info) > max_info_length:
				max_info_length = len(info)

			whitespace = " " * (max_info_length - len(info) + 4)

			print(f"  {info}{whitespace}{read_i:,}  ", end='\r')

			for key in list(locus_d.keys()):

				if locus_d[key].check_passed(read.reference_start):
					locus_d[key].write_to_file()

				del locus_d[key]


	for key in list(locus_d.keys()):

		if locus_d[key].check_passed(read.reference_start):
			locus_d[key].write_to_file()

		del locus_d[key]






	missed_rg     = set(rgs)
	missed_length = set(range(15,31))
	missed_strand = set(["-","+", "*"])
	missed_loci   = set(['unannotated', 'unaligned'])




	for name, locus in loci + [('unannotated', "*"), ("unaligned", "*")]:

		count_line = [name, locus]

		c = Counter()
		with open(deep_counts_file, 'a') as deepf:


			for rg in rgs:
				try:
					cond = rev_conditions[rg]
				except KeyError:
					cond = 'None'

				for length in range(15,31):
					for strand in {"+", "-", "*"}:
						count = deep_c[(name, rg, length, strand)]
						c[rg] += 1

						if count > 0 or include_zeros:
							if count == 0 and strand == "*":
								continue

							print(name, cond, rg, length, strand, count, sep='\t', file=deepf)

							try:
								missed_strand.remove(strand)
							except KeyError:
								pass

							try:
								missed_rg.remove(rg)
							except KeyError:
								pass

							try:
								missed_length.remove(length)
							except KeyError:
								pass



		with open(counts_file, 'a') as countf:
			count_line += [c[rg] for rg in rgs]
			print("\t".join(map(str, count_line)), file=countf)



	with open(deep_counts_file, 'a') as outf:
		for rg in missed_rg:
			cond = rev_conditions[rg]
			print(name, cond, rg, '15', '-', '0', sep='\t', file=outf)
		for strand in missed_strand:
			cond = rev_conditions[rg]
			print(name, cond, rg, '15', strand, '0', sep='\t', file=outf)
		for length in missed_length:
			print(name, cond, rg, length, '-', '0', sep='\t', file=outf)





	counts_file.rename(Path(output_directory, 'counts', f'{name_str}counts.txt'))
	deep_counts_file.rename(Path(output_directory, 'counts', f'{name_str}deepcounts.txt'))
	if params['reanalyze']:
		analysis_file.rename(Path(output_directory, 'counts', f'{name_str}loci.txt')) 


	rmtree(temp_dir)

	print()




