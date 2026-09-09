# Simple quantification of features

from .generics import *

from shutil import rmtree
from collections import deque
import gzip




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

@optgroup.option("-an", "--annotation_files", 
	required=False, 
	type=click.UNPROCESSED, callback=validate_glob_path,
	multiple=True,
	help='Locus annotations to count in gff, gff3, gtf, bed, or tabular format. Tabular requires contig:start-stop and locus_name in the first two columns (tab-delimited, "#" escape char). Defaults to find all gff files associated with annotations.')

# @optgroup.option("-n", "--name", 
# 	required=False,
# 	type=str,
# 	help="Name for resulting counts and analysis. Required if not using the default annotation file.")

# @optgroup.option("-f", "--features", 
# 	required=False, 
# 	multiple=True,
# 	help='A filter to only analyze certain features of an input gff/gtf annotations (field 3). Permits all features by default. Multiple features can be included as space delimited entries. Does not apply for other annotation file types.')

# @optgroup.option("--reanalyze", 
# 	default=False,
# 	is_flag=True,
# 	help="Flag calling for the reanalysis of annotated features. Produces a file similar to the 'loci.txt' output of tradeoff, which includes major sRNA locus dimensions.")

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
	ic.check(['alignment_file'])

	output_directory     = ic.output_directory
	alignment_file       = ic.inputs['alignment_file']
	conditions           = ic.inputs['conditions']
	annotation_files     = ic.inputs['annotation_files']
	project_name         = ic.inputs['project_name']

	include_zeros        = params['include_zeros']
	# name                 = params['name']

	Path(output_directory, "counts").mkdir(parents=True, exist_ok=True)

	chromosomes, libraries = get_chromosomes(alignment_file)

	if not conditions:
		print("Warning: no conditions defined. \nUsing library names as conditions (supply these with -c)")

		conditions = {'default': []}
		for library in libraries:
			conditions['default'].append(library)

	rev_conditions = {}
	for k,vs in conditions.items():
		for v in vs:
			rev_conditions[v] = k


	chrom_depth_c = get_global_depth(alignment_file, aggregate_by=['rg','chrom'])


	aligned_read_count = sum(chrom_depth_c.values())


	try:
		annotation_files = ic.inputs['annotation_files']
	except:
		annotation_files = None

	if not annotation_files:
		annotation_files = [f for f in output_directory.glob("*/loci.gff3")]

		if not annotation_files:
			print("Error: No annotation files found! have you run an annotation? (yasma tradeoff)")
			print("Searches for annotations by default at */loci.gff3")
			print("Provide other annotations with -an [annotation_file_path]")
			sys.exit()

	if not annotation_files:
		sys.exit("Error: Annotation files not found. Are they supplied correctly?")


	print()
	print(f'counting annotations from:')

	pass_count = 0
	for af in annotation_files:
		good_ann = False
		if af.is_file():
			with open(af, 'r') as f:
				for line in f:
					if not line.startswith("#"):
						good_ann = True
						pass_count += 1
						break

		print(f"  {good_ann}\t{af}")

	if pass_count == 0:
		sys.exit("Error: none of the input annotations are valid! Check the files, maybe something went wrong with annotation step.")

	annotation_files = [a for a in annotation_files if a.is_file()]

	if not annotation_files:
		sys.exit("Error: All annotation paths are not findable...")
	print()

	c = Counter()
	features = list()



	with pysam.AlignmentFile(alignment_file, "rb") as bamf:
		header = bamf.header.to_dict()
		contigs = [entry['SN'] for entry in header['SQ']]



	annotation_names = set()
	locus_d = dict()
	coord_d = dict()
	used_names = list()

	for annotation_file in annotation_files:
		annotation_file = Path(annotation_file)
		annotation_name = annotation_file.parts[-2]
		annotation_names.add(annotation_name)


		with open(annotation_file, 'r') as f:
			if annotation_file.suffix == '.txt':
				f.readline()

			for i,line in enumerate(f):

				if len(line) > 1 and not line.startswith("#"):
					line = line.strip().split("\t")
					# print(line)

					if annotation_file.suffix == ".txt":
						coords, name = line[:2]
						coords = coords.replace("..", '-')

					elif annotation_file.suffix in ['.gff', '.gff2', '.gff3']:
						coords  = f"{line[0]}:{line[3]}-{line[4]}"
						name    = line[8].split(";")[0].split("=")[-1].strip('"')
						feature = line[2]

						# if feature_filter and feature in feature_filter:
						# 	continue

					elif annotation_file.suffix == '.gtf':
						coords  = f"{line[0]}:{line[3]}-{line[4]}"
						name    = line[8].split(";")[0].split()[-1].strip('"')
						feature = line[2]

						# if feature_filter and feature in feature_filter:
						# 	continue

					elif annotation_file.suffix == '.bed':
						name   = f"bed_{i}"
						coords = f"{line[0]}:{line[1]}-{line[2]}"

					else:
						print(f'file.suffix "{annotation_file.suffix}" not expected in annotation file...')
						sys.exit()


					contig   = coords.split(':')[0]
					contig_i = contigs.index(contig)
					start    = int(coords.split(':')[1].split("-")[0])
					end      = int(coords.split(':')[1].split("-")[1])


					if name in used_names:
						used_names.append(name)
						name = f"{name}_{used_names.count(name)}"
					else:
						used_names.append(name)

					features.append((contig, contig_i, start, end, name, annotation_name))

					try:
						locus_d[annotation_name].append(name)
					except KeyError:
						locus_d[annotation_name] = [name]

					coord_d[(annotation_name, name)] = coords

	# print(annotation_files)
	# print(annotation_names)
	# sys.exit()

	c      = Counter()	
	deep_c = Counter()
	# five_c = Counter()

	read_i = 0
	with pysam.AlignmentFile(alignment_file, "rb") as bamf:
		header = bamf.header.to_dict()

		features.sort(key=lambda x: x[2]) # sorting by left-most position
		features.sort(key=lambda x: x[1]) # sorting by contig order

		features = deque(features)


		last_contig = ''

		for read in bamf:

			if read_i % 100000 == 0:
				p = round(read_i / aligned_read_count * 100, 1)
				print(f"  read_i: {read_i:,}   ({p}%)   ", end='\r')

				# if read_i > 10000000:
				# 	break


			read_i += 1

			size    = read.query_length
			strand  = "-" if read.is_reverse else "+"
			library = read.get_tag("RG")
			seq     = read.get_forward_sequence()
			fivep   = seq[0]
			if seq[0] == 'T':
				fivep = 'U'


			if read.is_unmapped:
				# c[(annotation_name,'unmapped', library)] += 1
				# deep_c[(annotation_name, 'unmapped', library, size, "*")] += 1
				# continue
				break

			# print()
			# print(f"{read.qname}  {read.reference_name}  {read.reference_start}")

			if read.reference_name != last_contig:
				reference_contig_i = contigs.index(read.reference_name)
			last_contig = read.reference_name

			i = 0
			unannotated = annotation_names.copy()

			while True:
				try:
					contig, contig_i, start, end, name, annotation_name = features[i]
				except IndexError:
					break

				# print(f"feature_coords: {annotation_name} {name} {contig}:{start}-{end}")

				if contig_i > reference_contig_i:
					# print('features are ahead of reads (contig)')
					break

				if contig_i < reference_contig_i:
					# print('eliminating passed feature (contig)')
					features.popleft()
					continue

				if start > read.reference_end:
					# print('features are ahead of reads (coord)')
					break

				if end < read.reference_start:
					# print('eliminating passed feature (coord)')
					features.popleft()
					continue

				unannotated.discard(annotation_name)

				c[(annotation_name, name, library)] += 1
				deep_c[(annotation_name, name, library, size, strand, fivep)] += 1
				# five_c[(annotation_name, name, library, fivep)] += 1

				i += 1

			for annotation_name in unannotated:
				name = 'unannotated'
				c[(annotation_name, name, library)] += 1
				deep_c[(annotation_name, name, library, size, strand, fivep)] += 1
				# five_c[(annotation_name, name, library, fivep)] += 1

	outputs      = dict()
	output_files = dict()

	for annotation_name in annotation_names:

		key = (annotation_name, 'counts')
		output_files[key] = Path(output_directory, 'counts', f'{annotation_name}.counts.temp')
		outputs[key]      = open(output_files[key], 'w')
		print('name', 'locus', "\t".join(libraries), sep='\t', file=outputs[key])

		key = (annotation_name, 'deepcounts')
		output_files[key] = Path(output_directory, 'counts', f'{annotation_name}.deepcounts.temp')
		outputs[key]      = gzip.open(output_files[key], 'wt')
		print('name', 'condition', 'library','length','strand', 'fivep','count', sep='\t', file=outputs[key])

		# key = (annotation_name, 'fivep')
		# output_files[key] = Path(output_directory, 'counts', f'{annotation_name}.fivep.temp')
		# outputs[key]      = open(output_files[key], 'w')
		# print('name', 'condition', 'rg','fivep', 'count', sep='\t', file=outputs[key])




	for annotation_name in annotation_names:
		loci = locus_d[annotation_name]
		# loci += ['unannotated','unmapped']
		loci += ['unannotated']

		for locus in loci:

			counts = []

			if locus == 'unmapped':
				strands = '*'
			else:
				strands = ['+', '-']
			for library in libraries:

				try:
					condition = rev_conditions[library]
				except KeyError:
					continue


				# for five in ['A','U','C','G']:
				# 	count = five_c[(annotation_name, locus, library, five)]
				# 	print(locus, condition, library, five, count, sep='\t', file=outputs[(annotation_name, 'fivep')], flush=True)


				counts.append(c[(annotation_name, locus, library)])


				for size in range(15,31):

					for strand in strands:

						for fivep in ['A','U','G','C']:

							count = deep_c[(annotation_name, locus, library, size, strand, fivep)]

							if count == 0 and not include_zeros:
								continue

							print(locus, condition, library, size, strand, count, fivep, sep='\t', file=outputs[(annotation_name, 'deepcounts')], flush=True)

			if locus != 'unannotated':
				coords = coord_d[(annotation_name, locus)]
			else:
				coords = '*'




			print(locus, coords, "\t".join(map(str, counts)), sep='\t', file=outputs[(annotation_name, 'counts')])


	for f in outputs.values():
		f.close()

	print()
	print()
	print("saving outputs:")
	for file_name in output_files.values():

		if "deepcounts" in file_name.name:
			suff = ".txt.gz"
		else:
			suff = ".txt"
		print(" ", file_name.with_suffix(suff))
		file_name.rename(file_name.with_suffix(suff))

	ic.inputs['annotation_files'] = annotation_files
	ic.write()



