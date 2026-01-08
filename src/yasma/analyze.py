

from .generics import *
from pprint import pprint


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

# @optgroup.option("-c", "--conditions", 
# 	required=False, 
# 	multiple=True,
# 	type=click.UNPROCESSED, callback=validate_condition,
# 	help='Values denoting condition groups (sets of replicate libraries) for projects with multiple tissues/treatments/genotypes. Can be entered here as space sparated duplexes, with the library base_name and condition groups delimited by a colon. E.g. SRR1111111:WT SRR1111112:WT SRR1111113:mut SRR1111114:mut')

@optgroup.option("-an", "--annotation_file", 
	required=True, 
	type=click.UNPROCESSED, callback=validate_path,
	multiple=False,
	help='A gff3 annotation to assess.')

@optgroup.option("-n", "--name", 
	required=True,
	type=str,
	help="Name for analysis output folder suffix `analysis_[name]`.")


# @optgroup.option("--ignore_unaligned",
# 	is_flag=False,
# 	help="Include to skip counting unaligned reads in deepcounts.txt. These are useful for some analyses, but it can be faster to ignore.")



def analyze(**params):
	"""Analyzes alignment based on a supplied annotation."""

	rc = requirementClass()
	rc.check()

	ic = inputClass(params)
	ic.check(['alignment_file'])

	output_directory     = ic.output_directory
	alignment_file       = ic.inputs['alignment_file']
	# conditions           = ic.inputs['conditions']
	project_name         = ic.inputs['project_name']

	name                 = params['name']
	annotation_file      = params['annotation_file']

	annotation_name = f"analysis_{name}"

	analysis_dir = Path(output_directory, annotation_name)
	analysis_dir.mkdir(parents=True, exist_ok=True)


	### Getting alignment metadata

	chromosomes, libraries = get_chromosomes(alignment_file)


	### getting alignment depths by library and chromosome

	chrom_depth_c = get_global_depth(alignment_file, aggregate_by=['rg','chrom'])
	aligned_read_count = sum(chrom_depth_c.values())
	print()
	print(f"  {aligned_read_count:,} aligned reads")



	## initializing outputs

	def init_gff(file_name):
		with open(file_name, 'w') as outf:
			print("##gff-version 3", file=outf)

			for chrom, chrom_length in chromosomes:
				print(f"##sequence-region   {chrom} 1 {chrom_length}", file=outf)

	gff_file = Path(analysis_dir, 'loci.gff3')
	init_gff(gff_file)

	results_file = Path(analysis_dir, 'loci.txt')
	with open(results_file, 'w') as outf:
		print("\t".join(assessClass().header), file=outf)


	reads_file = Path(analysis_dir, 'reads.txt')
	with open(reads_file, 'w') as outf:
		print("cluster\tseq\trank\tdepth\trpm\tlocus_prop", file=outf)


	def top_reads_save(read_c, file, name):

		cum_count = 0
		top_reads = read_c.most_common(100)
		for rank, read in enumerate(top_reads):

			seq, dep = read
			rpm = round(dep / aligned_read_count * 1000000, 4)

			cum_count += dep

			loc_prop = round(cum_count / sum(read_c.values()), 4)

			with open(file, 'a') as outf:
				print(name, seq, rank, dep, rpm, loc_prop, file=outf, sep='\t')

				if loc_prop >= 0.3:
					break

	## analyzing from gff3

	resf = open(results_file, 'a')
	gfff = open(gff_file, 'a')



	# print(annotation_file)
	print()
	locus_i = 0
	last_contig = ''
	sizecall_summary = Counter()
	bamf = pysam.AlignmentFile(alignment_file,'rb')
	with open(annotation_file, 'r') as f:

		for line in f:
			if line.startswith("#"):
				continue

			line = line.strip().split('\t')

			contig = line[0]
			start  = int(line[3])
			stop   = int(line[4])
			ID     = line[8].split(";")[0].lstrip("ID=")

			seqs    = Counter()
			sizes   = sizeClass()
			strands = Counter()

			if contig != last_contig:
				last_stop = 0

			for read in bamf.fetch(contig=contig, start=start, stop=stop):
				if read.is_unmapped:
					continue

				lib = read.get_tag("RG")
				# if lib not in annotation_libraries:
				# 	continue

				sam_seq = read.get_forward_sequence().replace("T","U")
				seqs[sam_seq] += 1

				strand = "+" if read.is_forward else "-"
				strands[strand] += 1
				sizes.update(read.qlen)

			# print(seqs)
			# print(sizes)
			# print(strands)

			locus_i += 1
			print(ID, end='\r', flush=True)
			locus = [ID, contig, start, stop]

			results_line, gff_line = assessClass().format(locus, seqs, strands, sizes, aligned_read_count, last_stop, project_name, annotation_name)

			# print(results_line)
			# print(gff_line)


			print("\t".join(map(str, results_line)), file=resf)
			print("\t".join(map(str, gff_line)), file=gfff)

			top_reads_save(seqs, reads_file, name)

			if sum(strands.values()) == 0:
				sizecall_summary["None"] += 1
			else:
				sizecall_summary[str(sizes)] += 1

			last_stop   = stop
			last_contig = contig

	bamf.close()
	resf.close()
	gfff.close()

	print()
	print()
	print(f'reanalyzed: {locus_i:,} loci')

	print()
	print("sizecalls found:")
	for k,v in sizecall_summary.most_common():

		print(f"  {k} {(9-len(k)) * ' '} {v:,}")

	# pprint(sizecall_summary.most_common())












