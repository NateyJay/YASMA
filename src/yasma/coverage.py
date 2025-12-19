
from .generics import *
from .track_generics import *




@cli.command(group='Calculation', help_priority=4)


@optgroup.group('\n  Basic options',
				help='')

@optgroup.option("-o", "--output_directory", 
	# default=f"Annotation_{round(time())}", 
	required=False,
	type=click.UNPROCESSED, callback=validate_outdir,
	help="Directory name for annotation output. Defaults to the current directory, with this directory name as the project name.")


@optgroup.option("-a", "--alignment_file", 
	required=False, 
	type=click.UNPROCESSED, callback=validate_path,
	help='Alignment file input (bam or cram).')


@optgroup.group('\n  Run options',
				help='')

@optgroup.option("-p", "--peaks", 
	required=False,
	type=str,
	multiple=True,
	default=["20-25"],
	help='Entry of size limits for a peak. Encoded as two integers separated by a dash (`-`), for example: 21-22 is an appropriate entry for plant miRNAS. Also accepts single-size peaks without dash. Multiple peaks may be identified, separated with spaces or by calling the option again. Default: 20-25')

# @optgroup.option("--summarize", 
# 	default=False,
# 	is_flag=True,
# 	help="Produces a summary of all sizes and annotations")

@optgroup.option("--force", 
	default=False,
	is_flag=True,
	help="Forces remaking bigwigs even if all expected are found")

def coverage(**params):
	"""Produces bigwig coverage files"""

	rc = requirementClass()
	# rc.add_samtools()
	rc.check()

	ic = inputClass(params)
	ic.check(['alignment_file'])

	output_directory        = ic.output_directory
	alignment_file          = ic.inputs["alignment_file"]
	project_name            = ic.inputs['project_name']
	# annotation_readgroups   = ic.inputs['annotation_readgroups']
	conditions              = ic.inputs['conditions']

	force = params['force']

	def process_peaks(unprocessed_peaks):
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

	peaks = process_peaks(list(params['peaks']))


	chromosomes, libraries = get_chromosomes(alignment_file)

	chrom_depth_c = get_global_depth(alignment_file, aggregate_by=['rg'])

	library_to_condition = {}
	for condition, srrs in conditions.items():
		for srr in srrs:
			library_to_condition[srr] = condition

	libraries = list(library_to_condition.keys())



	rpm_d = {}
	for rg, depth in chrom_depth_c.items():
		rpm_d[rg] = 1 / depth * 1000000


	# keys = list(chrom_depth_c.keys())
	# for key in keys:
	# 	if key[0] in libraries:
	# 		chrom_depth_c[key[1]] += chrom_depth_c[key]

	# 	del chrom_depth_c[key]



	# aligned_depth = sum(chrom_depth_c.values())


	cov_dir = Path(output_directory, 'coverage')
	cov_dir.mkdir(parents=True, exist_ok=True)


	def check_done(force):
		for size in sizes:
			for strand in strands:
				for condition in conditions.keys():
					file = Path(cov_dir,f'{condition}_{size}{strand}.bw')
					if not file.is_file():
						return()

		print("All expected output files are already found.")
		if force:
			print("  force=True -> running anyways")
		else:
			print("  stopping (override with --force)")


	bw_d = {}

	sizes = peaks + ["non"]
	strands = ['+', "-"]

	for size in sizes:
		for strand in strands:
			for condition in conditions.keys():
				bw_d[(condition, size, strand)] = bigwigClass(Path(cov_dir, f'{condition}_{size}{strand}.temp.bw'), total_reads=None, chromosomes=chromosomes, strand= strand)





	bamf = pysam.AlignmentFile(alignment_file)
	
	for chrom_count, chrom_and_length in enumerate(chromosomes):



		chrom, chrom_length = chrom_and_length
		print(f"{chrom_count+1} / {len(chromosomes)}")
		print(f"chrom: {chrom} -> {chrom_length:,} bp")



		for key in bw_d.keys():
			bw_d[key].reset(chrom_length)



		for i,read in enumerate(bamf.fetch(contig=chrom)):

			# perc_out = perc.get_percent(i)
			# if perc_out:
			# 	print(f"   reading position depths ..... {perc_out}%", end='\r', flush=True)
			
			if read.is_unmapped:
				continue

			library = read.get_tag("RG")

			if library not in libraries:
				continue

			if read.is_forward:
				strand = "+"
			else:
				strand = "-"

			length    = read.query_length
			position  = read.reference_start
			condition = library_to_condition[library]
			val       = round(rpm_d[library], 4)


			if length in peaks:
				size = length
			else:
				size = 'non'

			bw_d[(condition, size, strand)].add(position, length, val)



		print()
		for condition, size, strand in bw_d.keys():
			name = f"{condition}_{size}{strand}"
			print(f"[{name}]", end='  ', flush=True)
			bw_d[(condition, size, strand)].rle(chrom)

		print()
		print()


	for key in bw_d.keys():
		bw_d[key].close()
		Path(bw_d[key].file).rename(bw_d[key].file.replace(".temp.bw", '.bw'))

	bamf.close()



