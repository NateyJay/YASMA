
import click
from click_option_group import optgroup

from .generics import *
from .cli import cli

from itertools import chain
from random import choices
import time





@cli.command(group='Processing', help_priority=3)



@optgroup.group('\n  Basic options',
				help='')


@optgroup.option('-tl', "--trimmed_libraries", 
	required=False, 
	type=click.UNPROCESSED, callback=validate_library_paths,
	multiple=True,
	help='Path to trimmed libraries. Accepts wildcards (*).')


@optgroup.option("-g", "--genome_file", 
	# default=f"Annotation_{round(time())}", 
	required=False,
	# type=click.Path(exists=True),
	type=click.UNPROCESSED, callback=validate_path,
	help='Genome or assembly which was used for the original alignment.')


@optgroup.option("-o", "--output_directory", 
	# default=f"Annotation_{round(time())}", 
	required=True,
	type=click.Path(),
	help="Directory name for annotation output.")


@optgroup.group('\n  Bowtie options',
				help='')

@optgroup.option('--cores',
	default=4,
	help='Number of cores to use for alignment with bowtie.')


@optgroup.option('--max_multi',
	default=50,
	help='The maximum number of possible mapping sites for a valid read.')


@optgroup.option('--max_random',
	default=3,
	help='Reads with no weighting will be unmapped if they exceed this number.')


@optgroup.option('--unique_locality',
	default=50,
	help='Window size in nucleotides for unique weighting.')

@optgroup.option('--compression',
	default='bam',
	type=click.Choice(['cram', 'bam']),
	help="Compression algorithm used for resulting alignment. Cram is more space efficient, but Bam is more robust/portable.")


@optgroup.option('--override', is_flag=True, default=False, help='Overrides config file changes without prompting.')






def align(**params):
	'''Aligner based on shortstack3'''


	rc = requirementClass()
	# rc.add_samtools()
	rc.add_bowtie()
	# rc.add_shortstack()
	# rc.add_rnafold()
	rc.check()

	ic = inputClass(params)
	ic.check(['trimmed_libraries', 'genome_file'])


	output_directory        = str(ic.output_directory)
	trimmed_libraries       = ic.inputs['trimmed_libraries']
	genome_file             = ic.inputs['genome_file']

	cores                   = params['cores']
	compression             = params['compression']
	max_multi               = params['max_multi']
	max_random              = params['max_random']
	locality                = params['unique_locality']

	# pprint(params)

	
	half_locality = round(locality/2)


	for tool, _, version in rc.reqs:
		if tool == 'bowtie':
			bowtie_version = float(version.replace(".",""))/100




	align_folder = Path(output_directory, 'align')
	align_folder.mkdir(parents=True, exist_ok=True)

	unsorted_bam = Path(align_folder, "alignment.unsorted.bam")
	sorted_bam = Path(align_folder, "alignment.bam")

	project_table = Path(align_folder, "project_stats.txt")
	library_table = Path(align_folder, "library_stats.txt")

	log_file = Path(output_directory,"align/log.txt")
	sys.stdout = Logger(log_file)



	bowtie_build_index = genome_file.with_suffix(".1.ebwt")

	if not bowtie_build_index.is_file():
		call = ['bowtie-build', genome_file, genome_file.with_suffix('')]
		print(f"bowtie index file not found '{bowtie_build_index}'")
		print(f"building de novo...")

		print(" ".join(map(str, call)))

		p = Popen(call, encoding=ENCODING, stdout=PIPE, stderr=PIPE)
		p.wait()


	def get_rg(lib):
		lib = Path(lib)
		extensions = "".join(lib.suffixes)

		while lib.suffix in {'.gz', '.zip', '.t', '.fastq', '.fq', '.fasta', '.fq', '.fna'}:
			lib = lib.with_suffix("")

		return str(lib.name)


	start = time.time()


	def get_lib_sizes():
		lib_sizes = []

		for lib in trimmed_libraries:
			if lib.suffix == ".gz":
				call = ['gzip', '-cd', lib]
				p0 = Popen(call, encoding=ENCODING, stdout=PIPE, stderr=PIPE)

				call = ['wc', '-l']
				p = Popen(call, encoding=ENCODING, stdin=p0.stdout, stdout=PIPE, stderr=PIPE)

			else:
				call = ['wc', '-l', lib] 
				p = Popen(call, encoding=ENCODING, stdout=PIPE, stderr=PIPE)

			line = p.stdout.readline()

			if ".fq" in lib.suffixes or ".fastq" in lib.suffixes:
				lib_sizes.append(round(int(line.strip().split()[0])/4))
			elif ".fa" in lib.suffixes or ".fasta" in lib.suffixes:
				lib_sizes.append(round(int(line.strip().split()[0])/2))

			p.wait()

		return(lib_sizes)

	print()
	print("Getting library sizes...", flush=True)
	lib_sizes = get_lib_sizes()
	total_reads = sum(lib_sizes)
	print("  library sizes:", lib_sizes, flush=True)
	print(f"  total_reads: {total_reads:,}", flush=True)
	print()




	def make_bam_header():
		header = dict()
		header['HD'] = {'VN': '1.0', 'SO':'coordinate'}

		header['SQ'] = []
		with open(genome_file.with_suffix(genome_file.suffix + ".fai"), 'r') as f:
			for line in f:

				ref, length, _, _, _ = line.strip().split()

				header['SQ'].append({'LN': int(length), 'SN': ref})


		header['RG'] = []
		for lib in trimmed_libraries:

			rg = get_rg(lib)
			header['RG'].append({'ID' : rg})

		return header

	header = make_bam_header()



	bamfile = pysam.AlignmentFile(unsorted_bam, "wb", header=header)





	def bowtie_generator(lib, mmap):

		bowtie_call = ['bowtie']

		if ".fa" in lib.suffixes or ".fasta" in lib.suffixes:
			bowtie_call.append('-f')
			suff = '.fa'

		elif ".fq" in lib.suffixes or ".fastq" in lib.suffixes:
			bowtie_call.append("-q")
			suff = '.fq'

		else:
			sys.exit(f'unknown library suffixes: {lib.suffixes}. Are you sure this is a library?')

		stem = lib.stem.rstrip(''.join(lib.suffixes))
		max_file = Path(align_folder, get_rg(lib) + '.max' + suff)


		if mmap == 'unique':
			bowtie_call += ['-v', '1', '-p', str(cores), '-S', '-m', '1', '--best', '--strata', '--max', str(max_file)]

		elif mmap == 'multi':
			lib = max_file
			bowtie_call += ['-v', '1', '-p', str(cores), '-S', '-a', '--best', '--strata']


		if bowtie_version >= 1.3:
			bowtie_call.append("-x")

		bowtie_call.append(str(genome_file.with_suffix('')))


		if ".gz" in lib.suffixes:

			call = ['gzip', '-cd', str(lib)]
			gzip = Popen(call, stdout=PIPE, encoding=ENCODING)

			bowtie_call.append("-")
			p = Popen(bowtie_call, encoding=ENCODING, stdout=PIPE, stderr=PIPE, stdin=gzip.stdout)

			# print(" ".join(call), "|", " ".join(bowtie_call))

		else:
			bowtie_call.append(str(lib))
			p = Popen(bowtie_call, encoding=ENCODING, stdout=PIPE, stderr=PIPE)

			# print(" ".join(bowtie_call))



		while True:
			line = p.stdout.readline().strip()

			if line.startswith("@"):
				continue

			if line == '':
				break


			a = pysam.AlignedSegment()
			try:
				a = a.fromstring(line, bamfile.header)
			except ValueError as err:
				print(err)
				print(f'call: {call}')
				print(f'line: {line}')
				raise 

			yield a

		p.wait()

		if mmap == 'multi':
			max_file.unlink()


	
	unique_d = dict()
	with open(genome_file.with_suffix(genome_file.suffix + ".fai"), 'r') as f:
		for line in f:

			ref, length, _, _, _ = line.strip().split()

			unique_d[ref] = [0] * int(length)





	def print_progress(read_i, map_c, current, done, status_message, terminal_only=False):
		read_p = round(read_i / total_reads * 100,1)

		counts = [map_c[c] for c in ['U','P','R','Q','H','N']]
		percs  = [round(c/total_reads*100,1) for c in counts]


		# to_print = f'\n  command:\n    {" ".join(call)}\n\n  libraries:\n'

		to_print = f"\n#### performing alignment ####\n\n"
		to_print += f'  libraries:\n'
		to_print += f'     u   m  \n'
		for lib in trimmed_libraries:
			rg = get_rg(lib)


			glyphs = []
			for step in ['unique','mmap']:

				if (rg, step) in done:
					glyphs.append("x")
				elif (rg, step) == current:
					glyphs.append("~")
				else:
					glyphs.append(" ")


			to_print += f"    [{glyphs[0]}] [{glyphs[1]}] {rg}     \n"

		to_print += f'''  
  status: {status_message}                   

  current read:\t{read_i} ({read_p}%)       
                                   maptag\tmapcat\t perc\treads
  (unique mappers) ............... XY:Z:U\tumap\t {percs[0]}%\t{counts[0]:,}         
  (mmap, weighted) ............... XY:Z:P\tmmap_wg\t {percs[1]}%\t{counts[1]:,}        
  (mmap, placed w/o weighting) ... XY:Z:R\tmmap_nw\t {percs[2]}%\t{counts[2]:,}            
  (nonmap, above rand_max) ....... XY:Z:Q\txmap_nw\t {percs[3]}%\t{counts[3]:,}          
  (nonmap, above max alignments) . XY:Z:H\txmap_ma\t {percs[4]}%\t{counts[4]:,}          
  (nonmap, no valid alignments ... XY:Z:N\txmap_nv\t {percs[5]}%\t{counts[5]:,}         
					'''

		if read_i > 1:
			sys.stdout.overwrite_lines(text=to_print.rstrip())


		sys.stdout.write(to_print + '\r' , terminal_only = terminal_only)
		sys.stdout.flush()





	map_c = Counter()
	lib_c = Counter()


	read_i = 0
	threshold_i = 0
	done_rgs = set()
	status_message = 'not defined'

	for lib in trimmed_libraries:
		rg = get_rg(lib)
		# print(rg, 'uniq + unal')

		gen = bowtie_generator(lib, mmap='unique')

		for a in gen:
			read_i += 1

			if a.is_mapped:

				left  = a.query_alignment_start - half_locality
				right = a.query_alignment_end   + half_locality

				for r in range(left, right):
					try:
						unique_d[a.reference_name][r] += 1
					except IndexError:
						pass

				a.set_tag("XY","U","Z")
				a.set_tag("XZ",1.0,'f')

			else:
				a.set_tag("XY","N","Z")
				a.set_tag("XZ",1.0,'f')


			map_c[a.get_tag("XY")] += 1
			lib_c[(rg, a.get_tag("XY"))] += 1

			a.set_tag("RG", rg, "Z")
			bamfile.write(a)

			if read_i >= threshold_i:
				threshold_i += 100000
				if threshold_i > total_reads:
					threshold_i = total_reads

				print_progress(read_i, map_c, (rg, 'unique'), done_rgs, status_message=f'{rg} unique mapping', terminal_only=True)


		done_rgs.add((rg, 'unique'))

	# pprint(map_c)

	for lib in trimmed_libraries:

		rg = get_rg(lib)
		# print(rg, 'mmap')

		gen = bowtie_generator(lib, mmap='multi')

		for a in gen:

			read_i += 1

			qname = a.query_name

			alignment_count = a.get_tag("XM")-1 ## Bowtie reports XM as +1 over the number of reported alignments

			if alignment_count > max_multi:

				a.set_tag("XY","H","Z")
				a.set_tag("XZ",0.0,'f')

				a.flag = 4
				a.reference_name = '*'
				a.reference_start = -1
				a.is_mapped = False

				## Clearing the other alignments for this read
				for r in range(alignment_count-1):
					next(gen)
		
			else:

				weight  = max([unique_d[a.reference_name][a.query_alignment_start], unique_d[a.reference_name][a.query_alignment_end]])
				weights = [weight]

				alns    = [a]

				# print("", read_count, a.query_name, f"{a.reference_name}:{a.query_alignment_start}", weight, sep='\t')
				


				for r in range(alignment_count-1):

					a = next(gen)

					# read_count += 1
					weight  = max([unique_d[a.reference_name][a.query_alignment_start], unique_d[a.reference_name][a.query_alignment_end]])
					weights.append(weight)

					# print("  ", read_count, a.query_name, f"{a.reference_name}:{a.query_alignment_start}", weight, sep='\t')
					alns.append(a)	


					if a.query_name != qname:
						# break
						print("qname mismatch!")
						print(f"expected: {qname}")
						print(f"found:    {a.query_name}")
						print(r+2, "<- alignment number")
						print(alignment_count, "<- total expected alignments")
						for a in alns:
							print("  ", a)
						print("WEIRD ERROR 1 - please report to nate!")
						sys.exit()



				same_weight = len(set(weights)) == 1

				if alignment_count > max_random and same_weight:

					a = alns[0]
					
					a.set_tag("XY","Q","Z")
					a.set_tag("XZ",round(1/len(weights),3), 'f')
					a.flag = 4
					a.reference_name = '*'
					a.reference_start = -1
					a.is_mapped = False

				else:

					if sum(weights) == 0:
						weights = [1] * len(weights)


					choice = choices(range(alignment_count), weights, k=1)[0]
					a = alns[choice]

					if same_weight:
						a.set_tag("XY","R","Z")
						XZ = round(1/sum(weights),3)

					else:
						a.set_tag("XY","P","Z")
						a.set_tag("XZ",round(weights[choice]/sum(weights),3), 'f')


			map_c[a.get_tag("XY")] += 1
			lib_c[(rg, a.get_tag("XY"))] += 1

			a.set_tag("RG", rg, "Z")
			bamfile.write(a)

			if read_i >= threshold_i:
				threshold_i += 100000
				if threshold_i > total_reads:
					threshold_i = total_reads

				print_progress(read_i, map_c, (rg, 'mmap'), done_rgs, status_message=f'{rg} mmap mapping', terminal_only=True)


		done_rgs.add((rg, 'mmap'))

	print_progress(read_i, map_c, (rg, 'mmap'), done_rgs, status_message='done', terminal_only=True)

	bamfile.close()





	with open(project_table, 'w') as outf:
		print("project\tumap\tmmap_wg\tmmap_nw\txmap_nw\txmap_ma\txmap_nv", file=outf)

		to_print = [ic.inputs['project_name']]
		to_print += [map_c[i] for i in ['U','P','R','Q','H','N']]

		print("\t".join(map(str, to_print)), file=outf)


	with open(library_table, 'w') as outf:
		print("project\tlibrary\tumap\tmmap_wg\tmmap_nw\txmap_nw\txmap_ma\txmap_nv", file=outf)

		for lib in trimmed_libraries:
			rg = get_rg(lib)
			to_print = [ic.inputs['project_name'], rg]
			to_print += [lib_c[(rg, i)] for i in ['U','P','R','Q','H','N']]

			print("\t".join(map(str, to_print)), file=outf)




	def print_elapsed(start):

		end = time.time()
		s_elapsed = round(end - start) % 60
		m_elapsed = round((end-start)/ 60) % 60
		h_elapsed = round((end-start)/ 60 / 60)
		print(f"time elapsed: {h_elapsed:02}:{m_elapsed:02}:{s_elapsed:02}")

	print()
	print_elapsed(start)


	start = time.time()
	print()
	print("Sorting...", flush=True)
	pysam.sort("-o", str(sorted_bam), str(unsorted_bam))

	print_elapsed(start)


	unsorted_bam.unlink()


	print()
	print("Writing table of abundance by library + reference + strand + length...", flush=True)
	make_depth_file(sorted_bam, verbose=False)


	ic.inputs['alignment_file'] = sorted_bam.absolute()
	ic.write()
	print()
	print("alignment complete!")


			# input()







