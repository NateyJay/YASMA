
import click
from click_option_group import optgroup

from .generics import *
from .cli import cli

import random
import time




@cli.command(group='Processing', help_priority=3)



@optgroup.group('\n  Basic options',
				help='')


@optgroup.option('-tl', "--trimmed_libraries", 
	required=False, 
	type=click.UNPROCESSED, callback=validate_glob_path,
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

@optgroup.option('-c', '--cores',
	default=4,
	help='Number of cores to use for alignment with bowtie.')


@optgroup.option('--compression',
	default='bam',
	type=click.Choice(['cram', 'bam']),
	help="Compression algorithm used for resulting alignment. Cram is more space efficient, but Bam is more robust/portable.")


@optgroup.option('--override', is_flag=True, default=False, help='Overrides config file changes without prompting.')




def nativealign(**params):
	'''Native aligner'''


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

	locality = 50
	half_locality = round(locality/2)

	max_multi  = 50
	max_random = 3



	align_folder = Path(output_directory, 'nativealign')
	align_folder.mkdir(parents=True, exist_ok=True)

	unsorted_bam = Path(align_folder, "alignment.unsorted.bam")
	sorted_bam = Path(align_folder, "alignment.bam")

	project_table = Path(align_folder, "project_stats.txt")
	library_table = Path(align_folder, "library_stats.txt")

	log_file = Path(output_directory,"nativealign/log.txt")
	sys.stdout = Logger(log_file)



	bowtie_build_index = genome_file.with_suffix(".1.ebwt")

	if not bowtie_build_index.is_file():
		call = ['bowtie-build', genome_file, genome_file.with_suffix('')]

		p = Popen(call, encoding=ENCODING, stdout=PIPE, stderr=PIPE)
		p.wait()


	def get_rg(lib):
		extensions = "".join(lib.suffixes)
		return str(lib.name.removesuffix(extensions))


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
			lib_sizes.append(round(int(line.strip())/4))
			p.wait()

		return(lib_sizes)

	print()
	print("Getting library sizes...")
	lib_sizes = get_lib_sizes()
	total_reads = sum(lib_sizes)
	print("  library sizes:", lib_sizes)
	print(f"  total_reads: {total_reads:,}")
	print()

	def get_unique_weighting():
		## unique alignment

		unique_d = dict()
		with open(genome_file.with_suffix(genome_file.suffix + ".fai"), 'r') as f:
			for line in f:

				ref, length, _, _, _ = line.strip().split()

				unique_d[ref] = [0] * int(length)


		for lib in trimmed_libraries:
			print(" ", str(lib))
			# print()


			call = ['bowtie', '-q', '-v', '1', '-p', str(cores), '-S', '-m', '1', '--best', '--strata', '-x', str(genome_file.with_suffix('')), str(lib)]
			# print(" ".join(call))

			p = Popen(call, encoding=ENCODING, stdout=PIPE, stderr=PIPE)

			for r in range(len(unique_d.keys())+2):
				p.stdout.readline()



			for line in p.stdout:

				line = line.strip().split('\t')
				# print(line)
				# input()

				flag = line[1]
				if flag != '4':

					ref  = line[2]
					loc  = int(line[3])
					size = int(line[5][:-1])

					for r in range(loc-half_locality, loc+half_locality+size):
						try:
							unique_d[ref][r] += 1
						except IndexError:
							pass


		return(unique_d)


	print("Getting unique neighborhood weighting... (ShortStack-U)")
	unique_d = get_unique_weighting()



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


	def print_progress(read_i, map_c, rg, done_rgs, terminal_only=False):
		read_p = round(read_i / total_reads * 100,1)

		counts = [map_c[c] for c in ['U','P','R','Q','H','N']]
		percs  = [round(c/read_i*100,1) for c in counts]


		to_print = '  libraries:\n'
		for lib in trimmed_libraries:
			if get_rg(lib) in done_rgs:
				done = 'x'
			elif rg == get_rg(lib):
				done = '~'
			else:
				done = ' '
			to_print += f"    [{done}] {get_rg(lib)}     \n"

		to_print += f'''  
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
			sys.stdout.overwrite_lines(text=to_print)


		sys.stdout.write(to_print + '\r', terminal_only = terminal_only)
		sys.stdout.flush()


	def do_weighted_alignment():
		read_i = 0
		threshold_i = 0
		lib_c = Counter()
		map_c = Counter()
		done_rgs = set()

		for lib in trimmed_libraries:
			rg = get_rg(lib)
			# print(rg)


			call = ['bowtie', '-q', '-v', '1', '-p', str(cores), '-S', '-a', '-m', str(max_multi), '--best', '--strata', '-x', str(genome_file.with_suffix('')), str(lib)]

			# print()
			# print(" ".join(call))
			# print()

			p = Popen(call, encoding=ENCODING, stdout=PIPE, stderr=PIPE)


			for r in range(len(unique_d.keys())+3):
				p.stdout.readline()


			for line in p.stdout:
				read_i += 1

				a = pysam.AlignedSegment()
				a= a.fromstring(line.strip(), bamfile.header)

				## some useful pysam properties
				# a.flag
				# a.reference_name
				# a.reference_length
				# a.get_tag("XM")
				# a.query_length
				# a.is_mapped
				# a.query_alignment_start
				# a.query_alignment_end




				if a.flag == 4:
					## non-mappers and excluded

					if a.get_tag("XM") == 0:
						a.set_tag("XY","N","Z")
					else:
						a.set_tag("XY","H","Z")

					a.set_tag("XZ",0.0,'f')

				else:
					alignment_count = a.get_tag("XM")-1


					if alignment_count == 1:
						a.set_tag("XY","U","Z")
						a.set_tag("XZ",1.0,'f')

					else:

						weight  = max([unique_d[a.reference_name][a.query_alignment_start], unique_d[a.reference_name][a.query_alignment_end]])

						weights = [weight]

						# alns.append((ref, loc))

						for r in range(alignment_count):
							if r == 0:
								alns = [a]
								weights = []

							else:
								line = p.stdout.readline().strip()
								a = a.fromstring(line, bamfile.header)
								try:
									unique_d[a.reference_name]
								except KeyError:
									print(line)
									print(a)
									print(a.reference_name)
									sys.exit("WEIRD ERROR CAPTURE - report to nate please!!!")
								alns.append(a)


							weight = max([unique_d[a.reference_name][a.query_alignment_start], unique_d[a.reference_name][a.query_alignment_end]])
							weights.append(weight)
							# alns.append((new_ref, new_loc))


						same_weight = len(set(weights)) == 1

						if alignment_count >= max_random and same_weight:
							a.set_tag("XY","Q","Z")
							a.set_tag("XZ",round(1/len(weights),3), 'f')


						else:

							if sum(weights) == 0:
								weights = [1] * len(weights)

							choice = random.choices(range(alignment_count), weights, k=1)[0]
							a = alns[choice]

							if same_weight:
								a.set_tag("XY","R","Z")
								XZ = round(1/sum(weights),3)

							else:
								a.set_tag("XY","P","Z")
								a.set_tag("XZ",round(weights[choice]/sum(weights),3), 'f')


							# print(choice)
							# print(weights)

				map_c[a.get_tag("XY")] += 1
				lib_c[(rg, a.get_tag("XY"))] += 1

				a.set_tag("RG", rg, "Z")
				bamfile.write(a)

				if read_i >= threshold_i:
					threshold_i += 100000
					if threshold_i > total_reads:
						threshold_i = total_reads

					print_progress(read_i, map_c, rg, done_rgs, terminal_only=True)

			done_rgs.add(rg)

		bamfile.close()

		print_progress(read_i, map_c, None, done_rgs)

		# print()
		# pprint(map_c)
		# pprint(lib_c)

		return(map_c, lib_c)

	print()
	print("Making alignment with multimapper placement...")
	map_c, lib_c = do_weighted_alignment()


	with open(project_table, 'w') as outf:
		print("project\tumap\tmmap_wg\tmmap_nw\txmap_nw\txmap_ma\txmap_nv", file=outf)

		to_print = [ic.inputs['project_name']]
		to_print += [map_c[i] for i in ['U','P','R','Q','H','N']]

		print("\t".join(map(str, to_print)), file=outf)


	with open(library_table, 'w') as outf:
		print("project\tlibrary\tumap\tmmap_wg\tmmap_nw\txmap_nw\txmap_ma\txmap_nv", file=outf)

		for lib in trimmed_libraries:
			rg = lib.stem
			to_print = [ic.inputs['project_name'], rg]
			to_print += [lib_c[(rg, i)] for i in ['U','P','R','Q','H','N']]

			print("\t".join(map(str, to_print)), file=outf)




	def print_elapsed(start):

		end = time.time()
		s_elapsed = round(end - start) % 60
		m_elapsed = round((end-start)/ 60) % 60
		h_elapsed = round((end-start)/ 60 / 60)
		print(f"time elapsed: {h_elapsed:02}:{m_elapsed:02}:{s_elapsed:02}")

	print_elapsed(start)


	start = time.time()
	print()
	print("Sorting...")
	pysam.sort("-o", str(sorted_bam), str(unsorted_bam))

	print_elapsed(start)


	os.remove(unsorted_bam)


	print()
	print("Writing table of abundance by library + reference + strand + length...")
	make_depth_file(sorted_bam, verbose=False)


	ic.inputs['alignment_file'] = sorted_bam.absolute()
	ic.write()
	print()
			# input()







