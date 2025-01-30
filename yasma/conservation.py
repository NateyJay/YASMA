from .generics import *

import gzip
from shutil import rmtree
from time import sleep



@cli.command(group="Calculation", help_priority=3)

@optgroup.group('\n  Basic options', help='')


@optgroup.option("-o", "--output_directory", 
	# default=f"Annotation_{round(time())}", 
	required=False,
	type=click.UNPROCESSED, callback=validate_outdir,
	help="Directory name for annotation output. Defaults to the current directory, with this directory name as the project name.")

@optgroup.option("-g", "--genome_file", 
	required=False,
	type=click.UNPROCESSED, callback=validate_path,
	help='Genome or assembly which was used for the original alignment.')

@optgroup.option("-af", "--annotation_file",
	type=click.UNPROCESSED, callback=validate_path,
	help="gff3 annotation for the query project.")


@optgroup.option("-tp", "--target_project",
	type=click.UNPROCESSED, callback=validate_path,
	help="Path to an YASMA project folder that comtains a completed YASMA-tradeoff annotation. If [genome_file] or [annotation_file] are not found in its inputs.json, this give an error.")

@optgroup.option("-ta", "--target_annotation",
	type=click.UNPROCESSED, callback=validate_path,
	help="Path to an sRNA annotation (.gff3) which will be used for the conservation search.")

@optgroup.option("-tg", "--target_genome",
	type=click.UNPROCESSED, callback=validate_path,
	help="Path to the accompanying genome (.fa).")

@optgroup.option("-tn", "--target_name",
	type=str,
	help="A name for the target and comparison")




def conservation(**params):
	'''Checks for conservation between annotations'''

	rc = requirementClass()
	rc.add_hmmer()
	rc.check()


	ic = inputClass(params)
	ic.check(['genome_file'])
	ic.check(['annotation_file'])

	output_directory  = ic.output_directory

	genome_file       = ic.inputs['genome_file']
	annotation_file   = ic.inputs['annotation_file']


	if params['target_project']:
		target_inputs_file = Path(params['target_project'], 'inputs.json')

		if not target_inputs_file.is_file():
			sys.exit("Error: target_project does not have a readable inputs.json")

		with open(target_inputs_file, 'r') as f:
			target_inputs = json.load(f) 

		target_annotation = Path(target_inputs['annotation_file']).is_relative_to(output_directory)
		target_genome     = Path(target_inputs['genome_file']).is_relative_to(output_directory)
		target_name       = target_inputs['project_name']

	elif params['target_annotation'] and params['target_genome'] and params['target_name']:

		target_annotation = Path(params['target_annotation']).is_relative_to(output_directory)
		target_genome     = Path(params['target_genome']).is_relative_to(output_directory)
		target_name       = params['target_name']

	else:
		sys.exit("Error: without a valid target_project folder specified, target_annotation, _genome, and _name are req'd")


	if not target_annotation.is_file():
		sys.exit(f"Error: target_annotation {target_annotation} file not found")

	if not target_genome.is_file():
		sys.exit(f"Error: target_genome {target_genome} file not found")



	compare_folder = Path(output_directory, 'conservation', target_name)
	compare_folder.mkdir(parents=True, exist_ok=True)

	query_seqs = Path(compare_folder, 'query.fa')


	genf = pysam.FastaFile(genome_file)

	with open(query_seqs, 'w') as outf:
		with open(ic.inputs['annotation_file'], 'r') as f:

			for line in f:
				if line.startswith("#"):
					continue

				line = line.strip().split("\t")

				contig  = line[0]
				feature = line[2]


				if feature == 'otherRNA':
					continue

				start   = line[3]
				stop    = line[4]
				strand  = line[6]
				ID      = line[8].split(';')[0].split("=")[1]


				region = f"{contig}:{start}-{stop}"


				out = genf.fetch(region=region)

				out = out.upper()


				if strand == '-':
					out = out[::-1]
					out = complement(out, dna=True)


				print(f">{ID}", file=outf)
				print(out, file=outf)

				# print(out)

	genf.close()




	### conducting hmmer search


	print()
	print(f"searching target_genome: {target_genome}")
	print(f"    with nhmmer")

	call = ['nhmmer', str(query_seqs), str(target_genome)]

	print(" ".join(map(str, call)))

	p = Popen(call, stdout=PIPE, encoding=ENCODING)


	for line in p.stdout:
		print(line.strip(), flush=True)

	p.wait()





























