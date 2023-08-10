


import sys
import os

import click
from click_option_group import optgroup

from pathlib import Path
from os.path import isfile, isdir
from collections import Counter#, deque
from pprint import pprint

# from random import sample

# import numpy as np
# from statistics import quantiles
# import math
import shutil
import re

from .generics import *
from .cli import cli

# from statistics import mean, median, StatisticsError

# from time import time



def get_log_details(output_directory):
	bpj_table = f"{output_directory}/align/project_stats.txt"
	bpjf = open(bpj_table, 'w')
	print("project\tXY:Z:N\tXY:Z:M\tXY:Z:O\tXY:Z:U\tXY:Z:R\tXY:Z:P", file=bpjf)


	# Create and open the "01out-srr_alignments.txt" file
	srr_table = f"{output_directory}/align/library_stats.txt"
	srrf = open(srr_table, 'w')
	print("project\tsrr\tunique\tmulti\tnon", file=srrf)



	# Construct the path to the corresponding log.txt
	path_to_log = f"{output_directory}/align/log.txt"
	# print("looking for:", path_to_log)
	# print("  -> ", isfile(path_to_log))
	
	


	# Open the log file
	with open(path_to_log, 'r') as f:
		for line in f:
			# print("->", line.strip())
			if "Completed. Results" in line:

				srr = line.split('/temp/')[-1].split('_readsorted.sam.gz')[0]
				# f.readline()
				# f.readline()
				unique = f.readline().split()[2]
				multi = f.readline().split()[2]
				non = f.readline().split()[2]
				
				

				# Print the values to the "01out-srr_alignments"
				print(output_directory, srr, unique, multi, non, sep='\t', file=srrf)

			elif "XY:Z:" in line:
				vals = []
				vals.append(line.strip().split()[-5])
				vals.append(f.readline().strip().split()[-5])
				vals.append(f.readline().strip().split()[-5])
				vals.append(f.readline().strip().split()[-5])
				vals.append(f.readline().strip().split()[-5])
				vals.append(f.readline().strip().split()[-5])

				vals = "\t".join(vals)

				# Print the values to the "01out-bioproject_alignments"
				print(output_directory, vals, sep='\t', file=bpjf)





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
	default='cram',
	type=click.Choice(['cram', 'bam']),
	help="Compression algorithm used for resulting alignment. Cram is more space efficient, but Bam is more robust/portable.")




def align(**params):
	'''Wrapper for alignment using ShortStack/bowtie.'''


	pprint(params)

	rc = requirementClass()
	rc.add_samtools()
	rc.add_bowtie()
	rc.add_shortstack()
	rc.add_rnafold()
	rc.check()

	ic = inputClass(params)
	ic.check(['trimmed_libraries', 'genome_file'])


	output_directory        = str(ic.output_directory)
	trimmed_libraries       = ic.inputs['trimmed_libraries']
	genome_file             = ic.inputs['genome_file']

	cores                   = params['cores']
	compression             = params['compression']




	Path(output_directory+ "/align/").mkdir(parents=True, exist_ok=True)

	temp_folder = output_directory + "/align/temp/"
	align_folder = output_directory + '/align/'

	if isdir(temp_folder):
		shutil.rmtree(temp_folder)

	trimmed_libraries = [t.replace(" ", "\ ") for t in trimmed_libraries]

	args = ["ShortStack", "--readfile"] + trimmed_libraries + ["--genomefile", genome_file, "--bowtie_cores", cores, "--align_only", "--cram", "--mmap", 'u', "--sort_mem", "200M", "--outdir", temp_folder]

	if compression == 'cram':
		args += ['--cram']


	args = list(map(str, args))

	print(" ".join(args))

	p = Popen(args, stdout=PIPE, stderr=PIPE, encoding=ENCODING)


	error_terms = ['FAILED. Aborting']
	for line in p.stderr:
		print(line.strip())

		for e in error_terms:
			if e in line:
				sys.exit("Error detected in ShortStack run!")
	p.wait()




	files = os.listdir(temp_folder)

	## Copying files to align directory

	for file in files:

		shutil.move(Path(temp_folder, file), Path(align_folder, file))


	shutil.rmtree(temp_folder)


	alignment_file = [f for f in files if f.endswith(compression)][0]
	# print("found:", alignment_file)


	os.rename(Path(align_folder, alignment_file), Path(align_folder, f'alignment.{compression}'))
	os.rename(Path(align_folder, 'Log.txt'), Path(align_folder, 'log.txt'))


	get_log_details(output_directory)

	ic.inputs['alignment_file'] = str(Path(output_directory, 'align', f"alignment.{compression}").absolute())
	ic.write()






