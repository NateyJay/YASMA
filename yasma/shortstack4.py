


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



from datetime import datetime


@cli.command(group='Ann. wrappers', help_priority=5)



@optgroup.group('\n  Basic options',
				help='')


@optgroup.option('-a', "--trimmed_libraries", 
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

@optgroup.option("-n", "--name", 
	# default=f"Annotation_{round(time())}", 
	default=None,
	required=False,
	type=str,
	help="Optional name alignment. Useful if comparing annotations.")


@optgroup.option('--subsample',
	help="Allows the user to subsample alignments for the annotation to a defined depth. Accepts an integer number of reads, which can be modified with a 10^3 prefix (ex. 10M).")

@optgroup.option('--subsample_seed',
	type=int,
	default=42,
	help="Seed value used for subsampling (default: 42)")

@optgroup.option('--override', is_flag=True, default=False, help='Overrides config file changes without prompting.')


def shortstack4(**params):
	'''Wrapper for annotation using ShortStack4.'''


	rc = requirementClass()
	rc.add_samtools()
	# rc.add_bowtie()
	# rc.add_shortstack()
	# rc.add_rnafold()
	rc.check()

	# pprint(params)
	ic = inputClass(params)
	ic.check(['alignment_file', 'genome_file'])


	output_directory        = str(ic.output_directory)
	alignment_file          = ic.inputs['alignment_file']
	genome_file             = ic.inputs['genome_file']

	target_depth            = params['subsample']
	seed                    = params['subsample_seed']
	name                    = params['name']



	chrom_depth_c = get_global_depth(output_directory, alignment_file, aggregate_by=['rg','chrom'])

	aligned_read_count = sum(chrom_depth_c.values())

	if target_depth:
		subsample = parse_subsample(target_depth, alignment_file, "bam", sum(chrom_depth_c.values()), seed=seed)

		perform_subsample(subsample)

		alignment_file = subsample.file


		dir_name = f'shortstack4_{subsample.string}{subsample.seed_string}'
	else:
		dir_name = f'shortstack4'

	if name:
		dir_name = f'shortstack4_{name}'



	temp_folder = Path(output_directory, dir_name, 'temp')
	annotation_folder = Path(output_directory, dir_name)
	annotation_folder.mkdir(parents=True, exist_ok=True)

	log_file = f"{output_directory}/{dir_name}/yasma_log.txt"
	sys.stdout = Logger(log_file)



	if isdir(temp_folder):
		shutil.rmtree(temp_folder)

	if alignment_file.suffix == '.cram':
		print('.cram file provided. Only bam is accepted for SS4.')
		sys.exit()

		print('.cram file provided... converting to .bam...')
		# new_alignment_file = alignment_file.with_suffix(".bam")

		# if not isfile(new_alignment_file):
		# 	with open(new_alignment_file, 'wb') as outf:
		# 		args = ['samtools', 'view', '-b', '-h', alignment_file]
		# 		p = Popen(args, stdout=outf)
		# 		p.wait()

		# # os.remove(alignment_file)
				
		# alignment_file = new_alignment_file



		# ic.inputs['alignment_file'] =  alignment_file
		# ic.write()


	args = ["ShortStack4", '--bamfile', alignment_file, "--genomefile", genome_file, "--outdir", temp_folder , '--threads', '4']

	args = list(map(str, args))

	print(" ".join(args))

	p = Popen(args)#, stdout=PIPE, stderr=PIPE, encoding=ENCODING)
	p.wait()

	# os.rename(Path(temp_folder, 'log.txt'), Path(temp_folder, 'shortstack_log.txt'))


	for file in os.listdir(temp_folder):
		Path(temp_folder, file).rename(Path(annotation_folder, file))

	shutil.rmtree(temp_folder)


	now = datetime.now()


	date_time = now.strftime("%Y/%m/%d, %H:%M:%S")
	print("Run completed:",date_time)	





