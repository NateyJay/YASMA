


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

@optgroup.option("-a", "--alignment_file", 
	required=False, 
	type=click.UNPROCESSED, callback=validate_path,
	help='Alignment file input (bam or cram).')

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




def shortstack3(**params):
	'''Wrapper for annotation using ShortStack3.'''

	module_title(__name__, version= 0.01)

	rc = requirementClass()
	rc.add_samtools()
	# rc.add_bowtie()
	# rc.add_shortstack()
	# rc.add_rnafold()
	rc.check()

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


		dir_name = f'shortstack3_{subsample.string}{subsample.seed_string}'
	else:
		dir_name = f'shortstack3'


	if name:
		dir_name = f'shortstack3_{name}'




	Path(output_directory, dir_name).mkdir(parents=True, exist_ok=True)

	log_file = f"{output_directory}/{dir_name}/yasma_log.txt"
	sys.stdout = Logger(log_file)

	temp_folder  = Path(output_directory, dir_name, "temp")
	annotation_folder = Path(output_directory, dir_name)


	if isdir(temp_folder):
		shutil.rmtree(temp_folder)




	args = ["ShortStack3"]

	if alignment_file.suffix == '.bam':
		args += ['--bamfile', alignment_file]
	elif alignment_file.suffix == '.cram':
		args += ['--cramfile', alignment_file]
	else:
		raise "alignment_file suffix not recognized..."

	args += ["--genomefile", genome_file, "--outdir", temp_folder]



	args = list(map(str, args))

	print(" ".join(args))

	p = Popen(args)#, stdout=PIPE, stderr=PIPE, encoding=ENCODING)
	p.wait()


	os.rename(Path(temp_folder, 'log.txt'), Path(temp_folder, 'shortstack_log.txt'))

	for file in os.listdir(temp_folder):
		Path(temp_folder, file).rename(Path(annotation_folder, file))

	shutil.rmtree(temp_folder)

	now = datetime.now()

	date_time = now.strftime("%Y/%m/%d, %H:%M:%S")
	print("Run completed:",date_time)	






