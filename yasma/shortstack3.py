


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





@cli.command(group='Wrappers', help_priority=5)



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




def shortstack3(**params):
	'''Wrapper for alignment using ShortStack/bowtie.'''


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

	cores                   = params['cores']
	compression             = params['compression']



	Path(output_directory+ "/shortstack3/").mkdir(parents=True, exist_ok=True)

	temp_folder = Path(output_directory, "shortstack3/temp/")
	align_folder = Path(output_directory, 'shortstack3/')


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


	for file in os.listdir(temp_folder):
		Path(temp_folder, file).rename(Path(align_folder, file))

	shutil.rmtree(temp_folder)






