


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




def shortstack4(**params):
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



	Path(output_directory+ "/shortstack4/").mkdir(parents=True, exist_ok=True)

	temp_folder = Path(output_directory, "shortstack4/temp/")
	align_folder = Path(output_directory, 'shortstack4/')


	if isdir(temp_folder):
		shutil.rmtree(temp_folder)

	if alignment_file.suffix == '.cram':
		print('.cram file provided... converting to .bam...')

		new_alignment_file = alignment_file.with_suffix(".bam")

		if not isfile(new_alignment_file):
			with open(new_alignment_file, 'wb') as outf:
				args = ['samtools', 'view', '-b', '-h', alignment_file]
				p = Popen(args, stdout=outf)
				p.wait()
				
		alignment_file = new_alignment_file


	args = ["ShortStack4", '--bamfile', alignment_file, "--genomefile", genome_file, "--outdir", temp_folder]

	args = list(map(str, args))

	print(" ".join(args))

	p = Popen(args)#, stdout=PIPE, stderr=PIPE, encoding=ENCODING)
	p.wait()


	for file in os.listdir(temp_folder):
		Path(temp_folder, file).rename(Path(align_folder, file))

	shutil.rmtree(temp_folder)






