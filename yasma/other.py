
import sys
import os

import click
from click_option_group import optgroup

from pathlib import Path, PurePath
from os.path import isfile, isdir
from collections import Counter#, deque
from pprint import pprint

# from random import sample

# import numpy as np
# from statistics import quantiles
# import math
import shutil
# import re

from .generics import *
from .cli import cli





@cli.command(group='Utilities', help_priority=5)



@optgroup.group('\n  Basic options',
				help='')


@optgroup.option("-o", "--output_directory", 
	# default=f"Annotation_{round(time())}", 
	required=True,
	type=click.Path(),
	help="Directory name for annotation output.")




def normalize_alignment_name(**params):
	'''Fixes old alignment file names.'''


	rc = requirementClass()
	rc.add_cutadapt()
	rc.check()


	ic = inputClass(params)



	output_directory        = ic.output_directory
	alignment_file          = ic.inputs['alignment_file']


	if not alignment_file:
		files = os.listdir(Path(output_directory, 'align'))
		if 'merged_alignments.cram' in files:
			alignment_file = Path(output_directory, 'align', 'merged_alignments.cram')
		elif "alignment.cram" in files:
			alignment_file = Path(output_directory, 'align', 'alignment.cram')

		params['alignment_file'] = alignment_file



	compression = alignment_file.suffix
	files = os.listdir(Path(output_directory, 'align'))
	files = [f for f in files if f.endswith(compression)]
	files = [f for f in files if f.startswith('merged_alignments')]

	for file in files:
		file = Path(output_directory, 'align', file)
		stem = file.stem.replace("merged_alignments", 'alignment')
		new_file = file.with_stem(stem)

		if "merged_alignments" in file.stem:
			print(f"alignment_file: {file}")
			print(f"           new: {new_file}")
			print()

			file.rename(new_file)

	alignment_file = Path(output_directory, 'align', 'alignment.cram')
	print(alignment_file)
	# sys.exit()
	ic.inputs['alignment_file'] = alignment_file.absolute()
	ic.write()

	# pprint(ic.inputs)



@cli.command(group='Utilities', help_priority=5)



@optgroup.group('\n  Basic options',
				help='')


@optgroup.option("-o", "--output_directory", 
	# default=f"Annotation_{round(time())}", 
	required=True,
	type=click.Path(),
	help="Directory name for annotation output.")

@optgroup.option("-a", "--alignment_file", 
	required=False, 
	type=click.UNPROCESSED, callback=validate_path,
	help='Alignment file input (bam or cram).')



def cram_to_bam(**params):
	'''Changes crams to bam alignments.'''


	rc = requirementClass()
	rc.add_samtools()
	rc.check()


	ic = inputClass(params)




	output_directory        = ic.output_directory
	alignment_file          = ic.inputs['alignment_file']



	print(alignment_file)
	# print(alignment_file.suffix)

	if alignment_file.suffix == '.cram':
		bam_file = alignment_file.with_suffix('.bam')
		print(f"  -> {bam_file}")
		# print(bam_file)
		# input()

		with open(bam_file, 'wb') as outf:
			call = ['samtools','view','-h','-b',alignment_file]
			p = Popen(call, stdout=outf)
			p.wait()


		os.remove(alignment_file)

		# sys.exit()
		ic.inputs['alignment_file'] = bam_file.absolute()
		ic.write()








