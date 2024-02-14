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

from tqdm import tqdm




@cli.command(group='Utilities', help_priority=5)



@optgroup.group('\n  Required',
				help='')


@optgroup.option("-o", "--output_directory", 
	# default=f"Annotation_{round(time())}", 
	required=True,
	type=click.Path(),
	help="Directory name for annotation output.")



@optgroup.group('\n  Optional',
				help='')
@optgroup.option('-f','--force', is_flag=True, default=False, help='Flag to force override of output_file (default does nothing when this file is found).')



def size_profile(**params):
	'''Convenience function for calculating aligned size profile.'''


	rc = requirementClass()
	rc.add_samtools()
	rc.check()


	ic = inputClass(params)



	output_directory        = ic.output_directory
	alignment_file          = ic.inputs['alignment_file']

	force = params['force']


	chromosomes, bam_rgs = get_chromosomes(alignment_file)
	chrom_depth_c = get_global_depth(output_directory, alignment_file, aggregate_by=['rg','chrom'])


	output_file = Path(output_directory, 'align', 'size_profile.txt')

	if not isfile(alignment_file):
		sys.exit(f"Alignment file {alignment_file} not found...")
	if isfile(output_file):
		if not force:
			sys.exit(f"Output file found: {output_file}\nstopping run. (override with -f flag)")
		else:
			print(f'Output file found: {output_file}\n***forcing override due to -f flag.')

	total = sum(chrom_depth_c.values())

	print()
	print(f"Total aligned reads {total:,}")

	print()
	print(f'Reading alignment...')
	pbar = tqdm(total=total)

	c = Counter()
	for line in samtools_view(alignment_file):
		pbar.update()
		strand, length, size, sam_pos, sam_chrom, rg, seq, read_id = line

		c[(length, rg)] += 1



	pbar.close()
	# pprint(c)

	chrom_depth_c = get_global_depth(output_directory, alignment_file, aggregate_by=['rg'])

	with open(output_file, 'w') as outf:
		print("project", "library", "size", 'depth', 'prop', sep='\t', file=outf)
		for rg in bam_rgs:

			for size in range(15, 36):
				depth = c[(size, rg)]
				prop = round(depth / chrom_depth_c[rg], 3)
				print(output_directory.name, rg, size, depth, prop, sep='\t', file=outf)




















