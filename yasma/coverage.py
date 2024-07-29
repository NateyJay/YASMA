
import sys
import os

import click
from click_option_group import optgroup

from pathlib import Path
from os.path import isfile, isdir
from collections import Counter#, deque
from pprint import pprint


from .generics import *
from .cli import cli



@cli.command(group='Calculation', help_priority=4)


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


def coverage(**params):
	"""Produces bigwig coverage files for use in jbrowse."""

	rc = requirementClass()
	rc.add_samtools()
	rc.check()

	ic = inputClass(params)
	ic.check(['alignment_file'])

	output_directory        = ic.output_directory
	alignment_file          = ic.inputs["alignment_file"]
	project_name            = ic.inputs['project_name']


	chromosomes, bam_rgs = get_chromosomes(alignment_file)

	chrom_depth_c = get_global_depth(alignment_file, aggregate_by=['rg','chrom'])
	aligned_depth = sum(chrom_depth_c.values())

	bw_file = alignment_file.with_suffix(".bw")

	tc = trackClass(bw_file, chromosomes)


	tc.process_bam(alignment_file, aligned_depth)
	tc.close()


