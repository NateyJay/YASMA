# sRNA prechecking module

import sys
import click
from click_option_group import optgroup

from pathlib import Path
from os.path import isfile, isdir
# from collections import Counter, deque


from .generics import *
from .cli import cli


@cli.command(name='readgroups', group='Utilities', help_priority=7)


@optgroup.group('\n  Basic options',
				help='')

@optgroup.option("-a", "--alignment_file", 
	required=False, 
	type=click.UNPROCESSED, callback=validate_path,
	help='Alignment file input (bam or cram).')

@optgroup.option("-o", "--output_directory", 
	required=True,
	type=click.Path(),
	help="Directory name for annotation output")

@optgroup.option('-ac', '--annotation_conditions', 
	required=False,
	multiple=True,
	default=None,
	help="List of conditions names which will be included in the annotation. Defaults to use all libraries, though this is likely not what you want if you have multiple groups.")


def list_readgroups(** params):
	"""Convenience function to list readgroups in an alignment file."""

	rc = requirementClass()
	rc.check()

	ic = inputClass(params)
	ic.check(['alignment_file'])

	output_directory     = str(ic.output_directory)
	alignment_file       = ic.inputs['alignment_file']
	conditions           = ic.inputs['conditions']

	chromosomes, rgs = get_chromosomes(alignment_file)

	conditions = reverse_conditions(conditions)

	print("Libraries (readgroups) found in file:")

	for rg in rgs:
		try:
			cond = f"\t-> {conditions[rg]}"
		except KeyError:
			cond = ''
		print("  ", rg, cond, sep='')

	missing_libs = []
	for lib, cond in conditions.items():

		if lib not in rgs:
			missing_libs.append(lib)

	if len(missing_libs) > 0:
		print()
		print(f"Warning: {len(missing_libs)} libraries found in conditions that are not present in the alignment")
		for lib in missing_libs:
			print("  " + lib, conditions[lib], sep='\t')





