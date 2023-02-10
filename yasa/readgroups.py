# sRNA prechecking module

import sys
import click
from pathlib import Path
from os.path import isfile, isdir
# from collections import Counter, deque


from .generics import *
from .cli import cli


@cli.command(name='readgroups', group='Utilities', help_priority=7)

@click.option("-a", "--alignment_file", 
	required=True, 
	type=click.Path(exists=True),
	help='Alignment file input (bam or cram).')


def list_readgroups(alignment_file):
	"""Convenience function to list readgroups in an alignment file."""

	chromosomes, rgs = get_chromosomes(alignment_file)

	for rg in rgs:
		print(rg)







