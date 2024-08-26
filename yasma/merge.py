
import sys
import os

import click
from click_option_group import optgroup

from pathlib import Path
from collections import Counter
from pprint import pprint

from .generics import *
from .cli import cli




@cli.command(group='Annotation', help_priority=2)


@optgroup.group('\n  Basic options',
				help='')

@optgroup.option('-m', '--merges', 
	required=True,
	multiple=True,
	help="List of condition names to be merged to a single annotation. These simply look for annotations in folders named tradeoff_[condition]. Required, because otherwise there is nothing to merge")

@optgroup.option("-o", "--output_directory", 
	# default=f"Annotation_{round(time())}", 
	required=True,
	type=click.Path(),
	help="Directory name for annotation output.")

@optgroup.option("-n", "--name", 
	# default=f"Annotation_{round(time())}", 
	required=False,
	type=str,
	default=None,
	help="Optional name. This is added to the merge folder name and is used if you want to do multiple different merges.")







def merge(**params):
	'''Tool for merging multiple alignments with bedtools.'''

	rc = requirementClass()
	rc.add_bedtools()
	rc.check()

	ic = inputClass(params)

	output_directory        = ic.output_directory

	name     = params['name']
	merges   = params['merges']


	if name:
		merge_dir = Path(output_directory, f"merge_{name}")
	else:
		merge_dir = Path(output_directory, "merge")


	merge_dir.mkdir(parents=True, exist_ok=True)

	log_file   = Path(merge_dir, "log.txt")
	sys.stdout = Logger(log_file)

	unsorted_file = Path(merge_dir, 'unsorted.gff3')
	sorted_file   = Path(merge_dir, 'sorted.gff3')
	merged_file   = Path(merge_dir, 'merged.bed')

	with open(unsorted_file, 'w') as outf:
		for m in merges:
			file = Path(output_directory, f'tradeoff_{m}', 'loci.gff3')
			print(file)

			with open(file, 'r') as f:
				for line in f:
					outf.write(line)

	command = ['bedtools', 'sort', '-i', str(unsorted_file)]
	print(" ".join(map(str, command)))
	with open(sorted_file, 'w') as outf:
		p = Popen(command, stdout=outf)
		p.wait()


	command = ['bedtools', 'merge', '-i', str(sorted_file), '-s', '-header']
	print(" ".join(map(str, command)))
	with open(merged_file, 'w') as outf:
		p = Popen(command, stdout=outf)
		p.wait()

	sorted_file.unlink()
	unsorted_file.unlink()


	## to add to this analysis:
	## gff3 details, including strand
	## details like length, n-merged, parents, 











