


import sys

import click
from click_option_group import optgroup

from pathlib import Path
from os.path import isfile, isdir
from collections import Counter#, deque
from pprint import pprint


from .generics import *
from .cli import cli

from shutil import copyfile






@cli.command(group='Preliminary', help_priority=0)



@optgroup.group('\n  Required',
				help='')


@optgroup.option("-o", "--output_directory", 
	# default=f"Annotation_{round(time())}", 
	required=True,
	type=click.Path(),
	help="Directory name for annotation output.")




@optgroup.group('\n  Inputs which may be logged',
				help='')

@optgroup.option("-a", "--alignment_file", 
	required=False, 
	default=None,
	type=click.UNPROCESSED, callback=validate_path,
	help='Alignment file input (bam or cram).')

# @optgroup.option('-r', '--annotation_readgroups', 
# 	required=False,
# 	default=None,
# 	multiple=True,
# 	# type=list,
# 	help="List of read groups (RGs, libraries) to be considered for the annotation. 'ALL' uses all readgroups for annotation, but often pertainent RGs will need to be specified individually.")



@optgroup.option("-g", "--genome_file", 
	# default=f"Annotation_{round(time())}", 
	required=False,
	default=None,
	type=click.UNPROCESSED, callback=validate_path,
	help='Genome or assembly which was used for the original alignment.')


@optgroup.option("-j", "--jbrowse_directory", 
	# default=f"Annotation_{round(time())}", 
	required=False,
	default=None,
	type=click.UNPROCESSED, callback=validate_path,
	help='A path to a working directory for a jbrowse2 instance.')



@optgroup.option("-ga", "--gene_annotation_file", 
	# default=f"Annotation_{round(time())}", 
	required=False,
	default=None,
	type=click.UNPROCESSED, callback=validate_path,
	help='Annotation file for genes (gff3) matching the included genome.')



@optgroup.option('-tl', "--trimmed_libraries", 
	required=False, 
	type=click.UNPROCESSED, callback=validate_glob_path,
	multiple=True,
	help='Path to trimmed libraries. Accepts wildcards (*).')


@optgroup.option("-ul", "--untrimmed_libraries", 
	required=False, 
	type=click.UNPROCESSED, callback=validate_glob_path,
	multiple=True,
	help='Path to untrimmed libraries. Accepts wildcards (*).')

@optgroup.option("-s", "--srrs", 
	required=False, 
	multiple=True,
	help='NCBI SRA codes for libraries. These will almost certainly start with SRR or ERR.')


@optgroup.option("-r", "--replicate_groups", 
	required=False, 
	multiple=True,
	type=click.UNPROCESSED, callback=validate_rep_group,
	help='Values denoting replicate groups (sets of replicate libraries) for projects with multiple conditions/tissues/treatments. Can be entered here as space sparated duplexes, with the library base_name and replicate_group delimited by a colon. E.g. SRR8280355:WT SRR8280356:WT SRR8280357:mut SRR8280358:mut')

def inputs(**params):
	'''A tool to log inputs, which will be referenced by later tools.'''


	rc = requirementClass()
	rc.add_samtools()
	rc.check()

	ic = inputClass(params)
	ic.check_chromosomes()

	ic.check_paired_end()

	print()
	print(color.BOLD + "Inputs:" + color.END)
	pprint(ic.inputs)





























