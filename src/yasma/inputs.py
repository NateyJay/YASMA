


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
	required=False,
	type=click.UNPROCESSED, callback=validate_outdir,
	help="Directory name for annotation output. Defaults to the current directory, with this directory name as the project name.")




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


@optgroup.group('\n  File input parameters',
				help='')

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

@optgroup.group('\n  Other shared parameters',
				help='')

@optgroup.option("-s", "--srrs", 
	required=False, 
	multiple=True,
	help='NCBI SRA codes for libraries. These will almost certainly start with SRR or ERR.')

@optgroup.option("-c", "--conditions", 
	required=False, 
	multiple=True,
	type=click.UNPROCESSED, callback=validate_condition,
	help='Values denoting condition groups (sets of replicate libraries) for projects with multiple tissues/treatments/genotypes. Can be entered here as space sparated duplexes, with the library base_name and condition groups delimited by a colon. E.g. SRR1111111:WT SRR1111112:WT SRR1111113:mut SRR1111114:mut')

@optgroup.option('-ac', '--annotation_conditions', 
	required=False,
	multiple=True,
	default=None,
	help="List of conditions names which will be included in the annotation. Defaults to use all libraries, though this is likely not what you want if you have multiple groups.")

@optgroup.option("--min_length",
	default = 15,
	help= 'Minimum allowed size for a trimmed read. (default 10)')


@optgroup.option("--max_length",
	default = 50,
	help= 'Maxiumum allowed size for a trimmed read. (default 50)')

def inputs(**params):
	'''Initialize a project and log inputs for later analyses'''


	rc = requirementClass()
	# rc.add_samtools()
	rc.check()

	ic = inputClass(params)
	ic.check_chromosomes()

	ic.check_paired_end()

	print()
	print("Inputs:")
	pprint(ic.inputs)





























