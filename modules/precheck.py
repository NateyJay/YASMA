# sRNA prechecking module

import sys
import click
from pathlib import Path
from os.path import isfile, isdir
from collections import Counter, deque


from modules.generics import *
from modules.cli import cli


@cli.command(group='Utilities', help_priority=6)

@click.option("-a", "--alignment_file", 
	required=True, 
	type=click.Path(exists=True),
	help='Alignment file input (bam or cram).')

@click.option('-r', '--annotation_readgroups', 
	required=True,
	multiple=True,
	help="List of read groups (RGs, libraries) to be considered for the annotation. 'ALL' uses all readgroups for annotation, but often pertainent RGs will need to be specified individually.")

@click.option("-o", "--output_directory", 
	default=f"Annotation_{round(time())}", 
	type=click.Path(),
	help="Directory name for annotation output")

@click.option("-f", "--force",
	is_flag=True,
	help='Force remake of supporting files')

def precheck(alignment_file, annotation_readgroups, output_directory, force):
	"""Runs precheck to identify likely dicer sizes."""

	Path(output_directory).mkdir(parents=True, exist_ok=True)


	log_file = f"{output_directory}/Log_precheck.txt"

	message = f"log_file is already exists ({log_file}). The annotator will not over-write by default (use --force to override). Be warned: this will trigger the overwrite of some files in this folder!"
	assert not isfile(log_file) or force, message

	sys.stdout = Logger(log_file)


	chromosomes, bam_rgs = get_chromosomes(alignment_file, output_directory)
	annotation_readgroups = check_rgs(annotation_readgroups, bam_rgs)



	depth_c = get_rg_depth(output_directory, alignment_file)


	len_c = Counter()
	for key in depth_c.keys():
		if key[0] in annotation_readgroups:
			len_c[key[1]] += depth_c[key]

	print(f"annotation_readgroups: {annotation_readgroups}")

	print()
	print("length\tp_depth\tp_highest\tabundance")

	dep = sum(len_c.values())
	highest = len_c.most_common(1)[0][1]

	for r in range(15, 31):
		print(r, 
			round(len_c[r] / dep, 4),
			round(len_c[r] / highest, 4),
			len_c[r],
			sep='\t')

	# pprint(depth_c)








