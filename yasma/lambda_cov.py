# Test utility to examine lambda values for poisson distribution

import sys
import click
from pathlib import Path
from os.path import isfile, isdir
from collections import Counter#, deque
from pprint import pprint

import numpy as np
# from statistics import quantiles
import math

from .generics import *
from .cli import cli



@cli.command(group='Testing', help_priority=4)

@click.option("-a", "--alignment_file", 
	required=True, 
	type=click.Path(exists=True),
	help='Alignment file input (bam or cram).')


@click.option("-o", "--output_directory", 
	default=f"Annotation_{round(time())}", 
	type=click.Path(),
	help="Directory name for annotation output.")

@click.option("--window",
	default=40,
	help="Window size (centered on position) for counting reads, used as k in the poisson model. Default 40 nt.")


def lambda_cov(alignment_file, output_directory, window):
	'''Test utility for interrogating lamdba methods for poisson.'''

	output_directory = output_directory.rstrip("/")

	Path(output_directory).mkdir(parents=True, exist_ok=True)


	chromosomes, bam_rgs = get_chromosomes(alignment_file, output_directory)


	lambda_d = get_lambda(alignment_file, chromosomes, window, output_directory)

	pprint(lambda_d)