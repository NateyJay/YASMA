
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


@optgroup.group('\n  Run options',
				help='')

@optgroup.option("-p", "--peaks", 
	required=True,
	type=str,
	multiple=True,
	help='Entry of size limits for a peak. Encoded as two integers separated by a dash (`-`), for example: 21-22 is an appropriate entry for plant miRNAS. Also accepts single-size peaks without dash. Multiple peaks may be identified, separated with spaces or by calling the option again.')



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
	annotation_readgroups   = ic.inputs['annotation_readgroups']

	peaks = list(params['peaks'])



	chromosomes, bam_rgs = get_chromosomes(alignment_file)
	annotation_readgroups = check_rgs(annotation_readgroups, bam_rgs)

	chrom_depth_c = get_global_depth(alignment_file, aggregate_by=['rg','chrom'])


	keys = list(chrom_depth_c.keys())
	for key in keys:
		if key[0] in annotation_readgroups:
			chrom_depth_c[key[1]] += chrom_depth_c[key]

		del chrom_depth_c[key]



	aligned_depth = sum(chrom_depth_c.values())

	cov_dir = Path(output_directory, 'coverage')
	cov_dir.mkdir(parents=True, exist_ok=True)


	peak_lookup = {}
	for i,peak in enumerate(peaks):
		peak = peak.split("-")
		peak = [int(p) for p in peak]

		if len(peak) == 1:
			peak_lookup[peak[0]] = i
		else:
			for p in range(min(peak), max(peak)+1):
				peak_lookup[p] = i




	bw_d = {}
	bw_d['all'] = bigwigClass(Path(cov_dir, 'all.bw'), aligned_depth, chromosomes, strand= "+", name='all')

	strands = ['+','-']
	for peak in ['all', 'other'] + peaks:
		for strand in strands:
			name = f"{peak}{strand}"
			bw_d[name] = bigwigClass(Path(cov_dir, f'{name}.bw'), aligned_depth, chromosomes, strand= strand, name=name)


	bamf = pysam.AlignmentFile(alignment_file)
	
	for chrom_count, chrom_and_length in enumerate(chromosomes):



		chrom, chrom_length = chrom_and_length
		print(f"{chrom_count+1} / {len(chromosomes)}")
		print(f"chrom: {chrom} -> {chrom_length:,} bp")



		for key in bw_d.keys():
			bw_d[key].reset(chrom_length)

		perc = percentageClass(1, chrom_depth_c[chrom])

		for i,read in enumerate(bamf.fetch(contig=chrom)):

			perc_out = perc.get_percent(i)
			if perc_out:
				print(f"   reading position depths ..... {perc_out}%", end='\r', flush=True)
			
			if read.is_unmapped:
				continue

			if read.get_tag("RG") not in annotation_readgroups:
				continue

			if read.is_forward:
				strand = "+"
			else:
				strand = "-"

			length   = read.query_length
			position = read.reference_start



			try:
				peak_name = f"{peaks[peak_lookup[length]]}{strand}"
			except KeyError:
				peak_name = f'other{strand}'


			bw_d['all'].add(position, length)
			bw_d[f'all{strand}'].add(position, length)
			bw_d[peak_name].add(position, length)


		print()
		for key in bw_d.keys():
			print(f"[{key}]", end='  ', flush=True)
			bw_d[key].rle(chrom)

		print()
		print()


	for key in bw_d.keys():
		bw_d[key].close()

	bamf.close()

	sys.exit()


	# bw_file = alignment_file.with_suffix(".bw")

	tc = trackClass(bw_file, chromosomes)


	tc.process_bam(alignment_file, aligned_depth)
	tc.close()


