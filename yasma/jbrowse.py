
import sys

import click
from click_option_group import optgroup

from pathlib import Path
from os.path import isfile, isdir
from collections import Counter#, deque
from pprint import pprint
from random import sample

import numpy as np
# from statistics import quantiles
import math

from .generics import *
from .cli import cli

from statistics import mean, median

from time import time

import json

from shutil import copyfile
from subprocess import PIPE, Popen







class bigwigClass():
	'''A class to handle producing rpm bigwig files from a counter object c[pos] = depth'''

	def __init__(self, file, total_reads, chromosomes, strand= "+", name=''):
		self.file = file
		# self.file = f"./{output_directory}/Coverages/{file}.wig"
		self.bw = pyBigWig.open(self.file, 'w')

		self.total_reads   = total_reads

		self.bw.addHeader(chromosomes)

		if strand == "+":
			self.strand = 1
		elif strand == "-":
			self.strand = -1 

		self.name= name


	def reset(self, chrom_length):
		# self.last_depth_pos = 1
		self.last_pos = 1
		# self.last_depth = 0
		# self.span = 0
		self.depths = [0] * chrom_length


	def add(self, pos, length):

		for r in range(pos, pos+length+1):
			self.depths[r] += 1

		self.last_pos = r

	def rle(self, chrom):
		last_depth = 0
		span = 1
		last_pos = 0

		starts = []
		ends   = []
		values = []

		# print(self.depths)
		# sys.exit()

		for pos, depth in enumerate(self.depths):

			if depth == last_depth:
				span += 1

			else:


				ends.append(pos)
				starts.append(last_pos)
				values.append(round(last_depth / self.total_reads * 1000000, 8) * self.strand)

				# print(self.name, chrom, last_pos, "->", pos, depth, sep='\t')

				last_depth = depth
				last_pos   = pos
				span       = 1

		if last_pos < pos:
			ends.append(pos)
			starts.append(last_pos)
			values.append(round(last_depth / self.total_reads * 1000000, 8) * self.strand)


		self.bw.addEntries([chrom] * len(values), starts, ends=ends, values=values)




@cli.command(group='Utilities', help_priority=3)



@optgroup.group('\n  Required options',
				help='')

@optgroup.option("-a", "--alignment_file", 
	required=False, 
	type=click.Path(exists=True),
	help='Alignment file input (bam or cram).')

@optgroup.option('-r', '--annotation_readgroups', 
	required=False,
	multiple=True,
	help="List of read groups (RGs, libraries) to be considered for the annotation. 'ALL' uses all readgroups for annotation, but often pertainent RGs will need to be specified individually.")

@optgroup.option("-o", "--output_directory", 
	# default=f"Annotation_{round(time())}", 
	required=True,
	type=click.Path(),
	help="Directory name for annotation output.")


@optgroup.option("-g", "--genome_file", 
	# default=f"Annotation_{round(time())}", 
	required=False,
	type=click.Path(exists=True),
	help='Genome or assembly which was used for the original alignment.')



@optgroup.group('\n  Options',
				help='')

@optgroup.option("-c", "--config_file", 
	# default=f"Annotation_{round(time())}", 
	required=False,
	type=click.Path(exists=True),
	help='An already made config to which we will add our entries.')

@optgroup.option("--min_size",
	default=20,
	help="Minimum read size for specific coverage treatment Default 20 nt.")

@optgroup.option("--max_size",
	default=25,
	help="Maximum read size for specific coverage treatment Default 25 nt.")

@optgroup.option("--force",
	default=False,
	help="Maximum read size for specific coverage treatment Default 25 nt.")



def jbrowse(**params):
	'''Tool to build coverage and config files for jbrowse2.'''

	ic = inputClass(params)
	ic.check(["alignment_file", "annotation_readgroups", "genome_file"])

	output_directory       = ic.inputs['output_directory']
	alignment_file         = ic.inputs["alignment_file"]
	annotation_readgroups  = ic.inputs['annotation_readgroups']
	genome_file            = ic.inputs['genome_file']
	input_config           = ic.inputs['config_file']

	min_size               = params['min_size']
	max_size               = params['max_size']
	force                  = params['force']

	genome_name = genome_file.rstrip(".gz").rstrip(".fa").split("/")[-1]

	cov_dir = output_directory+ f"/jbrowse/{genome_name}/{output_directory}_coverages"
	ann_dir = output_directory+ f"/jbrowse/{genome_name}/{output_directory}_annotations"
	deb_dir = output_directory+ f"/jbrowse/{genome_name}/{output_directory}_debug"

	for directory in [cov_dir, ann_dir, deb_dir]:
		Path(directory).mkdir(parents=True, exist_ok=True)


	chromosomes, bam_rgs = get_chromosomes(alignment_file, output_directory)
	annotation_readgroups = check_rgs(annotation_readgroups, bam_rgs)

	chrom_depth_c = get_global_depth(output_directory, alignment_file, aggregate_by=['rg','chrom'])

	keys = list(chrom_depth_c.keys())
	for key in keys:
		if key[0] in annotation_readgroups:
			chrom_depth_c[key[1]] += chrom_depth_c[key]

		del chrom_depth_c[key]

	total_reads = sum(chrom_depth_c.values())


	color_d = {
		"non" : "grey",
		20 : "lightblue",
		21 : "blue",
		22 : "mediumseagreen",
		23 : "orange",
		24 : "tomato",
		25 : "darkred"
	}

	sizes = [r for r in range(min_size, max_size + 1)] + ["non"]

	for r in sizes:
		if r not in color_d:
			color_d[r] = 'pink'


	def make_counter_dict():
		d = {}

		for r in sizes:
			for s in ['-','+']:

				d[f"{r}{s}"] = Counter()

		return(d)


	def make_plugins_object():
		d = {"name": "MyNoBuildPlugin",
			"umdLoc": {"uri": "myplugin.js"}
			}
			
		
		return(d)


	def make_assembly_object():
		d = {
				"name": genome_name,
				"sequence": {
				"type": "ReferenceSequenceTrack",
				"trackId": f"{genome_name}-ReferenceSequenceTrack",
				"adapter": {
						"type": "IndexedFastaAdapter",
						"fastaLocation": {
							"uri": f"{genome_name}/{genome_name}.fa",
							"locationType": "UriLocation"
					},
						"faiLocation": {
							"uri": f"{genome_name}/{genome_name}.fa.fai",
							"locationType": "UriLocation"
						}
					}
				}
			}
		return(d)

	def make_annotation_track(name):

		track_id = f"{output_directory}_{name}"

		d = { "type": "FeatureTrack",
				"trackId": track_id,
				"name": f"{output_directory}_{name}",
				"adapter": {
					"type": "Gff3Adapter",
					"gffLocation": {
						"uri": f"{genome_name}/{output_directory}_annotations/{name}.gff3",
						"locationType": "UriLocation"
					}
				},
				"assemblyNames": [
					genome_name
				],
				"displays": [
					{
						"type": "LinearBasicDisplay",
						"displayId": f"{output_directory}_{name}",
						"renderer": {
							"type": "SvgFeatureRenderer",
							"color1": "jexl:customColor(feature)"
						}
					}
				]
			}
		return(d)


	def make_bigwig_track():

		track_id = f"{output_directory}_Coverage"

		d = {
			"type": "MultiQuantitativeTrack",
			"trackId": track_id,
			"name": f"{output_directory}_Coverage",
			"assemblyNames": [
				genome_name
				],
			"adapter": {
				"type": "MultiWiggleAdapter",
				"subadapters": []
				},
			"displays": 
				[
					{
						"type": "MultiLinearWiggleDisplay",
						"displayId": f"{output_directory}_coverage",
						"defaultRendering": "multiline",
						"renderers": {
							"MultiXYPlotRenderer": {
								"type": "MultiXYPlotRenderer",
								"bicolorPivot": "none",
								"filled": False
								},
							"MultiRowLineRenderer": {
								"type": "MultiRowLineRenderer",
								"color": "rgb(202, 178, 202)"
								}
						}
					}
				]
			}



		for size in sizes:
			for strand in ['-', '+']:
				key = str(size) + strand

				subadapter = {
								"name": key,
								"type": "BigWigAdapter",
								"bigWigLocation": {
									"name": f"{key}.bigwig",
									"locationType": "UriLocation",
									"uri": f"{genome_name}/{output_directory}_coverages/{key}.bigwig"
									},
								"color": color_d[size]
							}
				d['adapter']['subadapters'].append(subadapter)

		return(d)



	def read_config(input_config):
		print()
		if not input_config:
			print("Making config.json de novo")
			return({})
		with open(input_config, 'r') as f:
			print(f"reading from input config.json: {input_config}")
			data = json.load(f)
		return(data)

	config_d = read_config(input_config)

	if 'plugins' not in config_d:
		config_d['plugins'] = []

	names = [c['name'] for c in config_d['plugins']]
	if 'MyNoBuildPlugin' not in names:
		print(f"  adding plugins: {color.BOLD}MyNoBuildPlugin{color.END}")
		config_d['plugins'].append(make_plugins_object())
	

	if 'assemblies' not in config_d:
		config_d['assemblies'] = []

	names = [c['name'] for c in config_d['assemblies']]
	if genome_name not in names:
		print(f"  adding assembly: {color.BOLD}{genome_name}{color.END}")
		config_d['assemblies'].append(make_assembly_object())
		

	if 'tracks' not in config_d:
		config_d['tracks'] = []

	names = [c['name'] for c in config_d['tracks']]

	if f"{output_directory}_Loci" not in names:
		print(f"  adding track: {color.BOLD}{output_directory}_Loci{color.END}")
		config_d['tracks'].append(make_annotation_track('Loci'))

	if f"{output_directory}_Regions" not in names:
		print(f"  adding track: {color.BOLD}{output_directory}_Regions{color.END}")
		config_d['tracks'].append(make_annotation_track('Regions'))

	if f"{output_directory}_Coverage" not in names:
		print(f"  adding track: {color.BOLD}{output_directory}_Coverage{color.END}")
		config_d['tracks'].append(make_bigwig_track())


	config_file = f'{output_directory}/jbrowse/config.json'
	with open(config_file, 'w') as outf:
		outf.write(json.dumps(config_d, indent=2))




	## Copying key files to jbrowse folder


	print()
	print("Copying data to jbrowse folder...")

	def copy_it(src, des):
		if force or not isfile(des):
			print(f"  {src.split('/')[-1]}")
			copyfile(src, des)

	copy_it(f"{output_directory}/Loci.gff3", f"{ann_dir}/Loci.gff3")
	copy_it(f"{output_directory}/peak/Regions.gff3", f"{ann_dir}/Regions.gff3")
	copy_it(genome_file, f"{output_directory}/jbrowse/{genome_name}/{genome_name}.fa")
	copy_it(genome_file+'.fai', f"{output_directory}/jbrowse/{genome_name}/{genome_name}.fa.fai")


	print()
	print(f"Calculating coverage for sizes {min_size} to {max_size}...")
	print()

	keys = []
	for r in sizes:
		for s in ['-','+']:
			key = f"{r}{s}"
			keys.append(key)




	coverage_already_made = False
	for key in keys:
		file = f"{cov_dir}/{key}.bigwig"
		if isfile(file):
			coverage_already_made = True

	if coverage_already_made and not force:
		sys.exit("Coverage files found -> skipping make! (override with --force)")



	bw_d = {}
	for key in keys:
		bw_d[key] = bigwigClass(f"{cov_dir}/{key}.bigwig", 
				total_reads=total_reads, 
				chromosomes=chromosomes,
				strand = key[-1],
				name=key)


	counter_d = {}

	for chrom_count, chrom_and_length in enumerate(chromosomes):



		chrom, chrom_length = chrom_and_length
		print(f"{chrom_count+1} / {len(chromosomes)}")
		print(f"chrom: {chrom} - {chrom_length:,} bp")

		print()


		for key in keys:
			bw_d[key].reset(chrom_length)



		perc = percentageClass(1, chrom_depth_c[chrom])

		reads = samtools_view(alignment_file, rgs=annotation_readgroups, locus=chrom)

		for i, read in enumerate(reads):

			strand, length, _, pos, _, _, _, _ = read
			# strand, length, size, sam_pos, sam_chrom, rg, seq, read_id)

			if min_size <= length <= max_size:
				key = str(length) + strand
			else:
				key = "non" + strand


			bw_d[key].add(pos, length)


			perc_out = perc.get_percent(i)
			if perc_out:
				print(f"   reading position depths ..... {perc_out}%", end='\r', flush=True)


			# bw_d[key].add(counter_d[key], pos, chrom)

		print()
		print("   writing: ", end='')
		for key in keys:
			print(key, end='  ', flush=True)
			bw_d[key].rle(chrom)

		print()
		print()

	for key in keys:
		bw_d[key].bw.close()









