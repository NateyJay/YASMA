
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



from datetime import datetime




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
			try:
				self.depths[r] += 1
			except IndexError:
				print(f"WARNING: position {r:,} is longer than the total chromosome length")
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




@cli.command(group='Calculation', help_priority=3)



@optgroup.group('\n  Required options',
				help='')

@optgroup.option("-a", "--alignment_file", 
	required=False, 
	type=click.UNPROCESSED, callback=validate_path,
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

@optgroup.option("--gene_annotation_file", 
	# default=f"Annotation_{round(time())}", 
	required=False,
	default=None,
	type=click.UNPROCESSED, callback=validate_path,
	help='Annotation file for genes (gff3) matching the included genome.')


@optgroup.option("-g", "--genome_file", 
	# default=f"Annotation_{round(time())}", 
	required=False,
	type=click.UNPROCESSED, callback=validate_path,
	help='Genome or assembly which was used for the original alignment.')


@optgroup.group("\n  Should include")

@optgroup.option("-j", "--jbrowse_directory", 
	# default=f"Annotation_{round(time())}", 
	required=False,
	default=None,
	type=click.UNPROCESSED, callback=validate_path,
	help='A path to a working directory for a jbrowse2 instance.')


@optgroup.group('\n  Options',
				help='')

@optgroup.option("--min_size",
	default=20,
	help="Minimum read size for specific coverage treatment Default 20 nt.")

@optgroup.option("--max_size",
	default=25,
	help="Maximum read size for specific coverage treatment Default 25 nt.")

@optgroup.option("--force",
	is_flag=True,
	default=False,
	help="Force remake of coverage files")

@optgroup.option("--overwrite_config",
	is_flag=True,
	default=False,
	help="Option to overwrite and make a new jbrowse config.json de novo.")

@optgroup.option('-x', '--remove_name', 
	required=False,
	multiple=True,
	help="Names of entries which should be removed from the config. These are synonymous with the output_directories of the runs used. Use this with caution!")


# @optgroup.option('--delete', 
# 	is_flag=True,
# 	default=False,
# 	help="Deletes source data of those files removed")



def jbrowse(**params):
	'''Tool to build coverage and config files for jbrowse2.'''

	ic = inputClass(params)
	ic.check(["alignment_file", "annotation_readgroups", "genome_file"])

	output_directory       = str(ic.output_directory)
	alignment_file         = ic.inputs["alignment_file"]
	annotation_readgroups  = ic.inputs['annotation_readgroups']
	genome_file            = ic.inputs['genome_file']
	jbrowse_directory      = ic.inputs['jbrowse_directory']
	gene_annotation_file   = ic.inputs["gene_annotation_file"]

	min_size               = params['min_size']
	max_size               = params['max_size']
	force                  = params['force']
	overwrite_config       = params['overwrite_config']
	removes                = params['remove_name']

	if jbrowse_directory:
		jbrowse_directory = jbrowse_directory.rstrip("/")
		input_config = f"{jbrowse_directory}/config.json"

		if not isfile(input_config):
			print("Warning: config.json in jbrowse_directory not found")
			print(" ->", input_config)
			print()
			print('  A new config will be generated de novo. Please check to confirm the jbrowse_directory exists and is the location of a jbrowse2 installation.')

			input_config = None

	else:
		input_config = None

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

	def make_annotation_track(name, uri):

		# track_id = f"{output_directory}_{name}"

		d = { "type": "FeatureTrack",
				"trackId": name,
				"name": name, #f"{output_directory}_{name}",
				"adapter": {
					"type": "Gff3Adapter",
					"gffLocation": {
						"uri": uri, #f"{genome_name}/{output_directory}_annotations/{name}.gff3",
						"locationType": "UriLocation"
					}
				},
				"assemblyNames": [
					genome_name
				],
				"displays": [
					{
						"type": "LinearBasicDisplay",
						"displayId": name, #f"{output_directory}_{name}",
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


	def backup_config():
		now = datetime.now()
		backup_dir = "/".join(input_config.split("/")[:-1]) + "/config_backups/" + now.strftime("%Y.%m.%d")
		Path(backup_dir).mkdir(parents=True, exist_ok=True)
		backup_file = backup_dir + "/config_" + str(round(time())) + ".json"
		copyfile(input_config, backup_file)

		print(f"  Backing up config: {backup_file}")


	def read_config(input_config, overwrite):
		print()
		if not input_config or overwrite:
			print("Making config.json de novo")
			return({})
		with open(input_config, 'r') as f:
			print(f"Reading from input config.json: {input_config}")
			data = json.load(f)
			backup_config()
		return(data)

	config_d = read_config(input_config, overwrite_config)








	if 'plugins' not in config_d:
		config_d['plugins'] = []

	names = [c['name'] for c in config_d['plugins']]
	if 'MyNoBuildPlugin' not in names:
		print(f"  adding plugins: {color.BOLD}MyNoBuildPlugin{color.END}")
		config_d['plugins'].append(make_plugins_object())
	


	## Assemblies

	if 'assemblies' not in config_d:
		config_d['assemblies'] = []

	names = [c['name'] for c in config_d['assemblies']]
	if genome_name not in names:
		print(f"  adding assembly: {color.BOLD}{genome_name}{color.END}")
		config_d['assemblies'].append(make_assembly_object())
		


	#### Tracks

	if 'tracks' not in config_d:
		config_d['tracks'] = []

	names = [c['name'] for c in config_d['tracks']]

	## Gene annotation
	if gene_annotation_file:
		name = gene_annotation_file.rstrip('.gff3').rstrip(".gff").split("/")[-1]
		if name not in names:
			print(f"  adding track: {color.BOLD}{name}.gff3{color.END}")
			at = make_annotation_track(
				name=name,
				uri=f"{genome_name}/{name}.gff3")
			config_d['tracks'].append(at)

	## sRNA loci
	if f"{output_directory}_Loci" not in names:
		print(f"  adding track: {color.BOLD}{output_directory}_Loci{color.END}")
		at = make_annotation_track(
			name=f'{output_directory}_Loci', 
			uri=f"{genome_name}/{output_directory}_annotations/loci.gff3")
		config_d['tracks'].append(at)

	## sRNA regions
	if f"{output_directory}_Regions" not in names:
		print(f"  adding track: {color.BOLD}{output_directory}_Regions{color.END}")
		at = make_annotation_track(
			name=f'{output_directory}_Regions', 
			uri=f"{genome_name}/{output_directory}_annotations/regions.gff3")
		config_d['tracks'].append(at)

	## Coverages
	if f"{output_directory}_Coverage" not in names:
		print(f"  adding track: {color.BOLD}{output_directory}_Coverage{color.END}")
		config_d['tracks'].append(make_bigwig_track())



	if len(removes) > 0:
		print()
		print("Deleting entries specified by option -x:")
		print(" ", removes)
		print()
		for name in removes:

			terms = [name + t for t in ['_Loci','_Regions','_Coverage']]
			tracks_to_remove = []

			# print(terms)

			for i,track in enumerate(config_d['tracks']):
				if track['name'] in terms:
					print(f"  {color.BOLD}x{color.END}", track['name'])
					tracks_to_remove.append(i)


			for i in tracks_to_remove[::-1]:
				del config_d['tracks'][i]
			# print(track[name + '_Coverage'])

		# pprint(type(config_d['tracks']))
		# pprint(len(config_d['tracks']))

		# print(name)


	config_file = f'{output_directory}/jbrowse/config.json'
	with open(config_file, 'w') as outf:
		outf.write(json.dumps(config_d, indent=2))


	# sys.exit()

	## Copying key files to jbrowse folder


	print()
	print("Copying data to jbrowse folder...")

	def copy_it(src, des):
		if force or not isfile(des):
			print(f"  {src.split('/')[-1]}")
			copyfile(src, des)

	copy_it(f"{output_directory}/peak/loci.gff3", f"{ann_dir}/loci.gff3")
	copy_it(f"{output_directory}/peak/regions.gff3", f"{ann_dir}/regions.gff3")
	copy_it(genome_file, f"{output_directory}/jbrowse/{genome_name}/{genome_name}.fa")
	copy_it(genome_file+'.fai', f"{output_directory}/jbrowse/{genome_name}/{genome_name}.fa.fai")
	copy_it(gene_annotation_file, f"{output_directory}/jbrowse/{genome_name}/{gene_annotation_file.split('/')[-1]}")


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
		print("Coverage files found -> skipping make! (override with --force)")

	else:

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
			print(f"chrom: {chrom} -> {chrom_length:,} bp")

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



			print()
			print("   writing: ", end='')
			for key in keys:
				print(key, end='  ', flush=True)
				bw_d[key].rle(chrom)

			print()
			print()

		for key in keys:
			bw_d[key].bw.close()



	if not isdir(jbrowse_directory):
		print("Error: jbrowse_directory not found")
		print(" ->", jbrowse_directory)
		print()
		print("Cannot copy files! Please confirm this directory exists, and reset if needed with 'yasma.py inputs --jbrowse_directory'")
		print()
		sys.exit()

	else:
		print()
		print('Copying to jbrowse_directory with rsync...')


		p = Popen(['rsync', '-arv', f'{output_directory}/jbrowse/', jbrowse_directory], 
			stdout=PIPE, encoding=ENCODING)

		for l in p.stdout:
			print(" ", l.strip())
		p.wait()















