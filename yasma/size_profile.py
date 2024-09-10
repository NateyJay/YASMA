import sys
import os

import click
from click_option_group import optgroup

from pathlib import Path, PurePath
from os.path import isfile, isdir
from collections import Counter#, deque
from pprint import pprint

from statistics import stdev, median

# from random import sample

# import numpy as np
# from statistics import quantiles
# import math
import shutil
# import re

from .generics import *
from .cli import cli





class peakClass():
	def __init__(self,
		project, 
		alignment_file, 
		libraries,
		min_size=15, 
		max_size=35, 
		candidate_threshold=1.0, 
		extension_threshold=0.5):

		self.project             = project
		self.alignment_file      = alignment_file
		self.sizes               = list(range(min_size,max_size))
		self.candidate_threshold = candidate_threshold
		self.extension_threshold = extension_threshold
		self.bam_rgs             = libraries

		self.master = dict()
		self.master['sizes'] = self.sizes

		self.calc_proportions()
		self.calc_statistics()
		self.call_peaks()

	def calc_proportions(self):

		rg_size_c = get_global_depth(self.alignment_file, aggregate_by=['rg','length'])
		rg_c = get_global_depth(self.alignment_file, aggregate_by=['rg'])

		# props = list()  ## a list of all proportions in order
		# prop_d = dict() ## a dictionary of proportions by size
		for rg in self.bam_rgs:
			self.master[rg] = list()

		self.master['prop'] = list()


		for size in self.master['sizes']:
			self.master['prop'].append(0)
			for rg in self.bam_rgs:
				count = rg_size_c[(rg, str(size))]
				prop  = count / rg_c[rg]
				# print(rg, size, count, prop, sep='\t')

				self.master[rg].append(prop)
				self.master['prop'][-1] += prop / len(rg_c)

	def calc_statistics(self):


		props = self.master['prop']

		non_zero_props = [p for p in props if p >= 0]

		sd  = stdev(non_zero_props)
		med = median(non_zero_props)

		self.sd  = sd
		self.med = med


		print()
		print("Basic stats:")
		print()
		print(f"  sd:  {round(sd,4)}")
		print(f"  med: {round(med,4)}")
		print()
		print(f"  zmed = (p - {round(med,4)}) / {round(sd,4)}")
		print()

		zprops = [(p - med) / sd for p in props]

		self.master['zprop'] = list()
		self.master['candidate'] = list()
		self.master['extension'] = list()
		self.master['peak'] = list()

		for z in zprops:
			self.master['zprop'].append(z)
			self.master['candidate'].append(z > self.candidate_threshold)
			self.master['extension'].append(z > self.extension_threshold)
			self.master['peak'].append(None)

	def call_peaks(self):
		props  = self.master['prop']
		zprops = self.master['zprop']


		candidates = [i for i,z in sorted(enumerate(zprops), key=lambda x:x[1], reverse=True) if z > self.candidate_threshold]
		extensions = [i for i,z in sorted(enumerate(zprops), key=lambda x:x[1], reverse=True) if z > self.extension_threshold]

		change_threshold = -50
		max_threshold = 50

		print()
		print("Peak finding log:")
		print()


		peak_i = 0
		for c in candidates:


			if self.master['peak'][c] is None:
				peak_i += 1

				self.master['peak'][c] = peak_i

				print()
				print('peak_i:', peak_i)
				print('candidate:',c)


				for direction in [" ==>", "<== "]:

					p_last = props[c]
					p_cand = props[c]

					if direction == " ==>":
						rang = range(c+1, len(props))
					else:
						rang = range(c-1, -1, -1)

					for r in rang:
						if r not in extensions:
							print(f'({r}) not a candidate')
							break
						# print(r)
						p_curr = props[r]
						p_change = round((p_curr - p_last) / p_last  * 100,1)
						p_max    = round(p_curr / p_cand * 100, 1)


						p_last = p_curr
						if p_max < max_threshold and p_change > change_threshold:
							print(f"({r}) below {max_threshold}% of max peak and above {change_threshold} change_threshold")
							break
						if p_change > 0:
							print(f'({r}) peak growing')
							break

						if self.master['peak'][r]:
							print(f"({r}) pos in peak")
							break

						# print(direction, r, round(p_curr,4), round(p_last,4), p_change, p_max, sep='\t')
						self.master['peak'][r] = peak_i

						if r == 0:
							for i,p in enumerate(self.master['peak']):
								if p == peak_i:
									self.master['peak'][i] = None

							peak_i -= 1

	def peak_table(self, out_file=False):

		props  = self.master['prop']
		zprops = self.master['zprop']
		peaks  = self.master['peak']

		candidates = [i for i,c in enumerate(self.master['candidate']) if c] 
		extensions = [i for i,e in enumerate(self.master['candidate']) if e] 



		


		if out_file:
			outf = open(out_file, 'w')
			print("project\ti\tsize\tprop\tzero\tzmed\tcand\thyst\tpeak",file=outf)


		print()
		print("Sizes in terms of peaks:")
		print()
		print("i\tsize\tprop\tzero\tzmed\tcand\thyst\tpeak")
		print("============================================================")
		for i,s in enumerate(self.sizes):

			print(i, s, round(props[i],4), props[i] ==  0, round(zprops[i],4), i in candidates, i in extensions, peaks[i], sep='\t')

			if out_file:
				print(self.project, i, s, round(props[i],4), props[i] ==  0, round(zprops[i],4), i in candidates, i in extensions, peaks[i], sep='\t', file=outf)


		if out_file:
			outf.close()


		# with open(alignment_file.with_suffix(".prop_summary.txt"), 'w') as outf:


		# 	print("i\tsize\tprop\tzero\tzmed\tcand\thyst\tpeak", file=outf)
		# 	for i,s in enumerate(sizes):
		# 		print(i, s, round(props[i],4), props[i] ==  0, round(zprops[i],4), i in candidates, i in extensions, peaks[i], sep='\t', file=outf)

	def summarize_peaks(self, out_file=None):

		sizes  = self.sizes
		props  = self.master['prop']
		zprops = self.master['zprop']
		peaks  = self.master['peak']
		candidates = [i for i,c in enumerate(self.master['candidate']) if c] 
		extensions = [i for i,e in enumerate(self.master['candidate']) if e] 

		max_peak = max([p for p in peaks if p])

		unplaced = 1.0
		unplaced_count = len(self.sizes)
		print()
		print("Peaks found:")
		print("")

		final_peaks = dict()

		peak_i_name = 1

		if out_file:
			outf = open(out_file, 'w')

		# print('peak','sizes','center','width','prop', 'm_prop', sep='\t', file=outf)
		print('peak','sizes','center','width','prop', 'avg_prop', sep='\t')

		if out_file:
			print('project','peak','sizes','center','width','prop', 'avg_prop', sep='\t', file=outf)

		print("==========================================")
		for peak_i in range(1, max_peak+1):


			positions  = [i for i,p in enumerate(peaks) if p == peak_i]
			peak_sizes = [sizes[p] for p in positions]
			cum_prop   = sum([props[p] for p in positions])
			width      = len(peak_sizes)


			unplaced  -= cum_prop
			unplaced_count -= width

			max_prop   = max([props[p] for p in positions])
			center     = [sizes[p] for p in positions if props[p] == max_prop][0]

			peak_name = f"peak{peak_i_name}"
			peak_i_name += 1

			for s in peak_sizes:
				final_peaks[s] = peak_name


			# print(peak_name, ",".join(map(str,peak_sizes)), center, width, round(cum_prop, 4), round(cum_prop/width, 4), sep='\t', file=outf)

			print(peak_name, ",".join(map(str,peak_sizes)), center, width, round(cum_prop, 4), round(cum_prop/width, 4), sep='\t')
			if out_file:
				print(self.project, peak_name, ",".join(map(str,peak_sizes)), center, width, round(cum_prop, 4), round(cum_prop/width, 4), sep='\t', file=outf)


		# print("none", '-','-',unplaced_count, round(unplaced,4), round(unplaced/unplaced_count,4), sep='\t', file=outf)
		print("none", '-','-',unplaced_count, round(unplaced,4), round(unplaced/unplaced_count,4), sep='\t')

		if (out_file):
			print(self.project, "none", '-','-',unplaced_count, round(unplaced,4), round(unplaced/unplaced_count,4), sep='\t', file=outf)
			outf.close()


	def plot_proportions(self, out_file = None):

		sizes  = self.sizes
		props  = self.master['prop']
		zprops = self.master['zprop']

		if out_file:
			outf = open(out_file, 'w')

		print()
		print("Plot:")
		print()
		print("size\tprop\tzmed\tpeak\t  0    5    10   15   20   25   30   35   40")
		print("==============================    |    |    |    |    |    |    |    |    |")

		print("size\tprop\tzmed\tpeak\t  0    5    10   15   20   25   30   35   40", file=outf)
		print("==============================    |    |    |    |    |    |    |    |    |", file=outf)

		for i,s in enumerate(sizes):
			# z = zprops[i]
			p = props[i]


			val = 0

			bar_string = '  '
			pch = "-"
			while True:
				if val > p:
					break

				z = (val - self.med) / self.sd

				if z > self.extension_threshold:
					pch = "•"
				if z > self.candidate_threshold:
					pch = "*"


				bar_string += pch
				val += 0.01

			peak_i = self.master['peak'][i]
			if peak_i:
				peak_name = f"peak{peak_i}"
			else:
				peak_name = ''

			print(s, round(p,3), round(zprops[i],3), peak_name, bar_string, sep='\t')
			print(s, round(p,3), round(zprops[i],3), peak_name, bar_string, sep='\t', file=outf)



		if out_file:
			outf.close()




@cli.command(group='Utilities', help_priority=5)



@optgroup.group('\n  Required',
				help='')


@optgroup.option("-o", "--output_directory", 
	# default=f"Annotation_{round(time())}", 
	required=True,
	type=click.Path(),
	help="Directory name for annotation output.")

@optgroup.option("-a", "--alignment_file", 
	required=False, 
	default=None,
	type=click.UNPROCESSED, callback=validate_path,
	help='Alignment file input (bam or cram).')

@optgroup.option("-c", "--conditions", 
	required=False, 
	multiple=True,
	type=click.UNPROCESSED, callback=validate_condition,
	help='Values denoting condition groups (sets of replicate libraries) for projects with multiple tissues/treatments/genotypes. Can be entered here as space sparated duplexes, with the library base_name and condition groups delimited by a colon. E.g. SRR1111111:WT SRR1111112:WT SRR1111113:mut SRR1111114:mut')

@optgroup.option('-ac', '--annotation_conditions', 
	required=False,
	multiple=True,
	default=None,
	help="List of conditions names which will be included in the profile. Defaults to use all libraries, though this is likely not what you want if you have multiple groups.")



@optgroup.group('\n  Optional',
				help='')
@optgroup.option('-f','--force', is_flag=True, default=False, help='Flag to force override of output_file (default does nothing when this file is found).')



def size_profile(**params):
	'''Convenience function for calculating aligned size profile.'''


	rc = requirementClass()
	# rc.add_samtools()
	rc.check()


	ic = inputClass(params)

	alignment_file          = ic.inputs['alignment_file']
	output_directory        = ic.output_directory
	conditions              = ic.inputs['conditions']
	annotation_conditions   = ic.inputs['annotation_conditions']


	chromosomes, bam_rgs = get_chromosomes(alignment_file)

	libraries = []

	if len(annotation_conditions) == 0:
		libraries = bam_rgs
	else:
		for a in annotation_conditions:
			libraries += conditions[a]

	pc = peakClass(ic.inputs['project_name'], alignment_file, libraries=libraries)

	pc.peak_table(alignment_file.with_suffix(".peak_table.txt"))
	pc.summarize_peaks(alignment_file.with_suffix(".peak_summary.txt"))
	pc.plot_proportions(alignment_file.with_suffix(".peak_plot.txt"))












