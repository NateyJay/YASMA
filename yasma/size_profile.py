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

from tqdm import tqdm




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


@optgroup.group('\n  Optional',
				help='')
@optgroup.option('-f','--force', is_flag=True, default=False, help='Flag to force override of output_file (default does nothing when this file is found).')



def size_profile(**params):
	'''Convenience function for calculating aligned size profile.'''


	rc = requirementClass()
	rc.add_samtools()
	rc.check()


	ic = inputClass(params)

	pprint(params)


	output_directory        = ic.output_directory
	alignment_file          = ic.inputs['alignment_file']

	# force = params['force']


	chromosomes, bam_rgs = get_chromosomes(alignment_file)
	rg_size_c = get_global_depth(output_directory, alignment_file, aggregate_by=['rg','length'])
	rg_c = get_global_depth(output_directory, alignment_file, aggregate_by=['rg'])


	

	props = list()  ## a list of all proportions in order
	prop_d = dict() ## a dictionary of proportions by size


	for rg in bam_rgs:
		for i,size in enumerate(range(15, 36)):
			count = rg_size_c[(rg, str(size))]
			prop  = count / rg_c[rg]
			# print(rg, size, count, prop, sep='\t')

			try:
				prop_d[size] += prop / len(rg_c)
			except KeyError:
				prop_d[size] = prop / len(rg_c)

			try:
				props[i] += prop / len(rg_c)
			except IndexError:
				props.append(prop / len(rg_c))




	## showing the averaging of all readgroups
	with open(alignment_file.with_suffix(".sizes.txt"), 'w') as outf:

		print("size", end = '', file=outf)
		for rg in bam_rgs:
			print('\t', rg, sep='', end='', file=outf)
		print('\taverage', file=outf)

		for size in range(15,31):
			print(size, end='\t', file=outf)

			for rg in bam_rgs:

				count = rg_size_c[(rg, str(size))]
				prop  = count / rg_c[rg]

				
				print(round(prop,4), end = '\t', file=outf)

			print(round(prop_d[size],4), file=outf)




	sizes = list(range(15,36))
	peaks = dict()
	for r in range(len(props)):
		peaks[r] = False

	non_zero_props = [p for p in props if p >= 0]

	sd  = stdev(non_zero_props)
	med = median(non_zero_props)

	print()
	print("Basic stats:")
	print()
	print(f"  sd:  {round(sd,4)}")
	print(f"  med: {round(med,4)}")
	print()
	print(f"  zmed = (p - {round(med,4)}) / {round(sd,4)}")
	print()

	zprops = [(p - med) / sd for p in props]
	candidates = [i[0] for i in sorted(enumerate(zprops), key=lambda x:x[1], reverse=True) if i[1] > 1]
	hysterics  = [i[0] for i in sorted(enumerate(zprops), key=lambda x:x[1], reverse=True) if i[1] > 0.5]

	# pprint(zprops)
	# pprint(candidates)


	change_threshold = -50
	max_threshold = 50

	print()
	print("Peak finding log:")
	print()

	peak_i = -1
	for c in candidates:


		if peaks[c] is False:
			peak_i += 1

			peaks[c] = peak_i

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
					if r not in hysterics:
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

					if peaks[r] is not False:
						print(f"({r}) pos in peak")
						break

					print(direction, r, round(p_curr,4), round(p_last,4), p_change, p_max, sep='\t')
					peaks[r] = peak_i



	print()
	print("Sizes in terms of peaks:")
	print()
	print("i\tsize\tprop\tzero\tzmed\tcand\thyst\tpeak")
	print("==========================================================")
	for i,s in enumerate(sizes):
		print(i, s, round(props[i],4), props[i] ==  0, round(zprops[i],4), i in candidates, i in hysterics, peaks[i], sep='\t')



	with open(alignment_file.with_suffix(".peak_prop.txt"), 'w') as outf:


		print("i\tsize\tprop\tzero\tzmed\tcand\thyst\tpeak", file=outf)
		for i,s in enumerate(sizes):
			print(i, s, round(props[i],4), props[i] ==  0, round(zprops[i],4), i in candidates, i in hysterics, peaks[i], sep='\t', file=outf)






	unplaced = 1.0
	print()
	print("Peaks found:")
	print("")

	final_peaks = dict()

	with open(alignment_file.with_suffix(".prop.txt"), 'w') as outf:

		print('peak','sizes','center','width','prop', sep='\t', file=outf)
		print('peak','sizes','center','width','prop', sep='\t')
		print("==========================================")
		for peak_i in range(max(peaks.values())+1):

			positions  = [k for k,v in peaks.items() if str(v) == str(peak_i)]
			peak_sizes = [sizes[p] for p in positions]
			cum_prop   = sum([props[p] for p in positions])
			width      = len(peak_sizes)
			unplaced  -= cum_prop

			max_prop   = max([props[p] for p in positions])
			center     = [sizes[p] for p in positions if props[p] == max_prop][0]

			peak_name = f"peak{peak_i}"

			for s in peak_sizes:
				final_peaks[s] = peak_name


			print(peak_name, ",".join(map(str,peak_sizes)), center, width, round(cum_prop, 4), sep='\t', file=outf)
			print(peak_name, ",".join(map(str,peak_sizes)), center, width, round(cum_prop, 4), sep='\t')

		print("none", '-','-','-', round(unplaced,4), sep='\t', file=outf)
		print("none", '-','-','-', round(unplaced,4), sep='\t')




	print()
	print("Plot:")
	print()
	print("size\tprop\tzmed\tpeak\t  0    5    10   15   20   25   30   35   40")
	print("==============================    |    |    |    |    |    |    |    |    |")

	for i,s in enumerate(sizes):
		# z = zprops[i]
		p = props[i]


		val = 0

		bar_string = '  '
		pch = "-"
		while True:
			if val > p:
				break

			z = (val - med) / sd

			if z > 0.5:
				pch = "•"
			if z > 1:
				pch = "*"


			bar_string += pch
			val += 0.01

		try:
			peak_name = final_peaks[s]
		except KeyError:
			peak_name = ''

		print(s, round(p,3), round(zprops[i],3), peak_name, bar_string, sep='\t')












