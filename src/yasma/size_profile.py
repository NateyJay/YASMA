

from .generics import *

from statistics import stdev, median






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
		self.libraries           = libraries

		self.min_size = min_size
		self.max_size = max_size

		self.master = dict()
		self.master['sizes'] = self.sizes
		self.master['mask']  = [False] * len(self.sizes)

		self.calc_proportions()
		self.calc_statistics()

		print('masking slopes...')
		mask_i = 0
		while True:
			pre_mask = " ".join(map(str, self.master['mask']))
			self.mask_slope()

			if any(self.master['mask']):
				self.calc_statistics()


			mask_i += 1
			print(f'  iteration {mask_i}')

			if " ".join(map(str, self.master['mask'])) == pre_mask:
				break


		try:

			print()
			print("Basic stats:")
			print()
			print(f"  sd:  {round(self.sd,4)}")
			print(f"  med: {round(self.med,4)}")
			print()
			print(f"  zmed = (p - {round(self.med,4)}) / {round(self.sd,4)}")
			print()

		except TypeError:
			sys.exit(f"Error: failed to calculate basic statistics. This is likely due to very low measured alignment rates.")


		self.call_peaks()


	def calc_proportions(self):

		rg_size_c = get_global_depth(self.alignment_file, aggregate_by=['rg','length'])
		rg_c = get_global_depth(self.alignment_file, aggregate_by=['rg'])

		# props = list()  ## a list of all proportions in order
		# prop_d = dict() ## a dictionary of proportions by size
		if self.libraries == 'all':
			self.libraries = list(set(rg_c.keys()))

		for rg in self.libraries:
			self.master[rg] = list()

		self.master['prop'] = list()


		for size in self.master['sizes']:
			self.master['prop'].append(0)
			for rg in self.libraries:
				count = rg_size_c[(rg, str(size))]
	
				try:
					prop  = count / rg_c[rg]
				except ZeroDivisionError:
					prop  = 0
				# print(rg, size, count, prop, sep='\t')

				self.master[rg].append(prop)
				self.master['prop'][-1] += prop / len(rg_c)


	def calc_statistics(self):


		props = self.master['prop']
		props = [p if not self.master['mask'][i] else 0 for i,p in enumerate(props)]

		non_zero_props = [p for p in props if p > 0]

		if len(non_zero_props) < 2:
			self.med = None
			self.sd  = None

			self.master['zprop']     = [0 for p in props]
			self.master['peak']      = [None for p in props]
			self.master['candidate'] = [False for p in props]
			self.master['extension'] = [False for p in props]
			return


		sd  = stdev(non_zero_props)
		med = median(non_zero_props)

		self.sd  = sd
		self.med = med



		zprops = [(p - med) / sd for p in props]

		self.master['zprop'] = list()
		self.master['candidate'] = list()
		self.master['extension'] = list()
		self.master['peak'] = list()

		for i,z in enumerate(zprops):
			p = props[i]
			self.master['zprop'].append(z)
			self.master['candidate'].append(z > self.candidate_threshold and p > 0.01)
			self.master['extension'].append(z > self.extension_threshold)
			self.master['peak'].append(None)



	def mask_slope(self):

		# print('masking slopes...')
		for i, zprop in enumerate(self.master['zprop']):

			prop = self.master['prop'][i]

			try:
				next_prop = self.master['prop'][i+1]
			except IndexError:
				next_prop = prop

			try:
				pchange = next_prop / prop
			except ZeroDivisionError:
				pchange = None


			# print(i, round(zprop,3), round(prop,3), round(next_prop,3), round(pchange,3), sep='\t')


			if zprop < 0 and pchange and pchange > 1:
				break

			elif zprop < self.candidate_threshold and pchange and pchange > 1.1:
				break

			else:
				self.master['mask'][i] = True


		## this also masks any leftward or rightward peaks which slope out of the window
		## these cannot be resolved in the window, and are therefore ignored.

		rs = [enumerate(self.master['extension']), reversed(list(enumerate(self.master['extension'])))]

		for rang in rs:
			masked_positions = []
			nucleated = False
			for r, e in rang:

				if self.master['candidate'][r]:
					nucleated = True
				# print(r, e, self.master['candidate'][r], nucleated)

				if not e:
					if not nucleated:
						masked_positions = []

					break

				masked_positions.append(r)

			for p in masked_positions:
				self.master['mask'][p] = True



		# for r in range(self.min_size)
		# 	for i,e in enumerate(extensions):
			




	def call_peaks(self):
		props  = self.master['prop']
		zprops = self.master['zprop']

		masked = [i for i,m in enumerate(self.master['mask']) if m]

		candidates = [i for i,z in sorted(enumerate(zprops), key=lambda x:x[1], reverse=True) if z > self.candidate_threshold]
		extensions = [i for i,z in sorted(enumerate(zprops), key=lambda x:x[1], reverse=True) if z > self.extension_threshold]

		candidates = [i for i in candidates if i not in masked]
		extensions = [i for i in extensions if i not in masked]

		candidates = [i for i in candidates if props[i] > 0.01]

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
						## Breaks if the position does not meet the minimum proportion for a peak based on median-k.
						if r not in extensions:
							print(f'({r}) not a candidate')
							break
							

						p_curr = props[r]
						p_change = round((p_curr - p_last) / p_last  * 100,1)
						p_max    = round(p_curr / p_cand * 100, 1)


						p_last = p_curr

						## This filter helps to separate close peaks. I have disabled it because I believe it tends to split peaks which are likely connected.

						# if p_max < max_threshold and p_change > change_threshold:
						# 	print(f"({r}) below {max_threshold}% of max peak and above {change_threshold} change_threshold")
						# 	break

						## This makes a peak cutoff if the peak increases (saying this is likely a different peak)
						# if p_change > 0:
						# 	print(f'({r}) peak growing')
						# 	break

						## Breaks if the peak is extending into an already established peak
						if self.master['peak'][r]:
							print(f"({r}) pos in peak")
							break

						# print(direction, r, round(p_curr,4), round(p_last,4), p_change, p_max, sep='\t')
						self.master['peak'][r] = peak_i


						# ## also masking any contiguous peaks that do not resolve within the leftward window
						# if r == 0:
						# 	for i,p in enumerate(self.master['peak']):
						# 		if p == peak_i:
						# 			self.master['peak'][i] = None
						# 			self.master['mask'][i] = True

						# 	peak_i -= 1


	def peak_table(self, out_file=False):

		props  = self.master['prop']
		zprops = self.master['zprop']
		peaks  = self.master['peak']
		masks  = self.master['mask']

		candidates = [i for i,c in enumerate(self.master['candidate']) if c] 
		extensions = [i for i,e in enumerate(self.master['extension']) if e] 


		if out_file:
			outf = open(out_file, 'w')
			print("project\ti\tsize\tprop\tzero\tzmed\tcand\thyst\tpeak\tmask",file=outf)


		print()
		print("Sizes in terms of peaks:")
		print()
		print("i\tsize\tprop\tzero\tzmed\tcand\thyst\tpeak\tmask")
		print("=====================================================================")
		for i,s in enumerate(self.sizes):

			print(i, s, round(props[i],4), props[i] ==  0, round(zprops[i],4), i in candidates, i in extensions, peaks[i], masks[i], sep='\t')

			if out_file:
				print(self.project, i, s, round(props[i],4), props[i] ==  0, round(zprops[i],4), i in candidates, i in extensions, peaks[i], masks[i], sep='\t', file=outf)


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
		extensions = [i for i,e in enumerate(self.master['extension']) if e] 


		if any(peaks):
			max_peak = max([p for p in peaks if p])
		else:
			max_peak = 0

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


				val += 0.01

				try:
					z = (val - self.med) / self.sd
				except TypeError:
					z = 0

				if z > self.extension_threshold:
					pch = "•"
				if z > self.candidate_threshold:
					pch = "*"


				bar_string += pch

			peak_i = self.master['peak'][i]
			if peak_i:
				peak_name = f"peak{peak_i}"
			else:
				peak_name = ''

			if self.master['mask'][i]:
				peak_name = 'mask'

			print(s, round(p,3), round(zprops[i],3), peak_name, bar_string, sep='\t')
			print(s, round(p,3), round(zprops[i],3), peak_name, bar_string, sep='\t', file=outf)



		if out_file:
			outf.close()




@cli.command(group='Utilities', help_priority=5)



@optgroup.group('\n  Required',
				help='')


@optgroup.option("-o", "--output_directory", 
	# default=f"Annotation_{round(time())}", 
	required=False,
	type=click.UNPROCESSED, callback=validate_outdir,
	help="Directory name for annotation output. Defaults to the current directory, with this directory name as the project name.")

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

@optgroup.option('--all', is_flag=True, default=False, help='Override any annotation conditions and use all libraries for size-profiling.')

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

	if not alignment_file:
		sys.exit("Error: [alignment_file] not specified (and therefore alignment.depth.txt not found)")

	depth_file = Path(alignment_file).with_suffix(".depth.txt")

	if not depth_file.is_file():
		sys.exit(f"Error: depth file {str(depth_file)} not found, cannot calculate peaks...")


	# chromosomes, bam_rgs = get_chromosomes(alignment_file)

	libraries = []

	if len(annotation_conditions) == 0 or params['all']:
		libraries = 'all'
	else:
		for a in annotation_conditions:
			try:
				libraries += conditions[a]
			except KeyError:
				print(f"KeyError: condition '{a}' not found in conditions:")
				print(list(conditions.keys()))
				sys.exit()

	if len(libraries) == 0:
		sys.exit("Error: no libraries included in size profile. Try yasma.py readgroups to confirm you have the correct conditions and can detect library readgroups")

	pc = peakClass(ic.inputs['project_name'], alignment_file, libraries=libraries)

	pc.peak_table(alignment_file.with_suffix(".peak_table.txt"))
	pc.summarize_peaks(alignment_file.with_suffix(".peak_summary.txt"))
	pc.plot_proportions(alignment_file.with_suffix(".peak_plot.txt"))












