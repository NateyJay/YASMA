# Genomic context module

import os
import click
from subprocess import PIPE, Popen, call
from pathlib import Path
from os.path import isfile, isdir
from pprint import pprint

from modules.generics import *


@click.command()

@click.option("-g", "--gene_annotation_file", 
	required=True, 
	type=click.Path(exists=True),
	help='Gene annotation file in gff3 format. Tested with NCBI annotation formats.')

@click.option("-o", "--output_directory", 
	default=f"Annotation_{round(time())}", 
	type=click.Path(),
	help="Directory name for annotation output")

@click.option("-f", "--force",
	is_flag=True,
	help='Force remake of supporting files')


def context(gene_annotation_file, output_directory, force):
	"""Compares annotations to identify cluster genomic context"""

	ann_file = f"{output_directory}/Annotation.gff3"

	if not isfile(ann_file):
		sys.exit(f"Annotation file {ann_file} does not exist. Must run 'annotate' module first and specify the same output folder.")



	# 1) make subfiles
	# 2) find closest mRNA
	# 3) find intersects


	def parse_attributes(x):
		if x == '.':
			return(None)
		d = {}
		x = x.strip().split(";")
		for xx in x:
			key, val = xx.split("=")
			d[key] = val
		return(d)

	def sort():
		temp_file = "temp_sort.gff3"
		call = ['bedtools', 'sort', '-i', gene_annotation_file]

		print("input gff not sorted - sorting now with:")
		print(" ".join(call))

		with open(temp_file, 'w') as outf:
			p = Popen(call, stdout=outf)
			p.wait()

		os.rename(temp_file, gene_annotation_file)


	def make_subsets():
		mRNA_file = gene_annotation_file.replace(".gff3", ".mRNA.gff3")
		exon_file = gene_annotation_file.replace(".gff3", ".exon.gff3")

		if not isfile(mRNA_file) or not isfile(exon_file):
			sort()

		if not isfile(mRNA_file):
			with open(mRNA_file, 'w') as outf:
				with open(gene_annotation_file, 'r') as f:
					for line in f:
						line = line.strip()

						if line.startswith("#"):
							print(line, file=outf)

						elif line.split('\t')[2] == 'mRNA':
							print(line, file=outf)

		
		if not isfile(exon_file):
			with open(exon_file, 'w') as outf:
				with open(gene_annotation_file, 'r') as f:
					for line in f:
						line = line.strip()

						if line.startswith("#"):
							print(line, file=outf)

						elif line.split('\t')[2] == 'exon':
							print(line, file=outf)

		return(mRNA_file, exon_file)

	mRNA_file, exon_file = make_subsets()



	# sys.exit()

	def closest():

		closest_d = {}

		call= ['bedtools', 'closest', '-a', ann_file, '-b', mRNA_file, '-d']

		p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')
		out, err = p.communicate()


		for o in out.strip().split("\n"):
			o = o.strip().split("\t")
			# print(o)

			s_strand = o[6]
			m_strand = o[15]

			sa = parse_attributes(o[8])
			ma = parse_attributes(o[17])

			# mRNA_id\ttranscript\tdistance

			d = {}

			d['mRNA_id'] = ma['ID']
			d['transcript'] = ma['orig_transcript_id']
			d['distance'] = int(o[-1])


			if s_strand == '.' or m_strand == '.':
				match = "?"
			elif s_strand == m_strand:
				match = "="
			else:
				match = '!'

			d['s_strand'] = s_strand
			d['m_strand'] = m_strand

			d['match'] = match

			closest_d[sa['ID']] = d

		return(closest_d)





	closest_d = closest()

	# for key, val in closest_d.items():
	# 	print(key, "\t".join(map(str, val)))



	def intersect():

		inter_d = {}

		call = ['bedtools', 'intersect', '-a', ann_file, '-b', exon_file, '-wao']


		p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')
		out, err = p.communicate()


		for o in out.strip().split("\n"):
			o = o.strip().split("\t")
			# print(o)
			# print(o[8])
			# print(o[17])
			overlap = int(o[-1])
			sa = parse_attributes(o[8])
			ma = parse_attributes(o[17])


			if ma:
				d = {}


				d['exon_id'] = ma['ID']
				d['overlap'] = int(o[-1])

				inter_d[sa['ID']] = d



			# sys.exit()
		return(inter_d)

	inter_d = intersect()

	with open(f"{output_directory}/GenomicContext.txt", 'w') as outf:
		print("\t".join(['cluster','transcript', 'mRNA_id', 'distance', 'exon_id', 'overlap', 's_strand','m_strand','match','category']), file=outf)

		with open(f"{output_directory}/Results.txt", 'r') as f:
			header = f.readline()

			for line in f:
				line = line.strip().split("\t")
				cluster = line[0]




				try:
					c_d = closest_d[cluster]
				except KeyError:
					c_d = {}


				try:
					i_d = inter_d[cluster]
				except KeyError:
					i_d = {}



				c_d.update(i_d)
				d = c_d

				# pprint(d)


				out_line = [cluster]
				for k in ['transcript', 'mRNA_id', 'distance', 'exon_id', 'overlap', 's_strand','m_strand','match']:

					try:
						out_line.append(d[k])
					except KeyError:
						d[k]= ''
						out_line.append("")


				if d['distance'] > 1000:
					category = 'intergenic'

				elif 0 < d['distance'] <= 1000:
					category = 'near-genic'

				elif d['overlap'] == '' and d['match'] == '=':
					category = 'intronic-sense'

				elif d['overlap'] == '' and d['match'] == '!':
					category = 'intronic-antisense'

				elif d['overlap'] == '' and d['match'] == '?':
					category = 'intronic-unstranded'

				elif d['overlap'] > 0 and d['match'] == '=':
					category = 'exonic-sense'

				elif d['overlap'] > 0 and d['match'] == '!':
					category = 'exonic-antisense'

				elif d['overlap'] > 0 and d['match'] == '?':
					category = 'exonic-unstranded'

				# elif d['match'] == '=':
				# 	category = 'intronic-sense'

				# elif d['match'] == '!':
				# 	category = 'intronic-antisense'

				# elif d['match'] == '?':
				# 	category = 'intronic-unstranded'

				else:
					print("WHAT ARE YOU?")
					pprint(d)
					sys.exit()

				out_line.append(category)



				# print(out_line)
				out_line = '\t'.join(map(str, out_line))
				# sys.exit()

				print(out_line, sep='\t')
				print(out_line, sep='\t', file=outf)







