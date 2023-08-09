# Genomic context module

import os
import click
from subprocess import PIPE, Popen, call
from pathlib import Path
from os.path import isfile, isdir
from pprint import pprint

from .generics import *
from .cli import cli

@cli.command(group='Calculation', help_priority=3)

@click.option("-g", "--gene_annotation_file", 
	required=False, 
	type=click.UNPROCESSED, callback=validate_path,
	help='Gene annotation file in gff3 format. Tested with NCBI annotation formats.')

@click.option("-o", "--output_directory", 
	required=True, 
	type=click.Path(),
	help="Directory name for annotation output")


@click.option("--intergenic_distance", 
	default=1000,
	help="Distance from a gene for a locus to be defined as 'intergenic'. Default 1000 bp.")

# @click.option("-m", "--method", 
# 	default="Poisson", 
# 	help="Annotator algorithm used (Poisson or Dicer)")

@click.option("-f", "--force",
	is_flag=True,
	help='Force remake of supporting files')


def context(**params):
	"""Compares annotations to identify cluster genomic context."""


	rc = requirementClass()
	rc.add_bedtools()
	rc.check()

	ic = inputClass(params)
	ic.check(['gene_annotation_file'])


	output_directory        = str(ic.output_directory)
	gene_annotation_file    = ic.inputs["gene_annotation_file"]

	force                   = params['force']
	intergenic_distance     = params['intergenic_distance']




	ann_file = f"{output_directory}/peak/loci.gff3"

	if not isfile(ann_file):
		sys.exit(f"Annotation file {ann_file} does not exist. Are you sure the output folder contains an annotation?")

	results_file = f"{output_directory}/peak/loci.txt"

	if not isfile(results_file):
		sys.exit(f"results file ({results_file}) not found")

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
		cds_file  = gene_annotation_file.replace(".gff3", ".cds.gff3")


		if not isfile(mRNA_file) or not isfile(exon_file) or not isfile(cds_file):
			sort()


		def write_file(file, key):
			lines_written = 0
			if not isfile(file):
				with open(file, 'w') as outf:
					with open(gene_annotation_file, 'r') as f:
						for line in f:
							line = line.strip()

							if line.startswith("#"):
								print(line, file=outf)

							elif line.split('\t')[2] == key:
								print(line, file=outf)

								lines_written += 1
								# print(lines_written)

				if lines_written < 100:
					print(f"WARNING: less than 100 lines written for {file} in the {key} sub-annotation. Is this a fully annotated NCBI-derived GFF3?")

		write_file(mRNA_file, 'mRNA')
		write_file(exon_file, 'exon')
		write_file(cds_file, 'CDS')
		
		return(mRNA_file, exon_file, cds_file)

	mRNA_file, exon_file, cds_file = make_subsets()


	def closest(file, ID):

		closest_d = {}

		call= ['bedtools', 'closest', '-a', ann_file, '-b', file, '-d']

		print(" ".join(call))
		p = Popen(call, stdout=PIPE, stderr=PIPE, encoding=ENCODING)
		out, err = p.communicate()


		for o in out.strip().split("\n"):

			o = o.strip().split("\t")

			s_strand = o[6]
			m_strand = o[15]

			sa = parse_attributes(o[8])
			ma = parse_attributes(o[17])

			# mRNA_id\ttranscript\tdistance

			d = {}

			d[f'{ID}_id'] = ma['ID']
			d['transcript'] = ma['orig_transcript_id']
			d[f'{ID}_distance'] = int(o[-1])


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

			# pprint(d)
			# sys.exit()

		return(closest_d)

	def intersect(file, ID):

		inter_d = {}

		call = ['bedtools', 'intersect', '-a', ann_file, '-b', file, '-wao', '-f', '0.1']


		p = Popen(call, stdout=PIPE, stderr=PIPE, encoding=ENCODING)
		out, err = p.communicate()
		# print(out)
		# sys.exit()

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


				d[f'{ID}_id'] = ma['ID']
				d[f'{ID}_overlap'] = int(o[-1])

				inter_d[sa['ID']] = d



			# sys.exit()
		return(inter_d)


	print("Finding intersections for...")
	print("   mRNAs")
	mRNA_d = closest(mRNA_file, ID='mRNA')
	print("   exons")
	exon_d = intersect(exon_file, ID='exon')
	print("   CDSs")
	cds_d = intersect(cds_file, ID='cds')


	Path(output_directory+ "/context/").mkdir(parents=True, exist_ok=True)
	output_file = f"{output_directory}/context/context.txt"

	print()
	print(f'Printing to...\n  {output_file}')

	with open(output_file, 'w') as outf:
		print("\t".join(['cluster','transcript', 'mRNA_id', 'mRNA_distance', 'exon_id', 'exon_overlap', 'cds_id', 'cds_overlap', 's_strand','m_strand','match','category']), file=outf)

		with open(results_file, 'r') as f:
			header = f.readline()

			for line in f:
				line = line.strip().split("\t")
				cluster = line[1]
				# print(cluster)

				try:
					m_d = mRNA_d[cluster]
				except KeyError:
					m_d = {}

				try:
					i_d = exon_d[cluster]
				except KeyError:
					i_d = {}

				try:
					c_d = cds_d[cluster]
				except KeyError:
					c_d = {}



				m_d.update(i_d)
				m_d.update(c_d)
				d = m_d


				out_line = [cluster]
				for k in ['transcript', 'mRNA_id', 'mRNA_distance', 'exon_id', 'exon_overlap', 'cds_id', 'cds_overlap', 's_strand','m_strand','match']:

					try:
						out_line.append(d[k])
					except KeyError:
						d[k]= ''
						out_line.append("")



				## finding gene relationship

				if d['mRNA_distance'] > intergenic_distance:
					gene_relationship = 'intergenic'

				elif 0 < d['mRNA_distance'] <= intergenic_distance:
					gene_relationship = 'near-genic'

				elif d['mRNA_distance'] == 0:
					gene_relationship = 'genic'

				else:
					pprint(d)
					sys.exit('gene relationship error')



				## finding strand relationship

				if gene_relationship != 'genic':
					stranding = "NA"

				elif d['match'] == '=':
					stranding = 'sense'

				elif d['match'] == '!':
					stranding = 'antisense'

				elif d['match'] == '?':
					stranding = 'unstranded'

				else:
					pprint(d)
					sys.exit('stranding error')


				## finding intragene context

				if gene_relationship != 'genic':
					intragene_context = "NA"

				elif d['exon_overlap'] == '':
					intragene_context = 'intronic'

				elif d['exon_overlap'] > 0 and d['cds_overlap'] == '':
					intragene_context = 'exon-UTR'

				elif d['exon_overlap'] > 0 and d['cds_overlap'] > 0:
					intragene_context = 'exon-CDS'

				else:
					pprint(d)
					sys.exit('intragene_context error')



				if gene_relationship != 'genic':
					category = gene_relationship

				else:
					category = f"{stranding}_{intragene_context}"




				out_line.append(category)



				# print(out_line)
				out_line = '\t'.join(map(str, out_line))
				# sys.exit()

				# print(out_line, sep='\t')
				print(out_line, sep='\t', file=outf)







