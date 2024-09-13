## SRA tools utility for downloading SRRs

import sys

import click
from click_option_group import optgroup

from pprint import pprint

from .generics import *
from .cli import cli

import gzip
from shutil import rmtree


@cli.command(group='Processing', help_priority=2)


@optgroup.group('\n  Basic options',
				help='')

@optgroup.option("-s", "--srrs", 
	required=False, 
	multiple=True,
	help='NCBI SRA codes for libraries. These will almost certainly start with SRR or ERR.')

@optgroup.option("-o", "--output_directory", 
	# default=f"Annotation_{round(time())}", 
	required=True,
	type=click.Path(),
	help="Directory name for annotation output.")


@optgroup.option('--unzipped', is_flag=True, default=False, help='Do not compress downloaded files (default is to compress)')



def download(**params):
	'''Tool to check untrimmed-libraries for 3' adapter content.'''

	pprint(params)
	rc = requirementClass()
	rc.add_sratools()
	rc.check()

	ic = inputClass(params)
	ic.check(['srrs'])


	output_directory  = ic.output_directory
	srrs              = list(ic.inputs['srrs'])


	untrimmed_dir = Path(output_directory, "untrimmed")
	untrimmed_dir.mkdir(parents=True, exist_ok=True)

	download_dir = Path(output_directory, "download")

	for srr in srrs:
		prefetch_file = Path(download_dir, srr, f"{srr}.sra.lock")

		if prefetch_file.is_file():
			rmtree(str(Path(download_dir, srr)))



	call = ['prefetch', "-O", str(download_dir)] + srrs

	print("calling: ", " ".join(call))

	p = Popen(call, encoding=ENCODING, stdout=PIPE)
	for line in p.stdout:
		print("  ", line.strip())
	p.wait()



	call = ['fasterq-dump'] + [str(Path(download_dir, s)) for s in srrs] + ['-O', str(untrimmed_dir)]

	print()
	print()
	print("calling: ", " ".join(call))

	p = Popen(call, encoding=ENCODING, stdout=PIPE)
	for line in p.stdout:
		print("  ", line.strip())
	p.wait()


	print()
	print()
	print("zipping...")

	untrimmed_libraries = []

	if not params['unzipped']:
		for i,srr in enumerate(srrs):

			try:
				Path(untrimmed_dir, f"{srr}_1.fastq").rename(Path(untrimmed_dir, f"{srr}.fastq"))
			except:
				pass

			try:
				Path(untrimmed_dir, f"{srr}_2.fastq").unlink()
			except:
				pass

			try:
				Path(untrimmed_dir, f"{srr}_3.fastq").unlink()
			except:
				pass

			unzipped_file = Path(untrimmed_dir, f"{srr}.fastq")
			zipped_file   = Path(untrimmed_dir, f"{srr}.fq.gz")
			untrimmed_libraries.append(zipped_file)

			print(f"  {i+1} of {len(srrs)}")
			print(f"  {unzipped_file} ->")
			print(f"        {zipped_file}")
			with open(unzipped_file, 'rb') as unzippedf:
				with gzip.open(zipped_file, 'wb') as zippedf:
					zippedf.writelines(unzippedf)

			unzipped_file.unlink()

	else:
		for i,srr in enumerate(srrs):

			unzipped_file = Path(untrimmed_dir, f"{srr}.fastq")
			untrimmed_libraries.append(zipped_file)


	print(f"writing untrimmed_libraries to inputs.json")

	ic.inputs['untrimmed_libraries'] = untrimmed_libraries
	ic.write()



















