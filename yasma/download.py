## SRA tools utility for downloading SRRs

from .generics import *

import gzip
from shutil import rmtree
from time import sleep


@cli.command(group='Processing', help_priority=1)


@optgroup.group('\n  Basic options',
				help='')

@optgroup.option("-s", "--srrs", 
	required=False, 
	multiple=True,
	help='NCBI SRA codes for libraries. These will almost certainly start with SRR or ERR.')

@optgroup.option("-o", "--output_directory", 
	# default=f"Annotation_{round(time())}", 
	required=False,
	type=click.UNPROCESSED, callback=validate_outdir,
	help="Directory name for annotation output. Defaults to the current directory, with this directory name as the project name.")


@optgroup.option('--include_quals', is_flag=True, default=False, help='Download libraries as .fastq format (default is only .fasta)')

@optgroup.option('--zipped/--unzipped', is_flag=True, default=False, help='Whether to compress downloaded files (default is uncompressed)')



def download(**params):
	'''Download libraries from the NCBI SRA using their SRR code'''

	rc = requirementClass()
	rc.add_sratools()
	rc.check()

	sratools_version = int([v[2][0] for v in rc.reqs if v[0] == 'fasterq-dump'][0])

	ic = inputClass(params)
	ic.check(['srrs'])


	output_directory  = ic.output_directory
	srrs              = list(ic.inputs['srrs'])



	untrimmed_dir = Path(output_directory, "untrimmed")
	untrimmed_dir.mkdir(parents=True, exist_ok=True)

	download_dir = Path(output_directory, "download")
	download_dir.mkdir(parents=True, exist_ok=True)

	log_file = Path(output_directory,"download/log.txt")
	sys.stdout = Logger(log_file)

	for srr in srrs:
		lock_file_1 = Path(download_dir, srr, f"{srr}.sra.lock")
		lock_file_2 = Path(download_dir, srr, f"{srr}.sralite.lock")

		if lock_file_1.is_file() or lock_file_2.is_file():
			rmtree(str(Path(download_dir, srr)))




	untrimmed_libraries = []

	for i, srr in enumerate(srrs):

		suffix = 'fastq'
		if not params['include_quals']:

			if sratools_version < 3:
				print('warning: fasterq-dump version is older than 3.x.x, and will only output as fastq')
			else:
				suffix = 'fasta'

		unzipped_file = Path(untrimmed_dir, f"{srr}.{suffix}")
		zipped_file   = Path(untrimmed_dir, f"{srr}.fq.gz")

		print(f"\n  downloading {i+1} of {len(srrs)}  ")

		if params['zipped'] and zipped_file.is_file():
			print(' ', zipped_file, 'found...')
			untrimmed_libraries.append(zipped_file)
			continue

		elif not params['zipped'] and unzipped_file.is_file():
			print(' ', unzipped_file, 'found...')
			untrimmed_libraries.append(unzipped_file)
			continue




		valid = False
		try_counter = 0

		while not valid:
			try_counter += 1
			call = ['prefetch', "-O", str(download_dir), srr]

			print(f"calling (attempt {try_counter}): ", " ".join(call))

			p = Popen(call, encoding=ENCODING, stdout=PIPE)
			for line in p.stdout:
				print("  ", line.strip())

				if f"'{srr}.lite is found locally":
					valid=True
				if f"'{srr}' is found locally" in line:
					valid=True
				if f"'{srr}' was downloaded successfully" in line:
					valid=True

			p.wait()

			if not valid:
				n = try_counter * 10
				print(f'prefetch failed. Waiting {n} seconds')
				sleep(n)

			if try_counter > 25:
				sys.exit(f"ERROR: could not prefetch {srr}")

		



		call = ['fasterq-dump'] + [str(Path(download_dir, srr)), '-O', str(untrimmed_dir)]
		if suffix == 'fasta':
			call += ['--fasta']


		print()
		print()
		print("calling: ", " ".join(call))

		p = Popen(call, encoding=ENCODING, stdout=PIPE)
		for line in p.stdout:
			print("  ", line.strip())
		p.wait()


		print()
		print()


		if not params['zipped']:
			untrimmed_libraries.append(unzipped_file)

		else:
			print("zipping...")
			try:
				Path(untrimmed_dir, f"{srr}_1.{suffix}").rename(Path(untrimmed_dir, f"{srr}.{suffix}"))
			except:
				pass

			try:
				Path(untrimmed_dir, f"{srr}_2.{suffix}").unlink()
			except:
				pass

			try:
				Path(untrimmed_dir, f"{srr}_3.{suffix}").unlink()
			except:
				pass

			untrimmed_libraries.append(zipped_file)

			print(f"  {unzipped_file} ->")
			print(f"        {zipped_file}")

			with open(unzipped_file, 'rb') as unzippedf:
				with gzip.open(zipped_file, 'wb') as zippedf:
					zippedf.writelines(unzippedf)

			unzipped_file.unlink()



	print(f"writing untrimmed_libraries to inputs.json")

	ic.inputs['untrimmed_libraries'] = untrimmed_libraries
	ic.write()



















