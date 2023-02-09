#!/usr/bin/env python3 

from modules.hairpin import *
from modules.annotate import *
from modules.precheck import *
from modules.context import *
from modules.count import *
from modules.poisson import *

version = "v0.7x"



CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
	pass


cli.add_command(precheck)
cli.add_command(annotate)
cli.add_command(hairpin)
cli.add_command(context)
cli.add_command(count)
cli.add_command(poisson)

if __name__ == '__main__':
	print()
	print(f"\033[1m-- annotator {version} --\033[0m")
	print()

	cli()








