#!/usr/bin/env python3 

from modules.hairpin import *
from modules.annotate import *
from modules.precheck import *
from modules.context import *
from modules.count import *

version = "v0.6x"

@click.group()
def cli():
	pass


cli.add_command(precheck)
cli.add_command(annotate)
cli.add_command(hairpin)
cli.add_command(context)
cli.add_command(count)

if __name__ == '__main__':
	print()
	print(f"\033[1m-- annotator {version} --\033[0m")
	print()

	cli()








