#!/usr/bin/env python3 

from modules.hairpin import *
from modules.annotate import *
from modules.precheck import *


@click.group()
def cli():
	pass


cli.add_command(precheck)
cli.add_command(annotate)
cli.add_command(hairpin)

if __name__ == '__main__':
	cli()








