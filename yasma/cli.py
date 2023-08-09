

import click
from click_option_group import optgroup

class customClickClass(click.Group):
	'''A custom class designed to add 4 features, based on SE posts from the developer of click. 
	1) adds grouping categorization to the help message: https://stackoverflow.com/questions/58745652/
	2) adds custom ordering to the help message: https://stackoverflow.com/questions/47972638/
	3) makes -h an acceptible flag for help
	4) makes -h the default flag for commands: https://stackoverflow.com/questions/50442401/'''

	def __init__(self, *args, **kwargs):

		context_settings = kwargs.setdefault('context_settings', {})
		if 'help_option_names' not in context_settings:
			context_settings['help_option_names'] = ['-h', '--help']
		self.help_flag = context_settings['help_option_names'][0]

		self.help_priorities = {}
		super(customClickClass, self).__init__(*args, **kwargs)

	def get_help(self, ctx):
		self.list_commands = self.list_commands_for_help
		return super(customClickClass, self).get_help(ctx)

	def list_commands_for_help(self, ctx):
		"""reorder the list of commands when listing the help"""
		commands = super(customClickClass, self).list_commands(ctx)
		return (c[1] for c in sorted(
			(self.help_priorities.get(command, 1), command)
			for command in commands))



	def command(self, *args, **kwargs):
		"""Gather the command help groups"""
		help_group = kwargs.pop('group', None)
		help_priority = kwargs.pop('help_priority', 1)
		help_priorities = self.help_priorities
		decorator = super(customClickClass, self).command(*args, **kwargs)

		def wrapper(f):
			cmd = decorator(f)
			cmd.help_group = help_group
			help_priorities[cmd.name] = help_priority
			return cmd

		return wrapper

	def format_commands(self, ctx, formatter):
		# Modified fom the base class method

		commands = []
		for subcommand in self.list_commands(ctx):
			cmd = self.get_command(ctx, subcommand)
			if not (cmd is None or cmd.hidden):
				commands.append((subcommand, cmd))

		if commands:
			longest = max(len(cmd[0]) for cmd in commands)
			# allow for 3 times the default spacing
			limit = formatter.width - 6 - longest

			groups = {}
			for subcommand, cmd in commands:
				help_str = cmd.get_short_help_str(limit)
				subcommand += ' ' * (longest - len(subcommand))
				groups.setdefault(
					cmd.help_group, []).append((subcommand, help_str))

			with formatter.section('Commands'):
				for group_name, rows in groups.items():
					with formatter.section(group_name):
						formatter.write_dl(rows)

	def parse_args(self, ctx, args):
		if len(args) <= 1:
			args.append(self.help_flag)

		glob_options = set(["-tl", '--trimmed_libraries', "-ul", '--untrimmed_libraries'])

		# print(args)
		last_is_glob = False
		new_args = []
		for arg in args:

			if arg.startswith("-"):
				new_args.append(arg)
				if arg in glob_options:
					last_is_glob = True
					new_args.append("")
				else:
					last_is_glob = False

			else:

				if last_is_glob:
					new_args[-1] += " " + arg
				else:
					new_args.append(arg)

		# print(new_args)

		# args = ['test-function', '-l', 'SRR3222443_trimmed.fastq SRR3222444_trimmed.fastq']


		return super(customClickClass, self).parse_args(ctx, new_args)

	


@click.group(cls=customClickClass)
def cli():
	"""YASMA (Yet Another small RNA Annotator)"""





