#!/usr/bin/env python3


"""Convenience wrapper for running annotator directly from source tree."""


import os
from git import Repo


git_dir = os.path.dirname(os.path.realpath(__file__))
git_repo = Repo(git_dir)
git_commit = git_repo.head.commit.tree


print()
print(f"git_dir:    {git_dir}")
print(f"git_commit: {git_commit}")
print()


from yasma.cli import cli

if __name__ == '__main__':
	cli()
