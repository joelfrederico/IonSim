#!/usr/bin/env python3

import subprocess
import shlex
import argparse


def runcommand(command, **kwargs):
    print("Running: " + command)
    subprocess.run(shlex.split(command), **kwargs)


def update(part):

    try:
        runcommand("bumpversion {} --no-commit --no-tag".format(part), check=True)
    except:
        return
    
    runcommand("autoreconf -fiv")

    runcommand("git add Makefile.in config.h.in~ configure")
    
    runcommand("git checkout configure.ac .bumpversion.cfg")

    runcommand("git commit -m \"Automatically autoreconf'ed\"")

    runcommand("bumpversion {}".format(part))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Manage bumpversion/autotools interface")
    parser.add_argument('part', nargs=1, choices=['major', 'minor', 'patch'])

    args = parser.parse_args()
    part = args.part[0]

    print(part)
    update(part)
