#!/usr/bin/env python3
import sys
import os
import logging
import argparse
import importlib
from unittest.mock import patch

""" main()
cfDNAFE main commandline program.
"""

VERSION = '0.1.0'

commands = [
    "bam",
    "fsc",
    "fsd",
    "fsr",
    "cnv",
    "wps",
    "ocf",
    "mutation",
    "meth"
    ]

def main():
    if len(sys.argv) < 2 or (len(sys.argv) == 2 and sys.argv[1] in ('-h', '--help')):
        print_help()
        return
    elif len(sys.argv) == 2 and sys.argv[1] == '--version':
        print('cfDNAFE version', VERSION)
        return

    parser = argparse.ArgumentParser(
        description='cfDNAFE tool',
        usage='cfDNAFE <command> [<args>]')
    parser.add_argument('command', help='Subcommand to run')
    args = parser.parse_args(sys.argv[1:2])
    try:
        with patch.object(sys, 'argv', sys.argv[1:]):
            importlib.import_module(args.command).main()
            
    except ModuleNotFoundError as e:
        print(args.command)
        print(str(e))
        if args.command not in str(e):
            raise e
        print_invalid_command(args.command)
        print_help()

    except ValueError as e:
        eprint(f'Invalid input argument\n{e}')
        return 1


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def print_invalid_command(command):
    eprint('Invalid command:', f'\033[01;31m{command}\033[00m')
    from difflib import get_close_matches
    closets = [x for x in get_close_matches(command, commands)]
    if closets:
        eprint(f'did you mean \033[01;32m{closets[0]}\033[00m?')

def print_help(command=None):
    print(f'A tool for extracting cfDNA feature, version {VERSION}')
    msg = '\nUsage: cfDNAFE <command> [<args>]'
    msg += '\nrun cfDNAFE <command> -h for more information'
    msg += '\nOptional commands:'
    msg += '\nbam, extracting end motif, breakpoint motif and Motif-Diversity Score'
    msg += '\nfsc, extracting fragment size coverage'
    msg += '\nfsd, extracting fragment size distribution'
    msg += '\nfsr, extracting fragment size ratio'
    msg += '\ncnv, extracting copy number variations'
    msg += '\nwps, extracting window protect score'
    msg += '\nocf, extracting orientation-aware cfDNA fragmentation'
    msg += '\nmutation, extracting 96 single base substitution mutation (SBS) profile and a mutation signature profile'
    msg += '\nmeth, extracting UXM fragment-level'
    print(msg)

    return 1


if __name__ == '__main__':
    main()









