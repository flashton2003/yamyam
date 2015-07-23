import argparse
import os
from fastq_to_place import align_fastqs, parse_sam
from __init__ import __version__, Package

def run_command(args):
    package = Package(args.place_me_package)
    args.sample_name = args.fastq1.split('.')[0]
    if args.command == 'place_me':
        # align_fastqs(args, package)
        parse_sam(args, package)

def main():
    parser = argparse.ArgumentParser(prog='phylo-placer')
    parser.add_argument("-v", "--version", help="Installed snapperdb version", action="version",
                        version="%(prog)s " + str(__version__))
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command')

    # # metavar is a variable name for use in help
    # # nargs = '+' means that fastqs will take mutliple files and return a list and return an error if not at least 1

    parser_place_me = subparsers.add_parser('place_me', help='Takes fastqs and a place me package, provides placement on tree.')
    parser_place_me.add_argument('-1', dest='fastq1', help='first fastq file', required = True)
    parser_place_me.add_argument('-2', dest='fastq2', help='second fastq file', default = None)
    parser_place_me.add_argument('-p', dest='place_me_package', required=True,
                                    help='The path to a place me package')
    parser_place_me.add_argument('-o', dest='output_dir', help='output directory - default to cwd', default = os.getcwd())
    parser_place_me.add_argument('-s', dest='sample_name', help='sample_name - default to fastq name')
    args = parser.parse_args()
    run_command(args)

if __name__ == '__main__':
    main()
