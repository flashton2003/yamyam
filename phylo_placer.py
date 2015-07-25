import argparse
import os
from fastq_to_place import align_fastqs, parse_sam, run_taxit, run_pplacer, run_guppy
from __init__ import __version__, Package

def run_command(args):
    package = Package(args.place_me_package)
    print args.fastq1
    args.sample_name = os.path.splitext(args.fastq1)[0]
    print args.sample_name
    sys.exit()
    if args.command == 'place_me':
        align_fastqs(args, package)
        parse_sam(args, package)
        run_taxit(args, package)
        run_pplacer(args, package)
        run_guppy(args, package)

def main():
    parser = argparse.ArgumentParser(prog='phylo-placer')
    parser.add_argument("-v", "--version", help="Installed snapperdb version", action="version",
                        version="%(prog)s " + str(__version__))
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command')

    # # metavar is a variable name for use in help
    # # nargs = '+' means that fastqs will take mutliple files and return a list and return an error if not at least 1

    parser_place_me = subparsers.add_parser('place_me', help='Takes fastqs and a place me package, provides placement on tree.')
    parser_place_me.add_argument('-1', dest='fastq1', help='first sequence file (if you give 1 seq file, assumes ONT)', required = True)
    parser_place_me.add_argument('-2', dest='fastq2', help='second sequence file (if you give 2 seq files, assumes illumina)', default = None)
    parser_place_me.add_argument('-p', dest='place_me_package', required=True,
                                    help='The path to a place me package')
    parser_place_me.add_argument('-o', dest='output_dir', help='output directory - default to cwd', default = os.getcwd())
    parser_place_me.add_argument('-s', dest='sample_name', help='sample_name - default to fastq name')
    args = parser.parse_args()
    run_command(args)

if __name__ == '__main__':
    main()
