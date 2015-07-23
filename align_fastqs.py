import os


def align_fastqs(args, package):
	if args.fastq2 != None:
		os.system('bwa mem -t 4 -x ont2d {0} {1} {2}'.format(package.ref_genome, args.fastq1, args.fastq2, args.seq_tech))
	elif args.fastq2 == None:
		os.system('bwa mem -t 4 -x ont2d {0} {1}'.format(package.ref_genome, args.fastq1))


def parse_sam(args):
	os.system('cat /home/ubuntu/typhi_illumina.sam | python ~/yamyam/get_alleles_from_sam.py typhi_nanopore yamyam-packages/salmonella-typhi/2015-07-22.ebg_13_snps.list.var_list blah output.txt illumina')
