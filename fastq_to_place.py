import os


def align_fastqs(args, package):
	if args.fastq2 != None:
		os.system('bwa mem -t 4 -x ont2d {0} {1} {2}'.format(package.ref_genome, args.fastq1, args.fastq2, args.seq_tech))
	elif args.fastq2 == None:
		os.system('bwa mem -t 4 -x ont2d {0} {1}'.format(package.ref_genome, args.fastq1))


def parse_sam(args, package):
	this_dir = os.path.dirname(os.path.realpath(__file__))
	os.system('cat {0}/{3}.sam | python {1}/get_alleles_from_sam.py {3} {2} blah.vcf {0}/{3}.fa'.format(args.output_dir, this_dir, package.var_list, args.sample_name))
