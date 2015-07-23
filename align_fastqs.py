import os


def align_fastqs(args, package):
	if args.fastq2 != None:
		os.system('bwa mem -t 4 -x ont2d {0} {1} {2}'.format(package.ref_genome, args.fastq1, args.fastq2))
	elif args.fastq2 == None:
		os.system('bwa mem -t 4 -x ont2d {0} {1}'.format(package.ref_genome, args.fastq1))
