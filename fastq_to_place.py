import os
import datetime

def align_fastqs(args, package):
	if args.fastq2 != None:
		os.system('bwa mem -t 4 -x ont2d {0} {1} {2} > {3}/{4}.sam'.format(package.ref_genome, args.fastq1, args.fastq2, args.output_dir, args.sample_name))
	elif args.fastq2 == None:
		os.system('bwa mem -t 4 -x ont2d {0} {1} > {3}/{4}.sam'.format(package.ref_genome, args.fastq1, args.output_dir, args.sample_name))


def parse_sam(args, package):
	this_dir = os.path.dirname(os.path.realpath(__file__))
	os.system('cat {0}/{1}.sam | python {2}/get_alleles_from_sam.py {3} {4} blah.vcf {5}/{6}.fa {7}'.format(args.output_dir, args.sample_name, this_dir, package.var_list, args.output_dir, args.sample_name, args.seq_tech))

def run_taxit(args, package):
	now = str(datetime.datetime.now())
	os.sytem('taxit create --stats-type FastTree -s {0} -f {1} -l genome -t {2} -P {3}/taxit'.format(package.fasttree_log, package.variant_fasta, package.fasttree_tree, args.output_dir))