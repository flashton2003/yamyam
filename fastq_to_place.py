import os
import datetime

def align_fastqs(args, package):
	if args.fastq2 != None:
		os.system('bwa mem -t 4 {0} {1} {2} > {3}/{4}.sam'.format(package.ref_genome, args.fastq1, args.fastq2, args.output_dir, args.sample_name))
	elif args.fastq2 == None:
		os.system('bwa mem -t 4 -x ont2d {0} {1} > {2}/{3}.sam'.format(package.ref_genome, args.fastq1, args.output_dir, args.sample_name))


def parse_sam(args, package):
	this_dir = os.path.dirname(os.path.realpath(__file__))
	os.system('cat {0}/{1}.sam | python {2}/get_alleles_from_sam.py {1} {3} {0}/{1}.var.fa'.format(args.output_dir, args.sample_name, this_dir, package.var_list))

def run_taxit(args, package):
	now = str(datetime.datetime.now())
	os.system('taxit create --stats-type FastTree -s {0} -f {1} -l genome -t {2} -P {3}/taxit'.format(package.fasttree_log, package.variant_fasta, package.fasttree_tree, args.output_dir))

def run_pplacer(args, package):
	add_one_in = os.path.splitext(package.variant_fasta.split('/')[-1])[0] + 'add_one_in'
	os.system('cat {0}/{1}.var.fa {2} > {0}/taxit/{3}.fa'.format(args.output_dir, args.sample_name, package.variant_fasta, add_one_in))
	os.system('pplacer -c {0}/taxit --out-dir {0}/taxit {0}/taxit/{1}.fa'.format(args.output_dir, add_one_in))

def run_guppy(args, package):
	add_one_in = os.path.splitext(package.variant_fasta.split('/')[-1])[0] + 'add_one_in'
	os.system('guppy tog --out-dir {0} {0}/taxit/{1}.jplace'.format(args.output_dir, add_one_in))
	
