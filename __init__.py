import os

__version__ = '0.0.1'

class Package:
	"""docstring for Package"""
	def __init__(self, path_to_package):
		self.parse_package(path_to_package)

	def parse_package(self, path_to_package):
		for f in os.listdir(path_to_package):
			if f.endswith(('ref.fa', 'ref.fasta', 'ref.fna')):
				self.ref_genome = os.path.join(path_to_package, f)
			if f.endswith('positions.txt'):
				self.var_list = os.path.join(path_to_package, f)
			if f.endswith('.log'):
				self.fasttree_log = os.path.join(path_to_package, f)
			if f.endswith(('.alleles.fa', '.alleles.fasta', '.alleles.fna')):
				self.variant_fasta = os.path.join(path_to_package, f)
			if f.endswith(('.tree', '.nwk')):
				self.fasttree_tree = os.path.join(path_to_package, f)
		if not hasattr(self, 'fasttree_tree'):
			print self.variant_fasta
			name = os.path.splitext(self.variant_fasta)[0]
			print 'fasttree -log {0}.log -gtr -nt {1} > {0}.tree'.format(name, self.variant_fasta)
			os.system('fasttree -log {0}.log -gtr -nt {1} > {0}.tree'.format(name, self.variant_fasta))
			self.fasttree_log = '%s.log' % (name)
			self.fasttree_tree = '%s.tree' % (name)

