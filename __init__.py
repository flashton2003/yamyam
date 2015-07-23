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
			if f.endswith('var_list'):
				self.var_list = os.path.join(path_to_package, f)
			if f.endswith('.log'):
				self.fasttree_log = os.path.join(path_to_package, f)
			if f.endswith('.var.fa'):
				self.variant_fasta = os.path.join(path_to_package, f)
			if f.endswith(('.tree', '.nwk')):
				self.fasttree_tree = os.path.join(path_to_package, f)

