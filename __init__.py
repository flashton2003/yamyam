import os

__version__ = '0.0.1'

class Package:
	"""docstring for Package"""
	def __init__(self, path_to_package):
		self.parse_package(path_to_package)

	

	def parse_package(self, path_to_package):
		for f in os.listdir(path_to_package):
			if f.endswith(('.fa', '.fasta', '.fna')):
				self.ref_genome = os.path.join(path_to_package, f)

		pass
