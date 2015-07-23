import numpy
from Bio import SeqIO
import sys
import pysam
from pstats import Stats
from cProfile import run
import re

### cat sd_0001_PAO1_5k.sam | python get_alleles_from_sam.py sample_name positions.txt

def read_positions(positions_file):
	# all positions are 0-based (Pythonic)
	positions = set(int(pos.strip().split(',')[1]) for pos in open(positions_file, 'r'))
	return positions

def init_array(positions):
	#  1 2 3 4 5 ..
	#A 0 0 0 0 0 ..
	#C 0 0 0 0 0 ..
	#G 0 0 0 0 0 ..
	#T 0 0 0 0 0 ..
	genome_size = max(positions) + 1
	array = numpy.zeros((6,genome_size), dtype=int)
	return array

class parseString(object):

	def __init__(self, ref, string):
		self.ref = ref.upper()
		self.string = string.upper()
		self.types = {'A':0,'G':0,'C':0,'T':0,'-':0,'*':0,'+':0,'X':[]}
		#self.types = {'A':0,'G':0,'C':0,'T':0,'-':[],'*':0,'+':[],'X':[]}
		self.process()
        
	def process(self):
		# remove end of read character
		self.string = self.string.replace('$','')
		while self.string != '':
			if self.string[0] == '^':
				# skip two characters when encountering '^' as it indicates
				# a read start mark and the read mapping quality
				self.string = self.string[2:]
			elif self.string[0] == '*':
				self.types['*'] += 1
				# skip to next character
				self.string = self.string[1:]
			elif self.string[0] in ['.',',']:
				# a reference base
				self.types[self.ref] += 1
				self.string = self.string[1:]
			elif self.string[0] == '+':
				# insertion 
				insertionLength = int(self.string[1]) #will fail if len > 9
				insertionSeq = self.string[2:2+ insertionLength]
				self.types['+'] += 1
				#self.types['+'].append(insertionSeq)
				self.string = self.string[2+insertionLength:]
			elif self.string[0] == '-':
				# deletion
				deletionLength = int(self.string[1])
				deletionSeq = self.string[2:2+deletionLength]
				self.types['-'] += 1
				#self.types['-'].append(deletionSeq)
				self.string = self.string[2+deletionLength:]
			elif self.types.has_key(self.string[0]):
				# one of the four bases
				self.types[self.string[0]] += 1
				self.string = self.string[1:]
			else:
				# unrecognized character
				self.types['X'].append(self.string[0])
				self.string = self.string[1:]
		return

def read_pileup(positions, array):
	base_dict = {'A':0, 'C':1, 'G':2, 'T':3, '+':4, '-':5}
	for line in sys.stdin:
		cols = line.strip('\n').split('\t')
		pos = int(cols[1])
		if pos in positions:
			#print pos, cols[2:5]
			ref = cols[2].upper()
			depth = cols[3]
			counts = parseString(ref, cols[4])
			#print pos, counts.types
			for base in ['A', 'C', 'G', 'T', '+', '-']:
				#print base, counts.types[base]
				array[(base_dict[base]), pos] = counts.types[base]
			

def read_sam(positions, array):
	base_dict = {'A':0, 'C':1, 'G':2, 'T':3}
	alignment = pysam.Samfile('-', 'r')
	for line in alignment:
		# ignore any unmapped reads
		if line.is_unmapped: continue
		# ignore any reads with indels
		if len([f for f, v in line.cigar if f != 0]):
			continue
		read_positions = set(xrange(line.pos, line.aend))
		isec = positions.intersection(read_positions)
		if isec:
			overlap = [(pos, line.seq[pos-line.pos]) for pos in isec]
			#quality = [(pos, ord(line.qual[pos-line.pos])-33) for pos in isec]
			if overlap:
				for each in overlap:
					if each[1] != 'N':
						array[(base_dict[each[1]], each[0])] += 1
	return array

def write_alleles(sample_name, positions, array):
	frag = ''
	cov = int(sys.argv[3])
	chastity_list = []
	for pos in sorted(positions):
		counts = tuple(array[:,pos])

		counts_sort = sorted(counts, reverse=True)
		#print pos, counts_sort
		if sum(counts[:4]) < cov:
			frag += '-'
			continue
		#elif counts[4]:
		#	frag += '-'
		#elif counts[5]:
		#	frag += '-'
		# chastity is greatest / (greatest + second greatest)
		elif all(counts[0] > base for base in counts[1:4]):
			chastity_list.append(float(counts_sort[0]) / sum(counts_sort[:2]))
			frag += 'A'
		elif all(counts[1] > base for base in counts[0:1] + counts[2:4]):
			chastity_list.append(float(counts_sort[0]) / sum(counts_sort[:2]))
			frag += 'C'
		elif all(counts[2] > base for base in counts[0:2] + counts[3:4]):
			chastity_list.append(float(counts_sort[0]) / sum(counts_sort[:2]))
			frag += 'G'
		elif all(counts[3] > base for base in counts[0:3]):
			chastity_list.append(float(counts_sort[0]) / sum(counts_sort[:2]))
			frag += 'T'
		else:
			frag += '-'

	if chastity_list:
		mean_chastity = sum(chastity_list) / (len(frag) - frag.count('-')) * 100.0
	else:
		mean_chastity = 0
	print >>sys.stderr, 'Sample: %s' %(sample_name)
	print >>sys.stderr, 'Total positions: %i' %(len(positions))
	print >>sys.stderr, 'Positions covered: %i (%.2f %%)' %(len(frag) - frag.count('-'), 100 - (float(frag.count('-')) / len(positions) * 100))
	print >>sys.stderr, 'Gaps: %i (%.2f %%)' %(frag.count('-'), float(frag.count('-')) / len(positions) * 100)
	print >>sys.stderr, 'Mean chastity: %.2f %%' %(mean_chastity)
	if mean_chastity < 90:
		print >>sys.stderr, 'CHASTITY WARNING: low chastity can severely affect accuracy of placement'
	print >>sys.stdout, '>%s_cov%i\n%s\n' %(sample_name, cov, frag),

if __name__ == '__main__':
	positions = read_positions(sys.argv[2])
	sample_name = sys.argv[1]
	array = init_array(positions)
	read_pileup(positions, array)
	write_alleles(sample_name, positions, array)