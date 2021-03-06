import numpy
from Bio import SeqIO
import sys
import pysam
# import vcf
from pstats import Stats
from collections import OrderedDict
from cProfile import run
from pstats import Stats

USE_CHASTITY = False

### cat sd_0001_PAO1_5k.sam | python get_alleles_from_sam.py sample_name positions.txt vcf_file

def get_args():
	if len(sys.argv) != 4:
		print 'incorrect number of args, useage is get_alleles_from_sam.py <sample_name> <positions.txt> <output_file>'
	else:
		return sys.argv[1], sys.argv[2], sys.argv[3]

def read_positions(positions_txt):
	# all positions are 0-based (Pythonic)
	positions = OrderedDict()
	for line in open(positions_txt, 'r'):
		cols = line.strip().split(',')
		if cols[0] not in positions.keys():
			positions[cols[0]] = set()
			positions[cols[0]].add(int(cols[1]) - 1)
		else:
			positions[cols[0]].add(int(cols[1]) - 1)
	return positions

def init_array(refs, position_file):
	#  1 2 3 4 5 ...
	#A 0 0 0 0 0 ...
	#C 0 0 0 0 0 ...
	#G 0 0 0 0 0 ...
	#T 0 0 0 0 0 ...
	array = {}
	positions = read_positions(position_file)
	for ref in positions.keys():
		array[ref] = numpy.zeros((4, max(positions[ref]) + 1), dtype=int)
	return array, positions

def read_sam(position_file):
	base_dict = {'A':0, 'C':1, 'G':2, 'T':3}
	alignment = pysam.Samfile('-', 'r')
	array, positions = init_array(alignment.references, position_file)
	for line in alignment:
		# ignore any unmapped reads
		if line.is_unmapped:
			continue
		chrom = alignment.getrname(line.tid)
		read_positions = set(xrange(line.pos, line.aend))
		try:
			isecs = positions[chrom].intersection(read_positions)
		except KeyError:
			continue
		if isecs:
			#overlap = [(pos, line.seq[pos-line.pos]) for pos in isec]
			#quality = [(pos, ord(line.qual[pos-line.pos])-33) for pos in isec]
			aligned_pairs = dict((ref, query) for (query, ref) in line.get_aligned_pairs())
			for isec in isecs:
				if aligned_pairs[isec]:
					read_base = line.seq[aligned_pairs[isec]]
					if read_base != 'N':
							array[chrom][(base_dict[read_base], isec)] += 1
	return array, positions

def write_alleles(array, positions, sample_name, output_file):
	
	
	frag = ''
	chastity_list = []
	gaps = 0
	#indel_positions = [(position.CHROM, position.POS - 1) for position 
	#			   in list(vcf.Reader(open(sys.argv[3], 'r'))) if position.is_indel]
	indel_positions = []
	for chrom in positions.keys():
		for pos in sorted(positions[chrom]):
			if (chrom, pos) in indel_positions:
				gaps += 1
				frag += '-'
			else:
				counts = tuple(array[chrom][:,pos])
				counts_sort = sorted(counts, reverse=True)
				if all(counts[0] > base for base in counts[1:4]):
					# chastity is greatest / (greatest + second greatest)
					if USE_CHASTITY:
						chastity_list.append(float(counts_sort[0]) / 
							     sum(counts_sort[:2]))
					frag += 'A'
				elif all(counts[1] > base for base in counts[0:1] + counts[2:4]):
					if USE_CHASTITY:
						chastity_list.append(float(counts_sort[0]) / 
							     sum(counts_sort[:2]))
					frag += 'C'
				elif all(counts[2] > base for base in counts[0:2] + counts[3:4]):
					if USE_CHASTITY:
						chastity_list.append(float(counts_sort[0]) / 
							     sum(counts_sort[:2]))
					frag += 'G'
				elif all(counts[3] > base for base in counts[0:3]):
					if USE_CHASTITY:
						chastity_list.append(float(counts_sort[0]) / 
							     sum(counts_sort[:2]))
					frag += 'T'
				else:
					gaps += 1
					frag += '-'

	total = sum([len(positions[chrom]) for chrom in positions.keys()])
	print 'Sample: %s' %(sample_name)
	print 'Total positions: %i' %(total)
	print 'Gaps: %i' %(gaps)
	print 'Positions covered: %.2f %%' %(100 - (float(gaps) / total * 100.0))
	print 'Sample: %s' %(sample_name)
	if USE_CHASTITY:
		mean_chastity = (sum(chastity_list) / (len(frag) - gaps)) * 100.0
		print 'Mean chastity: %.2f %%' %(mean_chastity)
		if mean_chastity < 90:
			print 'CHASTITY WARNING: Mixed samples can severely affect accuracy of placement'

	with open(output_file, 'w') as file_out:
		print >>file_out, '>%s_new\n%s\n' %(sample_name, frag),


if __name__ == '__main__':
	sample_name, position_file, output_file = get_args()
	array, positions = read_sam(position_file)
	# run('write_alleles(array, positions)', 'stats')
	#stats = Stats('stats')
	write_alleles(array, positions, sample_name, output_file)



