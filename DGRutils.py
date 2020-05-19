from Bio import SeqIO
from Bio.Seq import Seq

import os, glob, time, csv

def is_fasta_file(filename):
	with open(filename, 'r') as f:
		try:
			fasta = SeqIO.parse(f, 'fasta')
			return any(fasta)
		except ValueError as error:
			return False

def is_fastq_file(filename):
	with open(filename, 'r') as f:
		try:
			fasta = SeqIO.parse(f, 'fastq')
			return any(fasta)
		except ValueError as error:
			return False

def is_fasta_or_fastq_file(filename):
	return any([is_fastq_file(filename), is_fasta_file(filename)])

'''	
	Input: target = sequence to find
		  reference_filename = filename containing reference sequence
	Output: the start and ending coordinates of target within reference sequence as well as the sequence number (for files that have multiple sequences) and the target sequence itself (reverse complemented if needed)
'''
def extract_coordinates(target, reference_filename):

	with open(reference_filename, 'r') as f:
		ref = SeqIO.parse(f, 'fasta')
		num = 0
		rev = False
		rev_target = str(Seq(target).reverse_complement())
		for sequence in ref:
			if sequence.seq.find(target) != -1 or sequence.seq.find(rev_target) != -1:
				start = sequence.seq.find(target)
				if start == -1:
					start = sequence.seq.find(rev_target)
					rev = True
				end = start + len(target)
				if end >= len(sequence.seq):
					end = len(sequence.seq) - 1
				if rev:
					return (start, end, num, rev_target)
				else:
					return (start, end, num, target)
			num += 1
	raise ValueError('target: %s could not be found in reference: %s' % (target, reference_filename))

# Wrapper for extract_coordinates, may expand this in future
def extract_sequence(target_filename, reference_filename):
	with open(target_filename, 'r') as f:
		results = SeqIO.parse(target_filename, 'fasta')
		target = str(next(results).seq)
		return extract_coordinates(target, reference_filename)

# Function that creates contig names that are sequentially ordered for easier reference later
def format_ref_genome_file(inputfile, outputfile):
	with open(outputfile, 'w') as f:
		data = SeqIO.parse(inputfile, 'fasta')
		i = 0
		for sequence in data:
			f.write(">Contig%i\n%s\n" % (i, str(sequence.seq)))
			i += 1

def cleanup_files(file_list):
	for temp in file_list:
		try:
			os.remove(temp)
		except:
			print('Could not remove %s' % (temp))

def cleanup_blast_db(file_list):
	for temp in file_list:
		blast_extensions = ['ndb', 'nhr', 'nin', 'not', 'nsq', 'ntf', 'nto']
		for ext in blast_extensions:
			try:
				os.remove('%s.%s' % (temp, ext))
			except:
				print('Could not remove %s.%s' % (temp, ext))

def cleanup(temp_files, temp_dbs):
	time.sleep(10)
	cleanup_files(temp_files)
	cleanup_blast_db(temp_dbs)

#Input: hmm output file in tblout format, protein sequence gff file from prodigal
#Output: List of contigs that meet criteria: example [['Contig#', start, end],]
def parse_hmmoutput_from_prodigal(inputfilename, protein_sequence_file, score_cutoff=20, aa_length_min=200, aa_length_max=550):
	#First read through protein sequence gff file to store information
	prodigal_info = {}
	with open(protein_sequence_file, 'r') as f:
		reader = csv.reader(f, delimiter='\t')
		for line in reader:
			if line[0][0] != '#':
				contig = line[0]
				start = int(line[3])
				end = int(line[4])
				protein_id = contig + '_' + line[8].split(';')[0].split('=')[1].split('_')[1]
				prodigal_info[protein_id] = [contig, start, end]

	#Next parse the hmmoutput file which is space delimited
	output = []
	with open(inputfilename, 'r') as f:
		for line in f:
			info = line.split()
			if info[0][0] != '#':
				score = float(info[5])
				if score >= score_cutoff:
					protein_id = info[2]
					contig, start, end = prodigal_info[protein_id]
					aa_length = abs(end - start) / 3
					if aa_length >= aa_length_min and aa_length <= aa_length_max:
						output.append([contig, start, end])

	return output

def extract_RT_regions(genome_file, RT, outputfile, window=20000):
	with open(outputfile, 'w') as writer:
		contig_name = RT[0]
		start = int(RT[1]) - window
		end = int(RT[2]) + window
		seq = extract_genomic_region(genome_file, start, end, contig_name)
		writer.write('>%s-RT_area\n%s\n' % (contig_name, seq))


#Given coordinates of a contig, extract the DNA sequence of between start and end
def extract_genomic_region(genome_file, start, end, contig_name):
	contigs = SeqIO.parse(genome_file, 'fasta')
	for contig in contigs:
		if contig.name == contig_name:
			if start < 0:
				start = 0
			if end > len(contig.seq):
				end = len(contig.seq)
			return str(contig.seq)[start:end]
	raise ValueError('Could not find contig: %s' % (contig_name))
