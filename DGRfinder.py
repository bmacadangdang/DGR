from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline

import numpy as np
import csv, glob, math, datetime, os, time
import pandas as pd
from pathlib import Path

import DGRutils as utils

def create_sliding_windows(input_file, output_file, window_size=200, step_size=50):

	output = open(output_file, 'w')

	records = SeqIO.parse(input_file, 'fasta')
	for item in records:
		contig = item.name
		seq = str(item.seq)

		num_chunks = ((len(seq)-window_size)//step_size)+1

		for x in range(0, num_chunks*step_size, step_size):
			y = x + window_size
			dna = seq[x:x+window_size]
			output.write('>%s_%i_to_%i\n%s\n' % (contig, x, y, dna))

	output.close()

#Input: two potential DNA sequences that could be VR or TR (from Blast)
#Output: False is no VR/TR found, otherwise return [VR, TR]
def determine_VR_TR_pair(seq1, seq2):
	if len(seq1) != len(seq2):
		return False

	if len(seq1) < 60:
		return False

	bases = ['A', 'T', 'C', 'G']
	seq1_counts = [0, 0, 0, 0]
	seq2_counts = [0, 0, 0, 0]
	seq1_mismatches = [0, 0, 0, 0]
	seq2_mismatches = [0, 0, 0, 0]

	seq1 = seq1.upper()
	seq2 = seq2.upper()

	for i in range(len(seq1)):
		if seq1[i] in bases and seq2[i] in bases:
			pos1 = bases.index(seq1[i])
			pos2 = bases.index(seq2[i])

			seq1_counts[pos1] += 1
			seq2_counts[pos2] += 1

			if seq1[i] != seq2[i]:
				seq1_mismatches[pos1] += 1
				seq2_mismatches[pos2] += 1

	posA = bases.index('A')
	posT = bases.index('T')

	varA1, varT1, varA2, varT2 = [0] * 4

	if seq1_mismatches[posA] > 0:
		varA1 = seq1_mismatches[posA] / sum(seq1_mismatches)
	if seq1_mismatches[posT] > 0:
		varT1 = seq1_mismatches[posT] / sum(seq1_mismatches)
	if seq2_mismatches[posA] > 0:
		varA2 = seq2_mismatches[posA] / sum(seq2_mismatches)
	if seq2_mismatches[posT] > 0:
		varT2 = seq2_mismatches[posT] / sum(seq2_mismatches)

	if seq1_mismatches[posA] >= 5 or seq1_mismatches[posT] >= 5 or seq2_mismatches[posA] >= 5 or seq2_mismatches[posT] >= 5:
		if varA1 > 0.8 or varT1 > 0.8:
			return [seq2, seq1]
		elif varA2 > 0.8 or varT2 > 0.8:
			return [seq1, seq2]
		else:
			return False
	else:
		return False

def create_RT_hmmprofile(hmmfile=None):
	hmmfolder = 'hmmprofiles'
	#Look for hmmfile, if not there, then try default
	if hmmfile != None:
		if Path(hmmfile).exists():
			hmm = hmmfile
		else:
			if Path('hmmprofiles/DGR_RTs_MSA.fa').exists():
				hmm = 'hmmprofiles/DGR_RTs_MSA.fa'
			else:
				raise OSError('HMM Profile file: %s does not exist.  Cannot find default file either.')
	else:
		if Path('hmmprofiles/DGR_RTs_MSA.fa').exists():
			hmm = 'hmmprofiles/DGR_RTs_MSA.fa'
		else:
			raise OSError('Default HMM profile hmmprofiles/DGR_RTs_MSA.fa does not exist.  Please download from github')

	hmm_files = ['hmmprofiles/DGR_RTs_MSA.h3m', 'hmmprofiles/DGR_RTs_MSA.h3i', 'hmmprofiles/DGR_RTs_MSA.h3f', 'hmmprofiles/DGR_RTs_MSA.h3p', 'hmmprofiles/DGR_RTs_MSA']
	for hf in hmm_files:
		if Path(hf).exists():
			os.system('rm %s' % (hf))
	os.system('hmmbuild hmmprofiles/DGR_RTs_MSA %s' % (hmm))
	os.system('hmmpress hmmprofiles/DGR_RTs_MSA')


def search_for_DGRs(contigfile, output_folder, hmmfile=None):
	create_RT_hmmprofile(hmmfile)

	#Create temp_directory in output_folder to store information
	temp_folder = '%s/temp' % (output_folder)
	os.system('rm -rf %s' % (temp_folder))
	time.sleep(3)
	Path(temp_folder).mkdir(parents=True, exist_ok=True)

	#Keep track of files that need to be deleted at the end of the analysis
	temp_files = []
	temp_blast_db = []

	#Get the name of the rawdata file and remove the file extension
	rawdata_name = '.'.join(contigfile.split('/')[-1].split('.')[:-1])

	#Create output folders
	#OO_Raw will contain all raw VR/TR pairs
	#01_Filtered will contain all unique VR/TR pairs
	#02_Processed will contain all unique VR/TR pairs that are contained within an ORF
	Path('%s/00_Raw' % (output_folder)).mkdir(parents=True, exist_ok=True)
	Path('%s/01_Filtered' % (output_folder)).mkdir(parents=True, exist_ok=True)
	Path('%s/02_Processed' % (output_folder)).mkdir(parents=True, exist_ok=True)

	#If the rawdata file is in fastq, it needs to be converted to fasta for blast, delete that file at the end
	#Otherwise just use the fastafile supplied and it does not need ot be deleted
	if utils.is_fastq_file(contigfile):
		print('Converting rawdata from FASTQ to FASTA')
		rawdata_fasta = '%s/%s.fa' % (temp_folder, rawdata_name)
		temp_files.append(rawdata_fasta)
		with open(rawdata_fasta, 'w') as f:
			data = SeqIO.parse(contigfile, 'fastq')
			i = 1
			for seq in data:
				f.write('>Sequence%i\n%s\n' % (i, str(seq.seq)))
				i += 1
	else:
		rawdata_fasta = contigfile

	#Format the reference genome file for easier searching later
	formatted_contigs = '%s/%s-formatted_contigs.fasta' % (temp_folder, rawdata_name)
	temp_files.append(formatted_contigs)
	utils.format_ref_genome_file(rawdata_fasta, formatted_contigs)


	#Find all ORFs in the rawdata and convert to protein sequences for input into hmmscan, using Prodigal anonymous mode (as recommended by software for metagenomes)
	protein_coordinates = '%s/%s-genes.gff' % (temp_folder, rawdata_name)
	protein_seqs = '%s/%s-protein_seqs.fa' % (temp_folder, rawdata_name)
	temp_files.append(protein_coordinates)
	temp_files.append(protein_seqs)
	print('Predicting genes with Prodigal.')
	os.system('prodigal -i %s -o %s -a %s -f gff -p meta -q' % (formatted_contigs, protein_coordinates, protein_seqs))

	#Using tblout is likely a better output file format, but will require further testing
	print('Predicting RTs with HMMER')
	hmmout = '%s/RT_hits.txt' % (output_folder)
	temp_files.append(hmmout)
	#os.system('hmmscan --domtblout %s/RT_hits_domtblout.txt hmmprofiles/DGR_RTs_MSA %s' % (output_folder, protein_seqs))
	os.system('hmmscan --tblout %s -o /dev/null hmmprofiles/DGR_RTs_MSA %s' % (hmmout, protein_seqs))

	#Filter RT hits (score > 20, protein length 200-550 AA long)
	#RT_hits should be a list of RTs in format: [[contig#, start, end], ]
	RT_hits = utils.parse_hmmoutput_from_prodigal(hmmout, protein_coordinates, score_cutoff=20, aa_length_min=200, aa_length_max=550)

	#Extract 20kb up and downstream of each RT in order to analyze for VR/TR pairs
	extract_window = 20000
	for RT_hit in RT_hits:
		RT_area_file = '%s/%s-RT_area.fa' % (temp_folder, rawdata_name)
		temp_files.append(RT_area_file)
		utils.extract_RT_regions(formatted_contigs, RT_hit, RT_area_file, window=extract_window)

		RT_contig = RT_hit[0]
		RT_start = int(RT_hit[1])
		RT_end = int(RT_hit[2])

		RT_area_start = RT_start - extract_window
		if RT_area_start < 0:
			RT_area_start = 0

		#Create sliding windows for each RT_area
		sliding_window_file = '%s/%s-sliding_window.fa' % (temp_folder, rawdata_name)
		temp_files.append(sliding_window_file)
		create_sliding_windows(RT_area_file, sliding_window_file)

		#Create blast database of each 200 bp window in order to blast against each other
		output_options = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'sseq', 'qseq']
		blast_out_str = ' '.join(output_options)
		blast_sliding_window_hits = '%s/%s-sliding_window-blast_output.txt' % (temp_folder, rawdata_name)
		temp_files.append(blast_sliding_window_hits)
		sliding_window_blastdb = '%s/%s_sliding_window_blastdb' % (temp_folder, rawdata_name)
		temp_blast_db.append(sliding_window_blastdb)

		cline = NcbimakeblastdbCommandline(dbtype='nucl', input_file=sliding_window_file, out=sliding_window_blastdb)
		cline()
		time.sleep(3)
		print('Blast for potential VR/TR pairs surrounding RT')
		cline = NcbiblastnCommandline(out=blast_sliding_window_hits, db=sliding_window_blastdb, query=sliding_window_file, outfmt='6 %s' % (blast_out_str) , word_size=8, reward=1, penalty=-1, evalue=1e-5, gapopen=6, gapextend=6, perc_identity=50)
		cline()

		#Parse blast output to determine if hits are true VR/TR pairs
		VR_raw_outfile = '%s/00_Raw/%s-RT%s_%i_%i_VRrawseqs.fa' % (output_folder, rawdata_name, RT_contig, RT_start, RT_end)
		TR_raw_outfile = '%s/00_Raw/%s-RT%s_%i_%i_TRrawseqs.fa' % (output_folder, rawdata_name, RT_contig, RT_start, RT_end)

		mismatch = output_options.index('mismatch')
		query_seq = output_options.index('qseq')
		sbjct_seq = output_options.index('sseq')
		query_name = output_options.index('qseqid')
		sbjct_name = output_options.index('sseqid')
		qstart = output_options.index('qstart')
		qend = output_options.index('qend')
		sstart = output_options.index('sstart')
		send = output_options.index('send')

		num_hits = 0
		with open(blast_sliding_window_hits, 'r') as blast_hits, open(VR_raw_outfile, 'w') as VR_writer, open(TR_raw_outfile, 'w') as TR_writer:
			blast_reader = csv.reader(blast_hits, delimiter='\t')
			for line in blast_reader:
				if int(line[mismatch]) > 0:
					pair = determine_VR_TR_pair(line[query_seq], line[sbjct_seq])
					if pair != False:
						num_hits += 1
						seq_length = len(line[query_seq])
						if pair[0] == line[query_seq]:
							VR_name = line[query_name]
							TR_name = line[sbjct_name]
							vstart_pos = int(line[qstart])
							if int(line[qend]) < vstart_pos:
								vstart_pos = int(line[qend])
							tstart_pos = int(line[sstart])
							if int(line[send]) < tstart_pos:
								tstart_pos = int(line[send])
						else:
							VR_name = line[sbjct_name]
							TR_name = line[query_name]
							tstart_pos = int(line[qstart])
							if int(line[qend]) < tstart_pos:
								tstart_pos = int(line[qend])
							vstart_pos = int(line[sstart])
							if int(line[send]) < vstart_pos:
								vstart_pos = int(line[send])

						vstart = RT_area_start + int(VR_name.split('_')[-3]) + vstart_pos-1
						tstart = RT_area_start + int(TR_name.split('_')[-3]) + tstart_pos-1

						VR_name = '%s-%s-%i_to_%i' % (rawdata_name, RT_contig, vstart, vstart + seq_length)
						TR_name = '%s-%s-%i_to_%i' % (rawdata_name, RT_contig, tstart, tstart + seq_length)

						VR_writer.write('>%s\n%s\n' % (VR_name, pair[0]))
						TR_writer.write('>%s\n%s\n' % (TR_name, pair[1]))

		#Combine overlapping VRs and TRs, there may be more than one
		if num_hits > 0:
			non_overlapping_coordinates_vr = []
			non_overlapping_coordinates_tr = []
			vr_hits = SeqIO.parse(VR_raw_outfile, 'fasta')
			tr_hits = SeqIO.parse(TR_raw_outfile, 'fasta')
			vr = next(vr_hits)
			parts = vr.name.split('-')
			coordinates = parts[-1].split('_')
			start = int(coordinates[0])
			end = int(coordinates[-1])
			non_overlapping_coordinates_vr.append([start, end])
			tr = next(tr_hits)
			parts = tr.name.split('-')
			coordinates = parts[-1].split('_')
			start = int(coordinates[0])
			end = int(coordinates[-1])
			non_overlapping_coordinates_tr.append([start, end])
			for vr,tr in zip(vr_hits, tr_hits):
				parts = vr.name.split('-')
				coordinates = parts[-1].split('_')
				vr_start = int(coordinates[0])
				vr_end = int(coordinates[-1])
				found = False
				for vr_noc,tr_noc in zip(non_overlapping_coordinates_vr, non_overlapping_coordinates_tr):
					if vr_start >= vr_noc[0] and vr_start <= vr_noc[1]:
						found = True
						if vr_end > vr_noc[1]:
							tr_noc[1] += (vr_end - vr_noc[1])
							vr_noc[1] = vr_end
					if vr_end >= vr_noc[0] and vr_end <= vr_noc[1]:
						found = True
						if vr_start < vr_noc[0]:
							tr_noc[0] -= (vr_noc[0] - vr_start)
							vr_noc[0] = vr_start
				if not found:
					non_overlapping_coordinates_vr.append([vr_start, vr_end])
					parts = tr.name.split('-')
					coordinates = parts[-1].split('_')
					tr_start = int(coordinates[0])
					tr_end = int(coordinates[-1])
					non_overlapping_coordinates_tr.append([tr_start, tr_end])


			#Save all non-overlapping VR/TR pairs to a file
			VR_unique_seqs = '%s/01_Filtered/%s-%s_%i_%i-VR_uniques.fa' % (output_folder, rawdata_name, RT_contig, RT_start, RT_end)
			TR_unique_seqs = '%s/01_Filtered/%s-%s_%i_%i-TR_uniques.fa' % (output_folder, rawdata_name, RT_contig, RT_start, RT_end)
			with open(VR_unique_seqs, 'w') as vr_writer, open(TR_unique_seqs, 'w') as tr_writer:
				for vr,tr in zip(non_overlapping_coordinates_vr, non_overlapping_coordinates_tr):
					sequences = SeqIO.parse(formatted_contigs, 'fasta')
					for sequence in sequences:
						if sequence.name == RT_contig:
							vr_seq = str(sequence.seq)[vr[0]:vr[1]]
							tr_seq = str(sequence.seq)[tr[0]: tr[1]]
							vr_writer.write('>%s-%s-%i_to_%i\n%s\n' % (rawdata_name, RT_contig, vr[0], vr[1], vr_seq))
							tr_writer.write('>%s-%s-%i_to_%i\n%s\n' % (rawdata_name, RT_contig, tr[0], tr[1], tr_seq))

			#Determine if VR exists within an ORF
			VR_final = '%s/02_Processed/%s-%s_%i_%i-VR.fasta' % (output_folder, rawdata_name, RT_contig, RT_start, RT_end)
			TR_final = '%s/02_Processed/%s-%s_%i_%i-TR.fasta' % (output_folder, rawdata_name, RT_contig, RT_start, RT_end)
			vrs = SeqIO.parse(VR_unique_seqs, 'fasta')
			trs = SeqIO.parse(TR_unique_seqs, 'fasta')
			for vr,tr in zip(vrs, trs):
				coordinates = vr.name.split('-')[-1].split('_')
				vr_start = int(coordinates[0])
				vr_end = int(coordinates[-1])
				with open(protein_coordinates, 'r') as pc, open(VR_final, 'w') as vr_writer, open(TR_final, 'w') as tr_writer:
					reader = csv.reader(pc, delimiter='\t')
					for line in reader:
						if line[0][0] != '#':
							if line[0] == RT_contig:
								protein_start = int(line[3])
								protein_end = int(line[4])
								if (vr_start >= protein_start and vr_start <= protein_end) or (vr_end >= protein_start and vr_end <= protein_end):
									vr_writer.write('>%s\n%s\n' % (vr.name, str(vr.seq)))
									tr_writer.write('>%s\n%s\n' % (tr.name, str(tr.seq)))



