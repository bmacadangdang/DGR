import argparse
from pathlib import Path

import DGRutils as utils
import DGRanalysis as analysis
import DGRfinder as finder
import DGRprocessor as processor

if __name__ == '__main__':
	# Argument parsing
	parser = argparse.ArgumentParser(description='DGR analysis from metagenomic data')
	parser.add_argument('-r1',dest='rawdatafile',required=True,
	                    help='Path to raw sequencing data file in fasta or fastq format')
	parser.add_argument('-g',dest='reference_genome',
	                    help='Path to reference genome file in fasta format')
	parser.add_argument('-m',dest='mode',choices=['dgrFinder', 'dgrAnalysis', 'dgrProcess'],required=True,
	                    help='Mode to run: dgrFinding or dgrAnalysis')
	parser.add_argument('-o',dest='out_dir',default='.',
	                    help='Path to output file directory (default: current directory)')
	parser.add_argument('-v',dest='VR',
						help='Path to file that contains VR sequence in fasta format [Required if dgrAnalysis mode]')
	parser.add_argument('-t',dest='TR',
						help='Path to file that contains TR seqeunce in fasta format [Required if dgrAnalsyis mode]')
	parser.add_argument('--vt_folder',nargs='?',dest='vt_folder', const=True, default=False,
						help='Set to true if VR and TR are a prefix for multiple files within the folder')
	parser.add_argument('-r2',dest='rawdatafile2',
						help='Path to raw sequencing datafile 2 if paired end reads')
	parser.add_argument('--qf', nargs='?',dest='qf', const=True, default=False,
						help='Option to quality filter the reads prior to processing')
	args = parser.parse_args()

	'''
	if args.mode == 'dgrAnalysis':
		#Checking if VR and TR are specified 
		if args.VR == None:
			raise OSError('VR needs to be specified if running dgrAnalysis')
		if not utils.is_fasta_file(args.VR):
			raise OSError('VR file needs to be a fasta file: %s' % (args.VR))
		if args.TR == None:
			raise OSError('TR needs to be specified if running dgrAnalysis')
		if not utils.is_fasta_file(args.TR):
			raise OSError('TR file needs to be a fasta file: %s' % (args.TR))

	'''
	#Need to add handling of gz files

	'''
	#Check raw data in correct format
	if not utils.is_fasta_or_fastq_file(args.rawdatafile):
		raise OSError('Raw data file is not fasta or fastq: %s' % (args.rawdatafile))

	#Check raw data 2 in correct format if defined
	if args.rawdatafile2 is not None:
		if not utils.is_fasta_or_fastq_file(args.rawdatafile2):
			raise OSError('Raw data file is not fasta or fastq: %s' % (args.rawdatafile2))

	if args.qf:
		#Add in quality filtering here
		pass

	'''

	if args.mode == 'dgrAnalysis' or args.mode == 'dgrFinder' or args.mode == 'dgrProcess':
		if not Path(args.out_dir).exists():
			Path(args.out_dir).mkdir(parents=True, exist_ok=True)

	if args.mode == 'dgrAnalysis':
		if args.reference_genome == None:
			raise ValueError('Reference genome is required for DGR Analysis')
		else:
			#Check if reference genome is correct format
			if not utils.is_fasta_file(args.reference_genome):
				raise OSError('Reference genome is not in fasta format: %s' % (args.reference_genome))
		analysis.process_determine_DGR_activity(args.rawdatafile, args.reference_genome, args.VR, args.TR, args.out_dir, args.rawdatafile2, args.vt_folder)
		#analysis.determine_DGR_activity_with_aligners(args.rawdatafile, args.reference_genome, args.VR, args.TR, args.out_dir)
	elif args.mode == 'dgrFinder':
		finder.search_for_DGRs(args.rawdatafile, args.out_dir)
	elif args.mode == 'dgrProcess':
		processor.process_raw_data(args.rawdatafile, args.rawdatafile2, args.out_dir)