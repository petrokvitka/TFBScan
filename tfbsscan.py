"""
TFBSscan.py produces the data to be used to join the footprint and motif information across genome

@author: Anastasiia Petrova
@contact: anastasiia.petrova(at)mpi-bn.mpg.de

"""

import argparse
import sys
import os
import re
import time
import multiprocessing
import logging
import subprocess
from Bio import SeqIO
import Bio.SeqIO.FastaIO as bio
import MOODS.scan
import MOODS.tools
import MOODS.parsers

logger = logging.getLogger('mytool')
logger.setLevel(logging.INFO)

formatter = logging.Formatter("%(asctime)s : %(message)s", "%Y-%m-%d %H:%M")

fh = logging.FileHandler('final_log.txt')
fh.setLevel(logging.INFO)
fh.setFormatter(formatter)
logger.addHandler(fh)

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)

#catch all the information about input and output files as well as information on the used tool (fimo or moods)
def parse_args():
	parser = argparse.ArgumentParser(prog = 'mytool_ver6', description = textwrap.dedent('''

		This script takes a list of motifs loaded from jaspar.genereg.net as a combined text file in MEME or .PFM format, a genome file in FASTA format and optionaly a .bed file (the one you want to be merged with the whole genome file) as input. If you want to merge a .bed file with the whole genome file, please enter --bed_file or -b bevor your .bed file. The tool will provide a file output_merge.fa, which you can also use for your research later on. If you already have a merged file, please give this one as genome file input. If there are several motifs in the input file, the tool will create a separate output file for each motif. Choose if you want to use fimo or moods with --use, this script uses by default fimo. Please note that the tool can not provide the calculation of q values with fimo due to the text mode that fimo needs to use. The tool sends merged genome file and motifs to fimo or moods, saves the sorted output for each of the given motifs as moods/fimo_output_[alternate name and id of the motif].txt in the output directory, then calculates the start and the end as real positions on the chromosom and writes this information in the ouput files. The columns in the output file are: chromosom, start, end, the name and score of TF. If a .bed file was given as input, the tool will also add the additional columns from it to the output. If the output file is empty, there were no machtes within given genome regions. Please note, if you want to have all intermediate output files, enter --clean nothing

		'''), epilog='That is what you need to make this script work for you. Enjoy it')

	required_arguments = parser.add_argument_group('required arguments')

	required_arguments.add_argument('-m', '--motifs', help='file in MEME format with mofits loaded from jaspar.genereg.net')
	required_arguments.add_argument('-g', '--genome', help='a whole genome file or regions of interest in FASTA format to be scanned with motifs')

	#all other arguments are optional
	parser.add_argument('-o', '--output_directory',  default='output', const='output', nargs='?', help='output directory, default ./output/')
	parser.add_argument('-b', '--bed_file',  nargs='?', help='a .bed file to be merged with the whole genome file to find regions of interest')
	parser.add_argument('--use', '--use_tool', default='fimo', const='fimo', nargs='?', choices=['fimo', 'moods'], help='choose the tool to work with, default tool is fimo')
	parser.add_argument('--clean', nargs='*', choices=['nothing', 'all', 'cut_motifs', 'fimo_output', 'merge_output', 'moods_output'], dest='cleans', help='choose the files you want to delete from the output directory, the default is deleting all the temporary files from the directory')
	parser.add_argument('--fimo', help='enter additional options for fimo using = inside "", for example fimo="--norc" to not score the reverse complement DNA strand. By default the --text mode is used and the calculation of the q values due to the --text mode is not possible')
	parser.add_argument('--cores', type=int, help='number of cores allowed to use by this tool, by default the tool uses 2 cores', default=2)
	parser.add_argument('-p', '--p_value', type=float, help='enter the p value, the default p value is 0.0001. Please note that if you enter the p value using --fimo="--thresh ..." as well, the one within --fimo call will be used', default=0.0001)
	parser.add_argument('--resolve_overlaps', action='store_true', help='delete overlaps with greater p value, by default no overlaps are deleted')
	parser.add_argument('--hide_info', action='store_true', help='while working with data write the information only into ./log.txt')
	parser.add_argument('--moods_bg', nargs='+', type=float, help='set the bg for moods, by default moods uses the bg is 0.25 0.25 0.25 0.25')

	args = parser.parse_args()
	return args

def check_directory(directory):
	if not os.path.exists(directory):
		os.makedirs(directory)
		logger.info('a new directory ' + directory + ' was created')
					
#merge the whole genome with the regions mentioned in .bed file
def merge(genome, bed_file, output_directory):
	logger.info('the merging of files ' + genome + ' and ' + bed_file + ' will end soon, the result file is output_merge.fa')
	output_merge = os.path.join(output_directory, "output_merge.fa")
	os.system("bedtools getfasta -fi " + genome + " -bed " + bed_file + " -fo " + output_merge)
	return output_merge

#split the motifs each in other file
def split_motifs(motifs, output_directory, usage):
	logger.info("the file with motifs " + motifs + " will be checked for motifs and if needed splitted in files each containing only one motif")
	
	first_line = subprocess.getoutput("head -1 " + motifs) #find the first line of the input file
	
	if usage == "moods":
		if first_line.startswith(">"):
			#the motif file probably has the .pfm format, try to read and split it
			splitted_motifs = read_pfm(motifs, output_directory)
		else: #maybe the file with motifs is in MEME format, so try to convert it
			logger.info("the input file has not the expected format, I will try to convert it to .pfm format")
			splitted_motifs = convert_meme_to_pfm(motifs, output_directory)

	elif usage == "fimo":
		if first_line.startswith("MEME version"):
			#the motifs file has probably the MEME format, try to read and split it
			splitted_motifs = read_meme(motifs, output_directory)
			
			#if the there was a convertion before, delete all the .pfm files as we don't need them
			for filename in os.listdir(output_directory):
				if filename.endswith(".pfm"):
					remove_file(os.path.join(output_directory, filename))

		else: #maybe the file with motifs is in .pfm format, so try to convert is
			logger.info("the input file has not the expected format, I will try to convert it to MEME format")
			splitted_motifs = convert_pfm_to_meme(motifs, output_directory)

	return splitted_motifs

def read_pfm(motifs, output_directory):
	splitted_motifs = [] #to save the names of files after splitting
	motif = [] #to save the motif itself, which will be written to the file

	with open(motifs) as read_file:
		lines = 0
		for line in read_file:
			#as the motif has first line with the name and 4 lines with information, if the 5th line is something else than the name of the next motif, the exit will be forced
			if lines == 5 and not line.startswith(">"):
				logger.info('please make sure that the file with motifs has a right format and the number of lines is right in the motif file')
				sys.exit()
			else:
				if line.startswith(">"):
					if 'written_file' in locals():
						written_file.write(''.join(motif))
						motif = []
						lines = 0
						written_file.close()

					motif_alternate_name = check_name(re.split(' ', line)[1].rstrip())
					motif_id = re.split(' ', line[1:])[0] #[1:] meands do not use the first character
					motif_name = os.path.join(output_directory, motif_alternate_name + '_' + motif_id + '.pfm')
							
					splitted_motifs.append(motif_name)
					written_file = open(motif_name, 'w')
						
			if lines >= 1 and lines <= 4: #one motif has 5 lines, the first consists the name, the next 4 - the information we need to proceed the data within moods
				motif.append(line)
					
			lines = lines + 1
	written_file.write(''.join(motif))
	written_file.close()

	return splitted_motifs

def read_meme(motifs, output_directory):
	splitted_motifs = [] #to save the names of files after splitting
	motif = [] #to save the motif itself, which will be written to the file
	head = [] #define a list for header, as fimo needs a header in each motif file it proceedes

	with open(motifs) as read_file:
		lines = 0
		for line in read_file:
			#make the head part
			if lines <= 8:
				if lines == 0 and not line.startswith("MEME version"):
					logger.info('please make sure that the file with motifs has a right format and the number of lines is right in the motif file')
					sys.exit()

				head.append(line)
			else:
				#search for motifs and save each to another file
				if line.startswith("MOTIF"):
					if 'written_file' in locals():
						written_file.write(''.join(motif))
						motif = []
						written_file.close()

					#the alternate name will be checked for validity and the invalid chars will be replaced with '_'
					if len(re.split(' ', line.rstrip())) == 3: #in the input motif file the motif id and the alternate name are splitted using the tab, otherwise they are splitted using _, but we do not want to change it if so
						motif_alternate_name = check_name(re.split(' ', line)[2].rstrip())
						motif_id = re.split(' ', line)[1]
						motif_name = os.path.join(output_directory, motif_alternate_name + '_' + motif_id + '.meme')
					else: 
						motif_alternate_name = check_name(re.split(' ', line)[1].rstrip())
						motif_name = os.path.join(output_directory, motif_alternate_name + '.meme')

					#make a list with all the motif names to know which files to iterate when fimo is called
					splitted_motifs.append(motif_name)

					written_file = open(motif_name, 'w')
					written_file.write(''.join(head))

				motif.append(line)

			lines = lines + 1
					
		#write the last motif
		written_file.write(''.join(motif))
		written_file.close()

	read_file.close()

	return splitted_motifs
	
def convert_meme_to_pfm(motifs, output_directory):
	#i can only convert the file to pfm if the motifs file is in MEME format

	splitted_motifs = [] #to save the names of files after splitting
	rows = [[] for row in range(4)]

	with open(motifs) as read_file:
		lines = 0
		for line in read_file:
			if lines == 0 and not line.startswith("MEME version"):
				logger.info('please make sure that the file with motifs has a right format and the number of lines is right in the motif file')
				sys.exit()
			else:
				#search for motifs and save each to another file
				if line.startswith("MOTIF"):
					
					if 'written_file' in locals():
						for row in rows:
							written_file.write('\t'.join(row) + '\n')

						rows = [[] for row in range(4)]

						written_file.close()
					
					#the alternate name will be checked for validity and the invalid chars will be replaced with '_'
					if len(re.split(' ', line.rstrip())) == 3: #in the input motif file the motif id and the alternate name are splitted using the tab, otherwise they are splitted using _, but we do not want to change it if so
						motif_alternate_name = check_name(re.split(' ', line)[2].rstrip())
						motif_id = re.split(' ', line)[1]
						motif_name = os.path.join(output_directory, motif_alternate_name + '_' + motif_id + '.pfm')

					else: 
						motif_alternate_name = check_name(re.split(' ', line)[1].rstrip())
						motif_name = os.path.join(output_directory, motif_alternate_name + '.pfm')
					
					#make a list with all the motif names to know which files to iterate when fimo is called
					splitted_motifs.append(motif_name)

					written_file = open(motif_name, 'w')
					
				elif line.startswith("letter-probability matrix"):
					columns = int(re.split(' ', re.split('w= ', line)[1])[0]) #find the number of columns from the line out of the motifs file
					nsites = int(re.split(' ', re.split('nsites= ', line)[1])[0]) #find the nsites to count the frequency count for .pfm file

				elif line.startswith(' '): #each line with information about frequency starts in MEME format with ' '
					for i in range(len(rows)):
						rows[i].append(str(round(float(re.findall(r'\S+', line)[i])*nsites))) #split the line, do not mention how much whitespaces are in between, multiply it with nsites and save it to the corresponding row

			lines = lines + 1
					
		#write the last motif
		for row in rows:
			written_file.write('\t'.join(row) + '\n')		

		written_file.close()
	read_file.close()

	return splitted_motifs

def convert_pfm_to_meme(motifs, output_directory):
	#i can only convert the file to meme, if motifs file is in .pfm format

	#first we need to split the pfm motifs as the jaspar2meme does not work with the files containing several motifs, but with the directory consisting files each with only one motif in pfm format
	pfm_motifs = read_pfm(motifs, output_directory)

	converted_name = os.path.join(output_directory, "converted_motifs.meme")

	os.system("jaspar2meme -pfm " + output_directory + " > " + converted_name)

	#need to call split motifs for meme file
	splitted_motifs = split_motifs(converted_name, output_directory, "fimo")

	remove_file(converted_name)
	return splitted_motifs

#if there are chars that are not allowed, they will be replaced with '_', to the possibly invalid names there will be added '_' at the beginning of the name
def check_name(name_to_test):
	badchars= re.compile(r'[^A-Za-z0-9_. ]+|^\.|\.$|^ | $|^$')
	badnames= re.compile(r'(aux|com[1-9]|con|lpt[1-9]|prn)(\.|$)')

	#replace all the chars that are not allowed with '_'
	name = badchars.sub('_', name_to_test)

	#check for the reserved by the os names
	if badnames.match(name):
		name = '_' + name
	return name

#use fimo to make a file
def call_fimo(fimo_data, p_value, one_motif, genome, output_directory):

	#make the filename for the fimo output
	fimo_output_file = os.path.join(output_directory, "fimo_output_" + os.path.splitext(os.path.basename(one_motif))[0] + ".txt")
	fimo_output_unsorted = os.path.join(output_directory, "fimo_output_unsorted_" + os.path.splitext(os.path.basename(one_motif))[0] + ".txt")

	#check if user needs special options for the fimo
	if fimo_data != None: 
		fimo_data = fimo_data + " --thresh " + str(p_value) + " "
	else: 
		fimo_data = "--thresh " + str(p_value) + " " #add the passed p value to the fimo options

	#call fimo for this motif and save the output to a temporary file		
	send_to_fimo = "fimo --text --no-qvalue " + fimo_data + one_motif + " " + genome + " > " + fimo_output_unsorted

	logger.info('fimo proceed the data using this call ' + send_to_fimo)

	fimo_stdout = subprocess.getoutput(send_to_fimo)
	
	#the normal output from fimo starts with Using motif ..., so print it from the logger, otherwise print what else fimo says
	if fimo_stdout.startswith("Using") and re.split('\n', fimo_stdout)[1]:
		logger.info('info from fimo: ' + re.split('\n', fimo_stdout)[0].rstrip())
		logger.info('info from fimo: ' + re.split('\n', fimo_stdout)[1].rstrip())
	else: #there were some problems with fimo, so we want to see what they were
		logger.info('info from fimo: ' + fimo_stdout)
	
	
	if not os.path.isfile(fimo_output_unsorted):
		logger.info('the usage of fimo was crashed, there is no required output file, the exit is forced')
		sys.exit()

	if os.stat(fimo_output_unsorted).st_size == 0: #if the output of fimo is empty
		fimo_output_unsorted = fimo_output_unsorted.replace('unsorted_', '')
		return fimo_output_unsorted
	else:
		#if the file was converted from pfm, the second column contains the positions, so we want to sort using this column, and not the next one
		second_line = subprocess.getoutput("sed -n '2{p;q}' " + fimo_output_unsorted)
		if re.split('\t', second_line)[2].startswith("chr"): #the re.split[1] is a ' ', so take the [2]
			os.system("cat " + fimo_output_unsorted + " | sort -k 2 -V > " + fimo_output_file)
		else:
			#we are sorting after the third column, which looks like chr1:123-126, -V means it will sort the digitals and not the strings
			os.system("cat " + fimo_output_unsorted + " | sort -k 3 -V > " + fimo_output_file)

		#make sure the output of fimo exists
		if not os.path.isfile(fimo_output_file):
			logger.info('the sorting of the output file from the fimo was crashed, the exit is forced')
			sys.exit()
		else:
			return fimo_output_file

def call_moods(one_motif, genome, output_directory, p_value, moods_bg):

	# setting standard parameters for moods
	pseudocount = 0.0001

	if moods_bg == None:
		bg = MOODS.tools.flat_bg(4)
	else:
		bg = tuple(moods_bg)

	logger.info("moods will work with the p_value " + str(p_value) + " and the bg " + str(bg))

	motif_name = os.path.basename(one_motif)

	moods_output_unsorted_name = os.path.join(output_directory, "moods_output_unsorted_" + os.path.splitext(motif_name)[0] + ".txt")
	moods_output_file_unsorted = open(moods_output_unsorted_name, 'w')

	moods_output_name = os.path.join(output_directory, "moods_output_" + os.path.splitext(motif_name)[0] + ".txt")
	moods_output_file = open(moods_output_name, 'w')

	matrix_names = [os.path.basename(one_motif)]

	matrices = []
	matrices_rc = []

	valid, matrix = pfm_to_log_odds(one_motif, bg, pseudocount)

	key_for_bed_dict = ''

	if valid:

		logger.info("please be patient, moods is working on the data")

		matrices.append(matrix)
		matrices_rc.append(MOODS.tools.reverse_complement(matrix,4))
		matrices_all = matrices + matrices_rc
		thresholds = [MOODS.tools.threshold_from_p(m, bg, p_value, 4) for m in matrices_all]

		scanner = MOODS.scan.Scanner(7)
		scanner.set_motifs(matrices_all, bg, thresholds)

		with open(genome) as handle:

			seq_iterator = bio.SimpleFastaParser(handle)
			
			for header, seq in seq_iterator:

				header_splitted = re.split(r':', header)

				if len(header_splitted) == 1: #if there are no positions given
					header = header + ":0-" #set the first position as 0 and split it once more
					header_splitted = re.split(r':', header)
					logger.info("moods works with " + header)
				else: #the given genome file is a file with peaks, so use the header of the peak as a key to search in the bed dictionary for additional information later on
					key_for_bed_dict = header

				chromosom = header_splitted[0]
				positions = re.split(r'-', header_splitted[-1])

				results = scanner.scan(seq)

				fr = results[:len(matrix_names)] #forward strand
				rr = results[len(matrix_names):] #reverse strand

				results = [[(r.pos, r.score, '+', ()) for r in fr[i]] + 
					[(r.pos, r.score, '-', ()) for r in rr[i]] for i in range(len(matrix_names))] #use + and - to indicate strand

				for (matrix, matrix_name, result) in zip(matrices, matrix_names, results):

					motif_id = re.split(r'_', matrix_name)[-1] #find the id of the given morif
					motif_alternate_name = matrix_name.replace(motif_id, '')[:-1] #the alternate name of the motif is the name of the file without id and with cutted last character, that is _

					if len(matrix) == 4:
						l = len(matrix[0])
					if len(matrix) == 16:
						l = len(matrix[0] + 1)
					for r in sorted(result, key=lambda r: r[0]):
						strand = r[2]
						pos = r[0]
						hitseq = seq[pos:pos+l] #sequence
						#score = r[1]
						score = format(r[1], '.15f') #round to 15 digits after floating point, already type str

						if key_for_bed_dict != '':
							start = pos + 1
							end = pos + len(hitseq)
							chromosom = key_for_bed_dict #instead of only the name of chromosom write the key to search in the bed_file					
						else:
							start = int(positions[0]) + pos + 1
							end = start + len(hitseq) - 1
						
						#moods_output_file_unsorted.write('\t'.join([motif_id, motif_alternate_name, chromosom, str(start), str(end), strand, str(score)]) + '\n')
						moods_output_file_unsorted.write('\t'.join([motif_id, motif_alternate_name, chromosom, str(start), str(end), strand, score]) + '\n')

		#now sort the output of moods
		os.system("cat " + moods_output_unsorted_name + " | sort -k 1 -V > " + moods_output_name)

		moods_output_file_unsorted.close()
		moods_output_file.close()

		return moods_output_name

	else:
		logger.info("The input for moods was not validated by the MOODS.parsers.pfm. Please check if it has the right format (note that the MOODS accepts only the old version of .pfm files, that is one without the header containing the name and id of the motif)")
		sys.exit()

#help function for the moods call, convert pfm to log odds
def pfm_to_log_odds(filename, bg, pseudocount):
	if pfm(filename):
		mat = MOODS.parsers.pfm_to_log_odds(filename, bg, pseudocount)
		if len(mat) != 4: #if something went wrong, the empty list will be returned
			return False, mat
		else:
			return True, mat
	else:
		logger.info('please make sure the motif file has a .pfm format needed for moods')
		sys.exit()

#help function for the moods call, check if the file is in a pfm format using moods
def pfm(filename):
	mat = MOODS.parsers.pfm(filename)
	if len(mat) != 4:
		return False
	else:
		return True

# calculate the real positions of TFs, if needed, resolve the overlaps, and write to the output file
def write_output_file(input_file, bed_dictionary, resolve_overlaps):

	if os.path.basename(input_file).startswith("moods"):
		name_without_moods_or_fimo = input_file.replace('moods_output_', '')
		used_tool = "moods"
	else:
		name_without_moods_or_fimo = input_file.replace('fimo_output_', '')
		used_tool = "fimo"

	output_file_name = os.path.splitext(name_without_moods_or_fimo)[0] + ".bed"

	logger.info('writing the output file ' + output_file_name)
	output_file = open(output_file_name, 'w')
	
	#calculate the real positions of TFs and write the information in the output file
	with open(input_file) as read_file:
		
		overlap = []
		printed_line = []
		last_line = []
		key_for_bed_dict = ''

		for line in read_file:
			if not line.startswith('#'):
				
				line_to_write = []
				
				line_array = re.split(r'\t', line.rstrip('\n'))
				
				chromosom_and_positions = re.split(r':', line_array[2])
				if len(chromosom_and_positions) == 1: #the whole genome was given, there is nothing to split
					chromosom = line_array[2]
					start = line_array[3]
					end = line_array[4]
				else:
					positions = re.split(r'-', chromosom_and_positions[-1])
					chromosom = chromosom_and_positions[0]
					start = str(int(positions[0]) + int(line_array[3]))
					end = str(int(positions[0]) + int(line_array[4]))	
					key_for_bed_dict = line_array[2] #use only if there is a bed_dictionary in input				

				#------- these are 5 needed columns to succesfully proceed the data

				name = os.path.splitext(os.path.basename(name_without_moods_or_fimo))[0]

				score = line_array[6]

				line_to_write.extend([chromosom, start, end, name, score])
				#------ here the additional information coule be added to the output file

				strand_inf = line_array[5]
				line_to_write.append(strand_inf)

				if used_tool == "fimo":
					p_value = line_array[7]
					line_to_write.append(p_value)

				#if the dictionary is not empty check for the information corresponding to these positions
				if bed_dictionary and key_for_bed_dict in bed_dictionary:
					line_to_write.append('\t'.join(bed_dictionary[key_for_bed_dict]))

				line_to_write.insert(0, "write") #insert a flag to know if the line should be written or not

				last_line = line_to_write #save the line in case it is the last line and if due to the check_overlap it could not be printed

				if resolve_overlaps:
					overlap, line_to_write, printed_line = check_overlap(line_to_write, overlap, printed_line, output_file)

				write_line_not_overlap(output_file, line_to_write)

		if not last_line[0].startswith('write'): #it could be that the write flag was deleted from the last_line so append it back 
			overlap.insert(0, "write")

		#if there is already any printed line, check if the saved last line was already printed. Otherwise print it
		if resolve_overlaps:
			if printed_line:
				if last_line[1] != printed_line[1] or last_line[2] != printed_line[2]:
					write_line_not_overlap(output_file, last_line)

	output_file.close()

def write_line_not_overlap(output_file, line_to_write):
	if line_to_write: #if line_to_write is not empty
		if line_to_write[0].startswith('write'): #if the line does not contain an overlap, it contains the flag "write" at the first position
			line_to_write.pop(0) #delete the flag
			output_file.write('\t'.join(line_to_write)  + '\n')

def check_overlap(line_to_write, overlap, printed_line, output_file):

	is_overlap = None

	if not overlap: #if the overlap list is empty
		is_overlap = False
		
	else: #if the overlap list is not empty

		if not overlap[0].startswith('write'): #it could be that the write flag was deleted from the overlap so append it back to make sure the next if clauses run right
			overlap.insert(0, "write")

		if line_to_write[1] == overlap[1] and float(line_to_write[2]) < float(overlap[3]): #the current line could overlap the previous because the start of the current line is smaller than the end of the previous one and they are both on the one chromosom
			#if p value in the current line is smaller than p value in the previous line (or score is bigger), save the current line as possible overlap for future 
			if float(line_to_write[5]) > float(overlap[5]):
				is_overlap = False
			else:
				#if the p value in the current line is greater than p value in the previous line or are these both p values the same, the current line will not be printed, but also it will not be saved
				line_to_write.pop(0)
				is_overlap = None
		else: #it is an other chromosom or the start at the current line is greater or the same as the end of the previous one
			if printed_line != overlap:
				is_overlap = True
			else:
				is_overlap = False

	if is_overlap == False:
		overlap = line_to_write #save the current line
		line_to_write.pop(0) #do not print anything due to deleting the flag ("write")
	elif is_overlap == True:
		if not overlap[0].startswith('write'):
			overlap.insert(0, "write")
		printed_line = overlap #the previous line is saved as temporary line to print it later on
		overlap = line_to_write #save the current line
		line_to_write = printed_line #print the temporary saved line

	return overlap, line_to_write, printed_line

def remove_file(file):
	if os.path.isfile(file):
		os.remove(file)

def clean_directory(cleans, output_directory, motif, tool_output_file):

	fimo_output_unsorted = os.path.join(tool_output_file.replace("fimo_output", "fimo_output_unsorted"))
	moods_output_unsorted = os.path.join(tool_output_file.replace("moods_output", "moods_output_unsorted"))
	
	for clean in cleans:
		if clean == 'all':
			remove_file(motif)
			remove_file(tool_output_file)
			remove_file(fimo_output_unsorted)
			remove_file(moods_output_unsorted)
		elif clean == 'fimo_output':
			if os.path.basename(tool_output_file).startswith("fimo"):
				remove_file(tool_output_file)
				remove_file(fimo_output_unsorted)
		elif clean == 'cut_motifs':
			remove_file(motif)
		elif clean == 'moods_output':
			if os.path.basename(tool_output_file).startswith("moods"):
				remove_file(tool_output_file)
				remove_file(fimo_output_unsorted)

def tool_make_output(usage, motif, genome, output_directory, cleans, p_value, bed_dictionary, fimo_data, resolve_overlaps, moods_bg):
	try:
		if usage == "fimo":
			tool_output_file = call_fimo(fimo_data, p_value, motif, genome, output_directory)
		elif usage == "moods":
			tool_output_file = call_moods(motif, genome, output_directory, p_value, moods_bg)

		output = write_output_file(tool_output_file, bed_dictionary, resolve_overlaps)
	finally:
		clean_directory(cleans, output_directory, motif, tool_output_file)

def multiprocess(motifs, genome, output_directory, cleans, fimo_data, p_value, bed_dictionary, cpu_count, resolve_overlaps, usage, moods_bg):

	if cleans == None:
		cleans = ['all']

	pool = multiprocessing.Pool(cpu_count) #by default is cpu_count 2

	length = len(motifs) #the number of the motifs to find the percentage of the job that was done
	step = 100/length #the percentage that should be added after the job with each single motif is done

	tasks = [] #here the jobs done by processes are saved

	for motif in motifs:
		tasks.append(pool.apply_async(tool_make_output, args = (usage, motif, genome, output_directory, cleans, p_value, bed_dictionary, fimo_data, resolve_overlaps, moods_bg, )))

	tasks_done = sum([task.ready() for task in tasks]) #the number of the processes that ended their job

	#check the number of the processes that are ready till the number of them reaches the number of processes started in the pool
	while tasks_done < len(tasks):
		#if the number of ready processes has changed, save the new number and print the percentage of the job done
		if sum([task.ready() for task in tasks]) != tasks_done:
			tasks_done = sum([task.ready() for task in tasks])
			sys.stdout.write("%-100s| %d%% \r" % (''*tasks_done, step*tasks_done))
			sys.stdout.flush()
			sys.stdout.write("\n")
		#update the number of ready processes each 0.05 seconds
		time.sleep(0.05)

	pool.close()
	pool.join() #make sure all the processes are done and exit
		
	#the processes should not delete the merged genome file
	#so make sure if this file is needed, otherwise delete it
	for clean in cleans:
		if clean == 'all' or clean == 'merge_output':
			for filename in os.listdir(output_directory):
				if filename == "output_merge.fa":
					remove_file(genome)

		if clean != 'nothing':
			logger.info('the directory ' + output_directory + ' was cleaned, only the required files were left')

def make_bed_dictionary(bed_file):

	bed_dictionary = {}

	with open(bed_file) as read_bed_file:
		for bed_line in read_bed_file:
			bed_line_array = re.split(r'\t', bed_line.rstrip('\n'))
			if bed_line_array[1].isdigit() and bed_line_array[2].isdigit() and int(bed_line_array[1]) <= int(bed_line_array[2]): #in the real bedfile the second column is a start position, and the third column is an end position, so we are checking if these are integers and if the start position is smaller than the end one
				key = bed_line_array[0] + ":" + bed_line_array[1] + "-" + bed_line_array[2]
				value = []
				for i in range(3, len(bed_line_array)):
					value.append(bed_line_array[i]) 

				bed_dictionary[key] = value
			else: #this is not a bed file, force exit
				logger.info('please make sure the input bed file has a right format, the problem occured on the line ' + bed_line)
				sys.exit()

	read_bed_file.close()

	return bed_dictionary

def is_fasta(check_fasta):
	if not os.path.isfile(check_fasta):
		logger.info('there is no file with genome, the exit is forced')
		sys.exit()
	else:
		# modified code from https://stackoverflow.com/questions/44293407/how-can-i-check-whether-a-given-file-is-fasta

	    with open(check_fasta, "r") as handle:
	        fasta = SeqIO.parse(handle, "fasta")
	        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file

def check_existing_input_files(args):
	if not args.motifs or not args.genome:
		logger.info('there is no satisfied input, please enter --help or -h to learn how to use this tool')
		sys.exit()
	elif not is_fasta(args.genome):
		logger.info('please make sure the input genome file has a fasta format')
		sys.exit()
	#check if the file with motifs exists
	elif not os.path.isfile(args.motifs):
		logger.info('there is no file with motifs, the exit is forced')
		sys.exit()
	#if the bedfile was given as input, check if this file exists
	elif args.bed_file:
		if not os.path.isfile(args.bed_file):
			logger.info('there is no such bed file ' + args.bed_file + ', the exit is forced')
			sys.exit()

def check_fimo_version():
	fimo_version = subprocess.getoutput("fimo --version") #works only on python 3.4
	fimo_version = int(fimo_version.replace('.', '')) #reaplace the . in the version to be able to compare it as int
	if fimo_version < 4120:
		logger.info('please make sure you are using fimo version 4.12.0 or the newer one')
		sys.exit()

def main():


	args = parse_args()
	
	if args.use == "fimo":
		check_fimo_version()

	#if user do not want to see the information about the status of jobs, remove the handler, that writes to the terminal
	if args.hide_info:
		#logger.disabled = True
		logger.removeHandler(ch)

	check_existing_input_files(args)

	#check if there is an existing directory that user gave as input, otherwise create this directory from the path provided from the user
	check_directory(args.output_directory)

	splitted_motifs = split_motifs(args.motifs, args.output_directory, args.use)
	
	#check if there is a .bed file to merge the genome file with. If so, merge them
	if args.bed_file:
		bed_dictionary = make_bed_dictionary(args.bed_file)
		args.genome = merge(args.genome, args.bed_file, args.output_directory)
	else:
		bed_dictionary = {}
		
	#if the usage is moods, call moods, otherwise call fimo
	multiprocess(splitted_motifs, args.genome, args.output_directory, args.cleans, args.fimo, args.p_value, bed_dictionary, args.cores, args.resolve_overlaps, args.use, args.moods_bg)
	
	for handler in logger.handlers:
		handler.close()
		logger.removeFilter(handler)
	
if __name__ == "__main__":
    main()