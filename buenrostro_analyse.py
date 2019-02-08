import argparse
import sys
import os
import re
import time
import logging
import subprocess
import timeit
import numpy as np
#import matplotlib.pylab as plt
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import stats

logger = logging.getLogger('buenrosto_analyse')
logger.setLevel(logging.INFO)

formatter = logging.Formatter("%(asctime)s : %(message)s", "%Y-%m-%d %H:%M")

fh = logging.FileHandler('buenrosto_log.txt')
fh.setLevel(logging.INFO)
fh.setFormatter(formatter)
logger.addHandler(fh)

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)

#catch all the information about input and output files as well as information on the used tool (fimo or moods)
def parse_args():
	parser = argparse.ArgumentParser(prog='analyse', description='This script analyses the mytool.py while comparing the output from the mytool.py with the given chipseq data', epilog='That is what you need to make this script work for you. Enjoy it')

	parser.add_argument('-m', '--motifs', help='a list of motifs that are looked for in the genome')
	parser.add_argument('-g', '--genome', help='a genome file')
	parser.add_argument('-c', '--chipseq', help='a directory with chipseq files to work with')
	parser.add_argument('-o', '--output_directory', default='output', const='output', nargs='?', help='please set the directory for the output files, if there is no path provided, the new directory named output will be created in the current directory. please note that you are required to create a new directory and do not use a directory which consists some files already!')
	parser.add_argument('-b', '--bed_file', nargs='?', help='a .bed file to be merged with the whole genome file')
	#parser.add_argument('--clean', nargs='*', choices=['nothing', 'all', 'cut_motifs', 'fimo_output', 'merge_output', 'moods_output'], dest='cleans', help='choose the files you want to delete from the output directory, the default is deleting all the temporary files from the directory')
	parser.add_argument('--encode', help='a file to encode jaspar names')

	args = parser.parse_args()
	return args

def check_directory(directory):
	if not os.path.exists(directory):
		logger.info("the new directory " + directory + " was created")
		os.makedirs(directory)

def make_plot(x_array, y_array, motifs_names, xlabel_name, ylabel_name, output_directory, figure_name, color_name):
		
	plt.xlim(0,1)
	plt.ylim(0,1)

	fig, ax = plt.subplots()
	plt.plot([0,1], '--', color = "black", linewidth = 0.5)
	plt.plot(x_array, y_array, 'o', color = color_name)

	plt.xlabel(xlabel_name)
	plt.ylabel(ylabel_name)

	for i,j,name in zip(x_array, y_array, motifs_names):
		ax.annotate('%s' %name, xy=(i,j), xytext=(1,0), textcoords='offset points', size='xx-small')

	fig.savefig(os.path.join(output_directory, figure_name))
	

	#------------------------ the old code

	#df = pd.DataFrame({"x": x_array, "y": y_array, 'group': motifs_names})

	#p = sns.regplot(data = df, x = "x", y = "y", fit_reg = True, marker = "o", color = color_name, scatter_kws = {'s':20}, line_kws={'linewidth': 0.5})
	#add annotations one by one with a loop
	#for line in range(0, df.shape[0]):
	#	p.text(df.x[line], df.y[line], df.group[line], horizontalalignment = 'left', size = 'xx-small', color = 'black', weight = 'light')

	#use this code to save the plot authomaticaly	
	#fig = p.get_figure()
	#fig.savefig(os.path.join(output_directory, figure_name))

	plt.show()

def make_mean_plot(arrays_names, means_array):

	plt.ylim(0,1)
	x = []
	for i in range(len(arrays_names)):
		x.append(i+1)

	plt.plot(x, means_array, 'g')
	#put the names of arrays on the x axis
	plt.xticks(x, arrays_names, rotation="vertical")
	#pad margins so that markers don't get clipped by the axes
	plt.margins(0.2)
	#tweak spacing to prevent clipping of tick-labels
	plt.subplots_adjust(bottom=0.50)
	plt.grid()
	plt.show()
			
def find_mean(array):
	return sum(array)/len(array)

def read_analyse_and_make_plots(analyse_output_name):
	with open(analyse_output_name) as file:
		lines = 0

		motifs_names = []

		fimo_bg_hg19 = []

		fimo_without_overlaps = []
		fimo_without_overlaps_1e_2 = []
		fimo_without_overlaps_1e_3 = []
		fimo_without_overlaps_1e_5 = []
		fimo_without_overlaps_1e_6 = []

		fimo_without_overlaps_bg = []
		fimo_without_overlaps_bg_1e_2 = []
		fimo_without_overlaps_bg_1e_3 = []
		fimo_without_overlaps_bg_1e_5 = []
		fimo_without_overlaps_bg_1e_6 = []

		fimo_with_overlaps = []
		
		moods_without_overlaps = []
		moods_without_overlaps_1e_2 = []
		moods_without_overlaps_1e_3 = []
		moods_without_overlaps_1e_5 = []
		moods_without_overlaps_1e_6 = []

		moods_without_overlaps_bg = []
		moods_without_overlaps_bg_1e_2 = []
		moods_without_overlaps_bg_1e_3 = []

		moods_without_overlaps_new_bg = []
		moods_without_overlaps_new_bg_1e_2 = []
		moods_without_overlaps_new_bg_1e_3 = []
		moods_without_overlaps_new_bg_1e_5 = []
		moods_without_overlaps_new_bg_1e_6 = []

		moods_with_overlaps = []
		

		for line in file:
			
			line_array = re.split('\t', line.rstrip())
			if lines >= 1: 

				motifs_names.append(line_array[0])

				fimo_bg_hg19.append(float(line_array[1]))

				fimo_without_overlaps.append(float(line_array[2]))
				fimo_without_overlaps_1e_2.append(float(line_array[3]))
				fimo_without_overlaps_1e_3.append(float(line_array[4]))
				fimo_without_overlaps_1e_5.append(float(line_array[5]))
				fimo_without_overlaps_1e_6.append(float(line_array[6]))

				fimo_without_overlaps_bg.append(float(line_array[7]))
				fimo_without_overlaps_bg_1e_2.append(float(line_array[8]))
				fimo_without_overlaps_bg_1e_3.append(float(line_array[9]))
				fimo_without_overlaps_bg_1e_5.append(float(line_array[10]))
				fimo_without_overlaps_bg_1e_6.append(float(line_array[11]))

				fimo_with_overlaps.append(float(line_array[12]))

				moods_without_overlaps.append(float(line_array[13]))
				moods_without_overlaps_1e_2.append(float(line_array[14]))
				moods_without_overlaps_1e_3.append(float(line_array[15]))
				moods_without_overlaps_1e_5.append(float(line_array[16]))
				moods_without_overlaps_1e_6.append(float(line_array[17]))

				moods_without_overlaps_bg.append(float(line_array[18]))
				moods_without_overlaps_bg_1e_2.append(float(line_array[19]))
				moods_without_overlaps_bg_1e_3.append(float(line_array[20]))

				moods_without_overlaps_new_bg.append(float(line_array[21]))
				moods_without_overlaps_new_bg_1e_2.append(float(line_array[22]))
				moods_without_overlaps_new_bg_1e_3.append(float(line_array[23]))
				moods_without_overlaps_new_bg_1e_5.append(float(line_array[24]))
				moods_without_overlaps_new_bg_1e_6.append(float(line_array[25]))

				moods_with_overlaps.append(float(line_array[26]))

			lines = lines + 1

	file.close()

	means_array = []
	arrays_names = []


	arrays_names.append("fimo_bg_hg19")
	means_array.append(find_mean(fimo_bg_hg19))

	arrays_names.append("fimo_without_overlaps_1e-2")
	means_array.append(find_mean(fimo_without_overlaps_1e_2))

	arrays_names.append("fimo_without_overlaps_1e-3")
	means_array.append(find_mean(fimo_without_overlaps_1e_3))

	arrays_names.append("fimo_without_overlaps_1e-4")
	means_array.append(find_mean(fimo_without_overlaps))

	arrays_names.append("fimo_without_overlaps_1e-5")
	means_array.append(find_mean(fimo_without_overlaps_1e_5))

	arrays_names.append("fimo_without_overlaps_1e-6")
	means_array.append(find_mean(fimo_without_overlaps_1e_6))

	arrays_names.append("fimo_without_overlaps_bg_1e-2")
	means_array.append(find_mean(fimo_without_overlaps_bg_1e_2))

	arrays_names.append("fimo_without_overlaps_bg_1e-3")
	means_array.append(find_mean(fimo_without_overlaps_bg_1e_3))

	arrays_names.append("fimo_without_overlaps_bg_1e-4")
	means_array.append(find_mean(fimo_without_overlaps_bg))

	arrays_names.append("fimo_without_overlaps_bg_1e-5")
	means_array.append(find_mean(fimo_without_overlaps_bg_1e_5))

	arrays_names.append("fimo_without_overlaps_bg_1e-6") 
	means_array.append(find_mean(fimo_without_overlaps_bg_1e_6))

	arrays_names.append("fimo_with_overlaps")
	means_array.append(find_mean(fimo_with_overlaps))



	
	arrays_names.append("moods_without_overlaps_1e-2")
	means_array.append(find_mean(moods_without_overlaps_1e_2))

	arrays_names.append("moods_without_overlaps_1e-3")
	means_array.append(find_mean(moods_without_overlaps_1e_3))

	arrays_names.append("moods_without_overlaps_1e-4")
	means_array.append(find_mean(moods_without_overlaps))

	arrays_names.append("moods_without_overlaps_1e-5")
	means_array.append(find_mean(moods_without_overlaps_1e_5))

	arrays_names.append("moods_without_overlaps_1e-6")
	means_array.append(find_mean(moods_without_overlaps_1e_6))


	"""
	arrays_names.append("moods_without_overlaps_bg_1e-2")
	means_array.append(find_mean(moods_without_overlaps_1e_2))
	arrays_names.append("moods_without_overlaps_bg_1e-3")
	means_array.append(find_mean(moods_without_overlaps_1e_3))
	arrays_names.append("moods_without_overlaps_bg_1e-4")
	means_array.append(find_mean(moods_without_overlaps_bg))
	"""

	arrays_names.append("moods_without_overlaps_bg_1e-2")
	means_array.append(find_mean(moods_without_overlaps_new_bg_1e_2))

	arrays_names.append("moods_without_overlaps_bg_1e-3")
	means_array.append(find_mean(moods_without_overlaps_new_bg_1e_3))

	arrays_names.append("moods_without_overlaps_bg_1e-4")
	means_array.append(find_mean(moods_without_overlaps_new_bg))

	arrays_names.append("moods_without_overlaps_bg_1e-5")
	means_array.append(find_mean(moods_without_overlaps_new_bg_1e_5))

	arrays_names.append("moods_without_overlaps_bg_1e-6")
	means_array.append(find_mean(moods_without_overlaps_new_bg_1e_6))

	arrays_names.append("moods_with_overlaps")
	means_array.append(find_mean(moods_with_overlaps))

	
	conditions_arrays = [[] for condition in range(6)]
	#conditions_arrays will consist of: fimo_with_overlaps, fimo_without_overlaps, fimo_without_overlaps_bg
	#									moods_with_overlaps, moods_without_overlaps, moods_without_overlaps_bg
	#fill the conditions_arrays with pvalues: 1e-2, 1e-3, 1e-4, 1e-5, 1e-6
	#use np.nan to mask the missing values
	cond_fimo_with_overlaps = [np.nan, np.nan, find_mean(fimo_with_overlaps), np.nan, np.nan]
	cond_fimo_without_overlaps = [find_mean(fimo_without_overlaps_1e_2), find_mean(fimo_without_overlaps_1e_3), find_mean(fimo_without_overlaps), find_mean(fimo_without_overlaps_1e_5), find_mean(fimo_without_overlaps_1e_6)]
	cond_fimo_without_overlaps_bg = [find_mean(fimo_without_overlaps_bg_1e_2), find_mean(fimo_without_overlaps_bg_1e_3), find_mean(fimo_without_overlaps_bg), find_mean(fimo_without_overlaps_1e_5), find_mean(fimo_without_overlaps_bg_1e_6)]

	cond_moods_with_overlaps = [np.nan, np.nan, find_mean(moods_with_overlaps), np.nan, np.nan]
	cond_moods_without_overlaps = [find_mean(moods_without_overlaps_1e_2), find_mean(moods_without_overlaps_1e_3), find_mean(moods_without_overlaps), find_mean(moods_without_overlaps_1e_5), find_mean(moods_without_overlaps_1e_6)]
	cond_moods_without_overlaps_bg = [find_mean(moods_without_overlaps_new_bg_1e_2), find_mean(moods_without_overlaps_new_bg_1e_3), find_mean(moods_without_overlaps_new_bg), find_mean(moods_without_overlaps_new_bg_1e_5), find_mean(moods_without_overlaps_new_bg_1e_6)]

	#make numpy matrix
	#data = np.stack((cond_fimo_with_overlaps, cond_moods_without_overlaps, cond_fimo_without_overlaps_bg, cond_moods_with_overlaps, cond_moods_without_overlaps, cond_moods_without_overlaps_bg))
	data = np.stack((cond_fimo_without_overlaps, cond_fimo_without_overlaps_bg, cond_moods_without_overlaps, cond_moods_without_overlaps_bg))

	column_labels = ('1e-2', '1e-3', '1e-4', '1e-5', '1e-6')
	#row_labels = ('fimo_with_overlaps', 'fimo_without_overlaps', 'fimo_without_overlaps_bg', 'moods_with_overlaps', 'moods_without_overlaps', 'moods_without_overlaps_bg')
	row_labels = ('fimo_without_overlaps', 'fimo_without_overlaps_bg', 'moods_without_overlaps', 'cond_moods_without_overlaps_bg')
	fig, ax = plt.subplots()

	#heatmap = ax.pcolor(data, cmap = plt.cm.cool,  vmin=np.nanmin(data), vmax=np.nanmax(data)) #mask the missing data using the white colour while plotting
	#heatmap.cmap.set_under('blue')

	heatmap = ax.pcolor(data, cmap = plt.cm.cool)

	#put the major ticks at the middle of each cell
	ax.set_xticks(np.arange(len(data) + 1) + 0.5)
	ax.set_yticks(np.arange(len(data)) + 0.5)

	#want a more natural, table-like display
	ax.invert_yaxis()
	ax.xaxis.tick_top()

	ax.set_xticklabels(column_labels, minor=False)
	ax.set_yticklabels(row_labels, minor=False)

	plt.colorbar(heatmap)

	plt.show()


	#---------------make plots

	make_plot(fimo_with_overlaps, moods_with_overlaps, motifs_names, "fimo_with_overlaps", "moods_with_overlaps", "buenrostro_analyse/sensitivity", "fimo_vs_moods_with_overlaps.png", "deepskyblue")

	make_plot(fimo_without_overlaps, moods_without_overlaps, motifs_names, "fimo_without_overlaps", "moods_without_overlaps", "buenrostro_analyse/sensitivity", "fimo_vs_moods_without_overlaps.png", "dodgerblue")

	make_plot(fimo_without_overlaps, fimo_without_overlaps_1e_3, motifs_names, "fimo_without_overlaps", "fimo_without_overlaps_1e-3", "buenrostro_analyse/sensitivity", "fimo_vs_fimo_1e-3_without_overlaps.png", "gold")

	make_plot(fimo_without_overlaps, fimo_without_overlaps_1e_5, motifs_names, "fimo_without_overlaps", "fimo_without_overlaps_1e-5", "buenrostro_analyse/sensitivity", "fimo_vs_fimo_1e-5_without_overlaps.png", "hotpink")

	make_plot(moods_without_overlaps, moods_without_overlaps_1e_3, motifs_names, "moods_without_overlaps", "moods_without_overlaps_1e-3", "buenrostro_analyse/sensitivity", "moods_vs_moods_1e-3_without_overlaps.png", "orange")
	
	make_plot(moods_without_overlaps, moods_without_overlaps_1e_5, motifs_names, "moods_without_overlaps", "moods_without_overlaps_1e-5", "buenrostro_analyse/sensitivity", "moods_vs_moods_1e-5_without_overlaps.png", "orchid")

	make_plot(fimo_without_overlaps_1e_3, moods_without_overlaps_1e_3, motifs_names, "fimo_without_overlaps_1e-3", "moods_without_overlaps_1e-3", "buenrostro_analyse/sensitivity", "fimo_1e-3_vs_moods_1e-3_without_overlaps.png", "mediumspringgreen")

	make_plot(fimo_without_overlaps_1e_5, moods_without_overlaps_1e_5, motifs_names, "fimo_without_overlaps_1e-5", "moods_without_overlaps_1e-5", "buenrostro_analyse/sensitivity", "fimo_1e-5_vs_moods_1e-5_without_overlaps.png", "aquamarine")

	make_plot(fimo_bg_hg19, fimo_with_overlaps, motifs_names, "fimo_bg_hg19", "fimo_with_overlaps", "buenrostro_analyse/sensitivity", "fimo_bg_hg19_vs_fimo_with_overlaps.png", "coral")

	make_plot(fimo_without_overlaps_bg, fimo_without_overlaps, motifs_names, "fimo_without_overlaps_bg", "fimo_without_overlaps", "buenrostro_analyse/sensitivity", "fimo_without_overlaps_bg_vs_fimo_without_overlaps", "orangered")

	make_plot(fimo_without_overlaps_1e_3, fimo_without_overlaps_1e_2, motifs_names, "fimo_without_overlaps_1e_3", "fimo_without_overlaps_1e_2", "buenrostro_analyse/sensitivity", "fimo_without_overlaps_1e_3_vs_fimo_without_overlaps_1e_2", "indianred")

	make_plot(fimo_without_overlaps_bg_1e_3, fimo_without_overlaps_1e_3, motifs_names, "fimo_without_overlaps_bg_1e_3", "fimo_without_overlaps_1e_3", "buenrostro_analyse/sensitivity", "fimo_without_overlaps_bg_1e_3_vs_fimo_without_overlaps_1e_3", "chocolate")

	make_plot(fimo_without_overlaps_bg_1e_2, fimo_without_overlaps_1e_2, motifs_names, "fimo_without_overlaps_bg_1e_2", "fimo_without_overlaps_1e_2", "buenrostro_analyse/sensitivity", "fimo_without_overlaps_bg_1e_2_vs_fimo_without_overlaps_1e_2", "tomato")

	make_plot(moods_without_overlaps, moods_without_overlaps_new_bg, motifs_names, "moods_without_overlaps", "moods_without_overlaps_new_bg", "buenrostro_analyse/sensitivity", "moods_without_overlaps_vs_moods_without_overlaps_new_bg", "springgreen")

	make_bar_plot_fimo(tuple(motifs_names), tuple(fimo_without_overlaps), tuple(fimo_without_overlaps_1e_3), tuple(fimo_without_overlaps_1e_5))
	make_bar_plot_moods(tuple(motifs_names), tuple(moods_without_overlaps), tuple(moods_without_overlaps_1e_3), tuple(moods_without_overlaps_1e_5))
	make_mean_plot(arrays_names, means_array)

def make_bar_plot_fimo(motifs_names, fimo_without_overlaps, fimo_without_overlaps_1e_3, fimo_without_overlaps_1e_5):
	
	N = len(motifs_names)

	ind = np.arange(N)    # the x locations for the groups
	width = 0.5      # the width of the bars: can also be len(x) sequence


	p1 = plt.bar(ind, fimo_without_overlaps_1e_3, width, color = 'paleturquoise', alpha = 0.5)
	p2 = plt.bar(ind, fimo_without_overlaps, width, color = 'royalblue', alpha = 0.5)
	p3 = plt.bar(ind, fimo_without_overlaps_1e_5, width, color = 'darkblue', alpha = 0.5)

	plt.ylabel('Scores')
	#plt.title('Scores by group and gender')
	plt.xticks(ind, motifs_names, rotation = "vertical", size='xx-small')
	plt.margins(0.2)
	plt.subplots_adjust(bottom=0.30)
	#plt.yticks(np.arange(0, 81, 10))
	plt.legend((p1[0], p2[0], p3[0]), ('fimo_without_overlaps_1e_3', 'fimo_without_overlaps', 'fimo_without_overlaps_1e_5'))

	plt.show()

def make_bar_plot_moods(motifs_names, moods_without_overlaps, moods_without_overlaps_1e_3, moods_without_overlaps_1e_5):

	N = len(motifs_names)

	ind = np.arange(N)    # the x locations for the groups
	width = 0.5      # the width of the bars: can also be len(x) sequence


	p1 = plt.bar(ind, moods_without_overlaps_1e_3, width, color = 'pink', alpha = 0.5)
	p2 = plt.bar(ind, moods_without_overlaps, width, color = 'hotpink', alpha = 0.5)
	p3 = plt.bar(ind, moods_without_overlaps_1e_5, width, color = 'crimson', alpha = 0.5)

	plt.ylabel('Scores')
	#plt.title('Scores by group and gender')
	plt.xticks(ind, motifs_names, rotation = "vertical", size='xx-small')
	plt.margins(0.2)
	plt.subplots_adjust(bottom=0.30)
	#plt.yticks(np.arange(0, 81, 10))
	plt.legend((p1[0], p2[0], p3[0]), ('moods_without_overlaps_1e_3', 'moods_without_overlaps', 'moods_without_overlaps_1e_5'))

	plt.show()

def main():

	# python buenrostro_analyse.py -m ../../../../PaperInPrep/TOBIAS/buenrostro_analysis/chipseq_GM12878_hg19/buenrostro_motifs.meme -g ../../../analysis_my_tool/chipseq/hg19.fasta -b buenrostro_union.bed -o buenrostro_analyse -c ../../../../PaperInPrep/TOBIAS/buenrostro_analysis/chipseq_GM12878_hg19/peaks --encode new_encode_jaspar_name.txt

	start_time = timeit.default_timer()
	
	args = parse_args()

	#check if the output_directory exists and if not, create one	
	check_directory(args.output_directory)
	
	#make the directories within the output_directory to save the results of mytool
	#output_dir_fimo_overlaps = os.path.join(args.output_directory, "fimo_with_overlaps")
	#output_dir_fimo_without_overlaps = os.path.join(args.output_directory, "fimo_without_overlaps")

	#output_dir_moods_overlaps = os.path.join(args.output_directory, "moods_with_overlaps")
	#output_dir_moods_without_overlaps = os.path.join(args.output_directory, "moods_without_overlaps")

	#output_dir_fimo_without_overlaps_1e_3 = os.path.join(args.output_directory, "fimo_without_overlaps_1e-3")
	#output_dir_moods_without_overlaps_1e_3 = os.path.join(args.output_directory, "moods_without_overlaps_1e-3")

	#output_dir_fimo_without_overlaps_1e_5 = os.path.join(args.output_directory, "fimo_without_overlaps_1e-5")
	#output_dir_moods_without_overlaps_1e_5 = os.path.join(args.output_directory, "moods_without_overlaps_1e-5")

	#output_dir_fimo_bg = os.path.join(args.output_directory, "fimo_bg_hg19")
	#output_dir_fimo_bg_without_overlaps = os.path.join(args.output_directory, "fimo_without_overlaps_bg")
	#output_dir_moods_bg_5 = os.path.join(args.output_directory, "moods_bg_5")
	
	#output_dir_fimo_1e_2 = os.path.join(args.output_directory, "fimo_without_overlaps_1e_2")
	#output_dir_moods_1e_2 = os.path.join(args.output_directory, "moods_without_overlaps_1e_2")

	#output_dir_fimo_bg_1e_3 = os.path.join(args.output_directory, "fimo_without_overlaps_bg_1e_3")
	#output_dir_moods_bg_1e_3 = os.path.join(args.output_directory, "moods_without_overlaps_bg_1e_3")

	#output_dir_fimo_bg_1e_2 = os.path.join(args.output_directory, "fimo_without_overlaps_bg_1e_2")
	#output_dir_moods_bg_1e_2 = os.path.join(args.output_directory, "moods_without_overlaps_bg_1e_2")

	#output_dir_fimo_bg_1e_5 = os.path.join(args.output_directory, "fimo_without_overlaps_bg_1e_5")

	#output_dir_fimo_1e_6 = os.path.join(args.output_directory, "fimo_without_overlaps_1e_6")
	#output_dir_fimo_bg_1e_6 = os.path.join(args.output_directory, "fimo_without_overlaps_bg_1e_6")

	#output_dir_moods_1e_6 = os.path.join(args.output_directory, "moods_without_overlaps_1e_6") 

	#output_dir_moods_new_bg = os.path.join(args.output_directory, "moods_without_overlaps_new_bg")
	#output_dir_moods_new_bg_1e_2 = os.path.join(args.output_directory, "moods_without_overlaps_new_bg_1e_2")
	#output_dir_moods_new_bg_1e_3 = os.path.join(args.output_directory, "moods_without_overlaps_new_bg_1e_3")
	#output_dir_moods_new_bg_1e_5 = os.path.join(args.output_directory, "moods_without_overlaps_new_bg_1e_5")
	#output_dir_moods_new_bg_1e_6 = os.path.join(args.output_directory, "moods_without_overlaps_new_bg_1e_6")

	#-----------------------------------------------------------------------------------------------------

	#call mytool to work with the data
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_fimo_overlaps + " -m " + args.motifs + " --use fimo")
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_fimo_without_overlaps + " -m " + args.motifs + " --use fimo --resolve_overlaps")
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_moods_overlaps + " -m " + args.motifs + " --use moods")
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_moods_without_overlaps + " -m " + args.motifs + " --use moods --resolve_overlaps")

	#default pvalue 1e-4 or 0.0001
	#test pvalue 1e-3 or 0.001
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_fimo_without_overlaps_1e_3 + " -m " + args.motifs + " --use fimo --resolve_overlaps -p 0.001")
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_moods_without_overlaps_1e_3 + " -m " + args.motifs + " --use moods --resolve_overlaps -p 0.001")

	#test pvalue 1e-5 or 0.00001
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_fimo_without_overlaps_1e_5 + " -m " + args.motifs + " --use fimo --resolve_overlaps -p 0.00001")
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_moods_without_overlaps_1e_5 + " -m " + args.motifs + " --use moods --resolve_overlaps -p 0.00001")

	#test another bg from the file found with help of fasta-get-markov
	#this bg has for g and c values 0,2045, so we will try to take the 5 for moods background as i think it will result the 0,2 frequency for g and c
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_fimo_bg_without_overlaps + " -m " + args.motifs + " --use fimo --fimo \"--bgfile markov_bg_hg19.txt\" --resolve_overlaps")
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_moods_bg_5 + " -m " + args.motifs + " --use moods --resolve_overlaps --moods_bg 5")
	
	#test pvalue 1e-2
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_fimo_1e_2 + " -m " + args.motifs + " --use fimo --resolve_overlaps -p 0.01")
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_moods_1e_2 + " -m " + args.motifs + " --use moods --resolve_overlaps -p 0.01")

	#test bg and pvalue 1e-3
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_fimo_bg_1e_3 + " -m " + args.motifs + " --use fimo --fimo \"--bgfile markov_bg_hg19.txt\" --resolve_overlaps -p 0.001")
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_moods_bg_1e_3 + " -m " + args.motifs + " --use moods --resolve_overlaps --moods_bg 5 -p 0.001")

	#test bg and pvalue 1e-2
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_fimo_bg_1e_2 + " -m " + args.motifs + " --use fimo --fimo \"--bgfile markov_bg_hg19.txt\" --resolve_overlaps -p 0.01")
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_moods_bg_1e_2 + " -m " + args.motifs + " --use moods --resolve_overlaps --moods_bg 5 -p 0.01")

	#test pvalue 1e-6 with and without bg
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_fimo_1e_6 + " -m " + args.motifs + " --use fimo --resolve_overlaps -p 0.000001")
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_fimo_bg_1e_6 + " -m " + args.motifs + " --use fimo --fimo \"--bgfile markov_bg_hg19.txt\" --resolve_overlaps -p 0.000001")

	#two missing fimo datas with bg
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_fimo_bg + " -m " + args.motifs + " --use fimo --fimo \"--bgfile markov_bg_hg19.txt\" --resolve_overlaps")
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_fimo_bg_1e_5 + " -m " + args.motifs + " --use fimo --fimo \"--bgfile markov_bg_hg19.txt\" --resolve_overlaps -p 0.00001")
	
	#i have hard changed the bg in moods and set it as a list, so now i dont need to set bg here
	#test moods with the bg setted directly in the mytool.py
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_moods_new_bg + " -m " + args.motifs + " --use moods --resolve_overlaps")
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_moods_new_bg_1e_2 + " -m " + args.motifs + " --use moods --resolve_overlaps -p 0.01")
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_moods_new_bg_1e_3 + " -m " + args.motifs + " --use moods --resolve_overlaps -p 0.001")
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_moods_new_bg_1e_5 + " -m " + args.motifs + " --use moods --resolve_overlaps -p 0.00001")
	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_moods_new_bg_1e_6 + " -m " + args.motifs + " --use moods --resolve_overlaps -p 0.000001")

	#os.system("python mytool.py -g " + args.genome + " --bed " + args.bed_file + " -o " + output_dir_moods_1e_6 + " -m " + args.motifs + " --use moods --resolve_overlaps -p 0.000001")
	#-------------------------------------------------------------------------------------------------------

	"""
	#for the given chipseq files find the regions within the given peaks
	if os.path.isfile(args.encode):
		with open(args.encode) as encode:
			for line in encode:
				motif_id = re.split('\t', line)[1]
				alternate_name = re.split('\t', line)[2].rstrip()
				name_to_search = os.path.join(args.chipseq, alternate_name + ".bed")
				
				if os.path.isfile(name_to_search):
					chipseq_sorted = os.path.join(args.output_directory, "sorted_" + alternate_name + "_" + motif_id + ".bed")
					os.system("cat " + name_to_search + " | sort -k1 -V > " + chipseq_sorted)
					#find the positions of chipseqfile within the .bed file with peaks
					chipseq_within_peaks = chipseq_sorted.replace("sorted", "within_peaks")
					os.system("bedtools intersect -a " + chipseq_sorted + " -b " + args.bed_file + " -u > " + chipseq_within_peaks) #-u report original entry from chipseq_sorted if at leas one overlap from bed_file was found
					if os.path.isfile(chipseq_sorted):
						os.remove(chipseq_sorted)
	"""
	#---------------make the analyse_output each column 
	"""
	logger.info("writing the analyse output file with all the data counted from the chipseq data")
	results_dict = {}
	results_dict2 = {}
	results_dict3 = {}
	header = ['motif']
	analyse_output_name = os.path.join(args.output_directory, "precision4.txt")
	analyse_output_name2 = os.path.join(args.output_directory, "f1score4.txt")
	analyse_output_name3 = os.path.join(args.output_directory, "sensitivity4.txt")
	analyse_output = open(analyse_output_name, 'w')
	analyse_output2 = open(analyse_output_name2, 'w')
	analyse_output3 = open(analyse_output_name3, 'w')
	for dirname, dirnames, filenames in os.walk(args.output_directory):
		for directory in dirnames:
			header.append(directory)
			logger.info(directory)
			for filename in os.listdir(os.path.join(args.output_directory, directory)):
				if filename.endswith(".bed"):
					#logger.info("working with dir " + directory)
					motif_name = os.path.splitext(filename)[0]
					search_chipseq = os.path.join(args.output_directory, "within_peaks_" + motif_name.upper() + ".bed")
					tool_output_file = os.path.join(args.output_directory, directory, filename) #get the path to acess the file
					
					if os.stat(tool_output_file).st_size == 0: #the file is empty if fimo has not found any data
						sensitivity = 0
						precision = 0
						f1score = 0
							
					else:
						tp = int(subprocess.getoutput("bedtools intersect -a " + search_chipseq + " -b " + tool_output_file + " -u | wc -l")) #-u report original entry in chipseq_within_peaks if at least one overlap was found in mytool output file
						fn = int(subprocess.getoutput("bedtools intersect -a " + search_chipseq + " -b " + tool_output_file + " -v | wc -l")) #-v report original entry fom chipseq_within_peaks where no corresponding entry in mytool output file was found
						fp = int(subprocess.getoutput("bedtools intersect -a " + tool_output_file + " -b " + search_chipseq + " -v | wc -l"))
						sensitivity = tp/(tp + fn)
						precision = tp/(tp + fp)
						f1score = (2*tp)/(2*tp + fp + fn)
					if not motif_name in results_dict.keys(): #if this motif_name has no entry in the dictionary yet
						
						results_dict3[motif_name] = [str(sensitivity)]
						results_dict[motif_name] = [str(precision)]
						results_dict2[motif_name] = [str(f1score)]
					else: #append the next sensitivity to the list of sensitivities corresponding to the motif
						results_dict3[motif_name].append(str(sensitivity))
						results_dict[motif_name].append(str(precision))
						results_dict2[motif_name].append(str(f1score))
					#logger.info("the sensitivity for " + motif_name + " using " + directory.replace("_", " ") + " is " + str(sensitivity) + "\n")
					
	#print the output file
	
	analyse_output.write('\t'.join(header) + '\n')
	for entry in results_dict:
		line_to_write = [entry]
		line_to_write.extend(results_dict[entry])
		analyse_output.write('\t'.join(line_to_write) + '\n')
	analyse_output.close()
	analyse_output2.write('\t'.join(header) + '\n')
	for entry in results_dict2:
		line_to_write = [entry]
		line_to_write.extend(results_dict2[entry])
		analyse_output2.write('\t'.join(line_to_write) + '\n')
	analyse_output2.close
	
	analyse_output3.write('\t'.join(header) + '\n')
	for entry in results_dict3:
		line_to_write = [entry]
		line_to_write.extend(results_dict3[entry])
		analyse_output3.write('\t'.join(line_to_write) + '\n')
	analyse_output3.close()	
	"""
	#---------------make arrays from the analyse_output

	#read the analyse_output.txt and save the information from it to corresponding arrays to plot them later
	#analyse_output_name = os.path.join(args.output_directory, "analyse_output4.txt")
	#read_analyse_and_make_plots(analyse_output_name)

	analyse_output_name = os.path.join(args.output_directory, "precision4.txt")
	analyse_output_name2 = os.path.join(args.output_directory, "f1score4.txt")
	analyse_output_name3 = os.path.join(args.output_directory, "sensitivity4.txt")

	read_analyse_and_make_plots(analyse_output_name3)
		
	
	end_time = timeit.default_timer()
	running_time = end_time - start_time
	logger.info("the running time is: " + str(running_time) + " seconds")

if __name__ == "__main__":
main()
