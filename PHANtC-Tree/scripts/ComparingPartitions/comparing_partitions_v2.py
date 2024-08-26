#!/usr/bin/env	python3


"""
Automated run of ComparingPartitions (May 2021) to identify regions of cluster stability (i.e. cg/wgMLST partition ranges in which cluster composition is similar)
original github repository: https://github.com/jacarrico/ComparingPartitions
By Veronica Mixao
@INSA
comparing_partitions_v2.py runs one of two options: 
A) comparing partitions of a single method to determine stability regions (where stability regions correspond to partition threshold ranges where 
at least 5 subsequent partitions have an Adjusted Wallace above 0.99) 
comparing_partitions_v2.py -i1 PARTITIONS -o1 0 -a stability -t OUTPUT_NAME -n 5 -thr 0.99 -log LOG
B) comparing partitions between two different methods to determine their congruence (STILL UNDER DEVELOPMENT... RUN AT YOUR OWN RISK)
comparing_partitions_v2.py -i1 PARTITIONS -i2 PARTITIONS2 -o1 0 -o2 0 -a between_methods -t OUTPUT_NAME -n 5 -thr 0.99 -log LOG
"""


import sys
import math
import pandas 
import numpy as np
import argparse
from metrics import getContTable, getContTableTotals, getMismatchMatrix, getAdjRand, getSimpsons, getWallace, getAdjustedWallace
import textwrap
import random


# functions	----------

def rm_redundant(matrix, log):
	""" determine the most discriminatory cut-off
	and remove redundant samples """
	
	# what is the most discriminatory column?
	
	max_size,main_col = 0,""
	for col in matrix.columns[1:]:
		size = len(set(matrix[col].values.tolist()))
		if size > max_size:
			max_size,main_col = size,col
	print("\tThe most discriminatory column found was " + main_col)
	print("\tThe most discriminatory column found was " + main_col, file = log)
	
	# get samples
	
	final_samples = []
	for cluster in set(matrix[main_col].values.tolist()):
		flt_data = matrix[matrix[main_col] == cluster]
		selected_sample = random.choice(flt_data[flt_data.columns[0]].values.tolist())
		final_samples.append(selected_sample)
	
	if len(final_samples) != max_size:
		print("\tSampling got more samples than clusters... something went wrong :-(")
		print("\tSampling got more samples than clusters... something went wrong :-(", file = log)
		sys.exit()
	
	# final matrix
	
	matrix = matrix[matrix[matrix.columns[0]].isin(final_samples)]
	
	return matrix	
			
	
def	order_columns(matrix):
	""" Order the columns """
	
	cols = matrix.columns.tolist()
	matrix_ordered = pandas.concat([matrix[cols[0]], matrix[reversed(cols[1:])]], axis=1, join="inner")
	
	return matrix_ordered
	

def	metrics4stability(matrix, tag, log):
	""" Comparing partitions of a single method """
	
	with open(tag + "_metrics.tsv", "w+") as out:
		
		#prepare output
		
		cols = ["col2", "col1", "AdjRand", "Simpson_col2", "Simpson_col1", "Wallace", "AdjWallace", "AdjWallaceIClow", "AdjWallaceIChigh", "col2_cluster", "col1_cluster"]
		print("\t".join(cols), file = out)
		
		# calculating metrics
		
		i = 0
			
		for col1 in matrix.columns:
			i += 1
			j = 0
			for col2 in matrix.columns:
				j += 1
				if col1 != col2:
					if j-i == 1:
						print("\t\t",col2," => ",col1)
						icol = matrix[col1].values.tolist()
						jcol = matrix[col2].values.tolist()
						cont = getContTable(icol, jcol)
						cont_2 = getContTable(jcol, icol)
						abcdn = getMismatchMatrix(cont, icol, jcol)
								
						sid1 = getSimpsons(icol)
						sid2 = getSimpsons(jcol)
						
						wallace = getWallace(abcdn[0], abcdn[1], abcdn[2])
						
						AdjustedWallace = getAdjustedWallace(cont_2, jcol, icol)
						print("\t\t\tAdjusted Wallace ",col2," vs ",col1,"\t%f\t%f\t%f" % AdjustedWallace)
						print("\t\t\tAdjusted Wallace ",col2," vs ",col1,"\t%f\t%f\t%f" % AdjustedWallace, file = log)
						
						print("\t".join([str(col2), str(col1), str(getAdjRand(icol,jcol)), str(sid2[0]), str(sid1[0]), str(wallace[0]), str(AdjustedWallace[0]), str(AdjustedWallace[1]), str(AdjustedWallace[2]), str(sid2[3]), str(sid1[3])]), file = out)


def	comparing_methods(matrix1, matrix2, tag, log):
	""" Comparing partitions of two methods """
	
	# defining dictionaries for df

	AdjWallace1 = {} #dictionary with all the Adjusted Wallace values from file1 to file2 (keys = file2)
	AdjWallace2 = {} #dictionary with all the Adjusted Wallace values from file2 to file1 (keys = file1)
	Wallace1 = {} #dictionary with all the Wallace values from file1 to file2 (keys = file2)
	Wallace2 = {} #dictionary with all the Wallace values from file2 to file1 (keys = file1)
	Simpsons1 = {} #dictionary with all the Simpsons values for file1 (keys = columns)
	Simpsons2 = {} #dictionary with all the Simpsons values for file2 (keys = columns)
	Rand = {} #dictionary with all the Rand metrics (keys = file2)
	final_metrics = {} #dictionary with the sum of Adjusted Wallace 1, Adjusted Wallace 2 and Adjusted Rand

	AdjWallace1["0_method_partition"] = []
	AdjWallace2["0_method_partition"] = []
	Wallace1["0_method_partition"] = []
	Wallace2["0_method_partition"] = []
	Rand["0_method_partition"] = [] 
	final_metrics["0_method_partition"] = [] 

	#determine all the columns for the dataframe
	
	for col1 in matrix1.columns:
		AdjWallace1["0_method_partition"].append(col1)
		Wallace1["0_method_partition"].append(col1)
		AdjWallace2["0_method_partition"].append(col1)
		Wallace2["0_method_partition"].append(col1)
		Rand["0_method_partition"].append(col1)
		final_metrics["0_method_partition"].append(col1)
		
	for col2 in matrix2.columns:
		AdjWallace1[col2] = []
		Wallace1[col2] = []
		Rand[col2] = []
		AdjWallace2[col2] = []
		Wallace2[col2] = []
		final_metrics[col2] = []

	#performing the comparisons
		
	for col1 in matrix1.columns:
		for col2 in matrix2.columns:
			print("\t\t",col1," ---------- ",col2)
			print("\t\t",col1," ---------- ",col2, file = log)
			icol = matrix1[col1].values.tolist()
			jcol = matrix2[col2].values.tolist()
			cont = getContTable(icol, jcol)
			cont_2 = getContTable(jcol, icol)
			abcdn = getMismatchMatrix(cont, icol, jcol)
			abcdn_2 = getMismatchMatrix(cont_2, jcol, icol)
				
			#getting rand metric
			rand_metric = getAdjRand(icol,jcol)
			Rand[col2].append(rand_metric)
				
			#getting simpson's metric
			sid1 = getSimpsons(icol)
			sid2 = getSimpsons(jcol)
				
			Simpsons1[col1] = sid1
			Simpsons2[col2] = sid2 
				
			#getting Wallace metric
			wallace = getWallace(abcdn[0], abcdn[1], abcdn[2])
				
			Wallace1[col2].append(wallace[0])
			Wallace2[col2].append(wallace[0])
						
			#getting Adjusted Wallace metric
			AdjustedWallace = getAdjustedWallace(cont, icol, jcol)
			AdjustedWallace_2 = getAdjustedWallace(cont_2, jcol, icol)
			
			AdjWallace1[col2].append(AdjustedWallace[0])
			AdjWallace2[col2].append(AdjustedWallace_2[0])
			
			#getting final metrics
			final_score = rand_metric + AdjustedWallace[0] + AdjustedWallace_2[0]
			
			final_metrics[col2].append(final_score)

	df_Rand = pandas.DataFrame(data = Rand)
	df_Simpsons1 = pandas.DataFrame(data = Simpsons1) #needs to be transposed
	df_Simpsons2 = pandas.DataFrame(data = Simpsons2) #needs to be transposed
	df_Wallace1 = pandas.DataFrame(data = Wallace1)
	df_Wallace2 = pandas.DataFrame(data = Wallace2)
	df_AdjWallace1 = pandas.DataFrame(data = AdjWallace1)
	df_AdjWallace2 = pandas.DataFrame(data = AdjWallace2)
	df_final_metrics = pandas.DataFrame(data = final_metrics)

	df_Simpsons1T = df_Simpsons1.T
	df_Simpsons2T = df_Simpsons2.T

	df_Rand.to_csv(tag + "_AdjustedRand.tsv", index = False, header= True, sep ="\t")
	df_Simpsons1T.to_csv(tag + "_Simpsons1.tsv", index = True, header= False, sep ="\t")
	df_Simpsons2T.to_csv(tag + "_Simpsons2.tsv", index = True, header= False, sep ="\t")
	df_Wallace1.to_csv(tag + "_Wallace1.tsv", index = False, header= True, sep ="\t")
	df_Wallace2.to_csv(tag + "_Wallace2.tsv", index = False, header= True, sep ="\t")
	df_AdjWallace1.to_csv(tag + "_AdjWallace1.tsv", index = False, header= True, sep ="\t")
	df_AdjWallace2.to_csv(tag + "_AdjWallace2.tsv", index = False, header= True, sep ="\t")
	df_final_metrics.to_csv(tag + "_final_score.tsv", index = False, header= True, sep ="\t")


def	get_stability(infile, thr, n, tag, log):
	""" Determine stability regions of a single method """
	
	print("\t\tDetermining stability regions with threshold: ", str(thr), " and minimum n of ", str(n))
	print("\t\tDetermining stability regions with threshold: ", str(thr), " and minimum n of ", str(n), file = log)
	
	with open(infile, "r") as f_open:
		with open(tag + "_stableRegions.tsv", "w+") as out:
			f = f_open.readlines()
			i = 0
			j = 0
			counter = 0
			stable_blocks = {}
			
			for line in f:
				i += 1
				if i > 1:
					l = line.split("\t")
					score = float(l[6])
					
					if score >= float(thr):
						info = l[0] + "->" + l[1]
						if j == 0:
							j += 1
							counter += 1
							stable_blocks["block_" + str(counter)] = []
						stable_blocks["block_" + str(counter)].append(info)
					else:
						j = 0
			
			print("#Stability regions of the Adjusted Wallace coefficient for ", infile, "\n#using a threshold of ", str(thr), " and minimum number of required observations of ", str(n), file = out)
			print("#block_id\tfirst_partition\tlast_partition\tlen_block", file = out)
			
			for block in stable_blocks.keys():
				if len(stable_blocks[block]) >= int(n):
					print("\t".join([block, stable_blocks[block][0], stable_blocks[block][-1], str(len(stable_blocks[block]))]), file = out)
					print("\t\t", block, stable_blocks[block][0], stable_blocks[block][-1], str(len(stable_blocks[block])))
					print("\t\t", block, stable_blocks[block][0], stable_blocks[block][-1], str(len(stable_blocks[block])), file = log)


# pipeline ----------
	
def main():
	
	# define options
	parser = argparse.ArgumentParser(prog="comparing_partitions_v2.py", formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent("""\
									###############################################################################             
									#                                                                             #
									#                         comparing_partitions_v2.py                          #
									#                                                                             #
									###############################################################################  
									                            
									comparing_partitions_v2 identifies regions of cluster stability (i.e. cg/wgMLST 
									partition threshold ranges in which cluster composition is similar)
									
									
									This is a modified version of:
									https://github.com/jacarrico/ComparingPartitions
									
									
									By default, for the stability analysis, this script searches for the most 
									discriminatory partition and only keeps one random samples per cluster. This 
									option can be reversed with the option '--keep-redundants'.
									comparing_partitions_v2.py runs one of two options: 
									
									A) comparing partitions of a single method to determine stability regions (where 
									stability regions correspond to partition threshold ranges where at least 5 
									subsequent partitions have an Adjusted Wallace above 0.99) 
									comparing_partitions_v2.py -i1 PARTITIONS -o1 0 -a stability -t OUTPUT_NAME -n 5 
									-thr 0.99 -log LOG
									B) comparing partitions between two different methods to determine their congruence
									comparing_partitions_v2.py -i1 PARTITIONS -i2 PARTITIONS2 -o1 0 -o2 0 -a 
									between_methods -t OUTPUT_NAME -n 5 -thr 0.99 -log LOG
									
									WARNING!! Cluster congruence analysis between different methods is still under 
									development in the frame of the One Health EJP BeOne project for a major 
									comparison of surveillance pipelines! So, for now, use it at your own risk!
									
									-------------------------------------------------------------------------------"""))
	
	parser.add_argument("-i1", "--input1", dest="input1", action= "store", required=True, help="[MANDATORY] Input matrix 1 (table with partitions)")
	parser.add_argument("-i2", "--input2", dest="input2", action= "store", help="Input matrix 2 (table with partitions)")
	parser.add_argument("-o1", "--order1", dest="order1", action= "store", default=0, required=True, help="Partitions order in matrix 1 (0: min -> max; 1: max -> min) [0]")
	parser.add_argument("-o2", "--order2", dest="order2", action= "store", help="Partitions order in matrix 2 (0: min -> max; 1: max -> min) [0]")
	parser.add_argument("-a", "--analysis", dest="analysis", action= "store", default="between_methods", help="Type of analysis (options: stability or between_methods) [between_methods]")
	parser.add_argument("-t", "--tag", dest="tag", action= "store", required=True, help="[MANDATORY] Tag for output file name")
	parser.add_argument("-n", "--n_obs", dest="n_obs", action="store", default=5, help="Minimum number of sequencial observations to consider an interval for method stability analysis [5]")
	parser.add_argument("-thr", "--threshold", dest="threshold", action= "store", default=0.99, help="Threshold of Adjusted Wallace score to consider an observation for method stability \
						analysis [0.99]")
	parser.add_argument("--keep-redundants", dest="keep_redundants", action= "store_true", help="Set ONLY if you want to keep all samples of each cluster of the most discriminatory partition\
						 (by default redundant samples are removed to avoid the influence of cluster size)")
	parser.add_argument("-log", "--log", dest="log", action= "store", default="log", help="Log file")
						
	args = parser.parse_args()
	
	
	# starting logs	----------

	log_name = args.log
	log = open(log_name, "a+")
	
	print("\n-------------------- comparing_partitions_v2.py --------------------\n")
	print("\n-------------------- comparing_partitions_v2.py --------------------\n", file = log)
	print(" ".join(sys.argv))
	print(" ".join(sys.argv), file = log)
	
	
	# pipeline comparing methods 	----------
	
	if args.analysis == "between_methods":
		if not args.input2 or not args.order2:
			print("Provided options are incompatible! 'between_method' analysis requires input1, input2, order1 and order2!")
			print("Provided options are incompatible! 'between_method' analysis requires input1, input2, order1 and order2!", file = log)
			sys.exit()
		else:
			matrix1 = pandas.read_table(args.input1)
			matrix2 = pandas.read_table(args.input2)
		print("Preparing for comparing methods...")
		print("Preparing for comparing methods...", file = log)
		
		# remove redundant samples
		
		if not args.keep_redundants: # needs to remove redundant samples 
			print("Removing redundant samples...")
			print("Removing redundant samples...", file = log)
			matrix1 = rm_redundant(matrix1, log)
			matrix2 = rm_redundant(matrix2, log)
		
		
		# ordering
		
		print("Ordering matrix...")
		print("Ordering matrix...", file = log)
		
		if int(args.order1) != 0:
			if int(args.order1)  == 1:
				ordered_matrix1 = order_columns(matrix1)
			else:
				print("\tInvalid order of matrix1!")
				print("\tInvalid order of matrix1!", file = log)
				sys.exit()
		else:
			ordered_matrix1 = matrix1
		
		if int(args.order2) != 0:
			if int(args.order2) == 1:
				ordered_matrix2 = order_columns(matrix2)
			else:
				print("\tInvalid order of matrix1!")
				print("\tInvalid order of matrix1!", file = log)
				sys.exit()
		else:
			ordered_matrix2 = matrix2
		
		# filter matrix
		
		print("Filtering matrix...")
		print("Filtering matrix...", file = log)
		
		filtered_matrix1_not_ordered = ordered_matrix1.loc[ordered_matrix1[ordered_matrix1.columns[0]].isin(ordered_matrix2[ordered_matrix2.columns[0]])]
		filtered_matrix2_not_ordered = ordered_matrix2.loc[ordered_matrix2[ordered_matrix2.columns[0]].isin(ordered_matrix1[ordered_matrix1.columns[0]])]
		
		filtered_matrix1 = filtered_matrix1_not_ordered.sort_values(by=filtered_matrix1_not_ordered.columns[0])
		filtered_matrix2 = filtered_matrix2_not_ordered.sort_values(by=filtered_matrix2_not_ordered.columns[0])
		
		# cleaning
		
		print("Cleaning matrix...")
		print("Cleaning matrix...", file = log)
		
		clean_matrix1 = filtered_matrix1.drop(filtered_matrix1.columns[0], axis=1)
		clean_matrix2 = filtered_matrix2.drop(filtered_matrix2.columns[0], axis=1)

		# removing columns in which all samples form a cluster 

		for col in clean_matrix1.columns:
			if len(pandas.unique(clean_matrix1[col])) == 1:
				clean_matrix1 = clean_matrix1.drop(col, axis=1)
		for col in clean_matrix2.columns:
			if len(pandas.unique(clean_matrix2[col])) == 1:
				clean_matrix2 = clean_matrix2.drop(col, axis=1)

		print("\tN partitions matrix 1:", len(clean_matrix1.columns))
		print("\tN partitions matrix 2:", len(clean_matrix2.columns))

		# comparing methods
		
		print("Comparing methods...")
		print("Comparing methods...", file = log)
		
		comparing_methods(clean_matrix1, clean_matrix2, args.tag, log)
		
			
	# pipeline stability
	
	elif args.analysis == "stability":
		if args.input2 or args.order2:
			print("Provided options are incompatible! 'stability' analysis only accepts input1 and order1!")
			print("Provided options are incompatible! 'stability' analysis only accepts input1 and order1!", file = log)
			sys.exit()
		else:
			matrix = pandas.read_table(args.input1)
		print("Preparing for stability analysis...")
		print("Preparing for stability analysis...", file = log)
		
		# remove redundant samples
		
		if not args.keep_redundants: # needs to remove redundant samples 
			print("Removing redundant samples...")
			print("Removing redundant samples...", file = log)
			matrix = rm_redundant(matrix, log)
			
		# ordering
		
		print("Ordering matrix...")
		print("Ordering matrix...", file = log)
		
		if int(args.order1) != 0:
			if int(args.order1)  == 1:
				ordered_matrix = order_columns(matrix)
			else:
				print("\tInvalid order!")
				print("\tInvalid order!", file = log)
				sys.exit()
		else:
			ordered_matrix = matrix
		
		# cleaning
		
		print("Cleaning matrix...")
		print("Cleaning matrix...", file = log)
		
		clean_matrix = ordered_matrix.drop(ordered_matrix.columns[0], axis=1)
			
		# calculating metrics
		
		print("Calculating metrics...")
		print("Calculating metrics...", file = log)
		
		metrics4stability(clean_matrix, args.tag, log)
		
		# get stability
		
		print("Get stability...")
		print("Get stability...", file = log)
		
		get_stability(args.tag + "_metrics.tsv", args.threshold, args.n_obs, args.tag, log)
	
	
	# other options
	
	else:
		print("Invalid analysis type!")
		print("Invalid analysis type!", file = log)
		sys.exit()
	
	
	print("\ncomparing_partitions_v2.py is done!")
	print("\ncomparing_partitions_v2.py is done!", file = log)
	
	log.close()

# running pipeline
	
if __name__ == '__main__':
	main()
