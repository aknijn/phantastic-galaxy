#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
############################################################################
# Istituto Superiore di Sanita'
# European Union Reference Laboratory (EU-RL) for Escherichia coli, including Verotoxigenic E. coli (VTEC)
# Developer: Arnold Knijn arnold.knijn@iss.it
############################################################################
"""

import sys, getopt
import os

def mlst_calls(call_file):
  # return a MLST call matrix, with samples on the rows and loci on the columns.
  # header row is excluded and last two columns (CC, clonal_complex) are excluded
  with open(call_file) as file_in:
    lines = []
    for line in file_in:
      lines.append(line.split('\t')[:-2])
  data_mlst = lines[1:]
  return data_mlst


def compare_alleles(allele_1, allele_2):
  # compare alleles, return 1 if different, 0 if equal or not to be compared
  if allele_1 == allele_2:
    comparison = 0
  elif allele_1 == '0' or allele_2 == '0':
    # allele not found
    comparison = 0
  else:
    comparison = 1
  return comparison


def mlst_distance(mlst):
  # a profile file was given, substitute the allele numbers with the allele sequence hashes
  rows = len(mlst)
  cols = len(mlst[0])
  D = [ [0]*(rows) for _ in range(rows) ]
  h = []
  for row in range(0, rows):
    h.append(mlst[row][0])
    for row2 in range(row+1, rows):
      dist = 0
      for col in range(1, cols):
        dist = dist + compare_alleles(mlst[row][col], mlst[row2][col])
      D[row][row2] = dist
      D[row2][row] = dist
  return D, h

  
def main(argv):
  input = ''
  strusage = 'mlst_distance.py -i <input.tsv> -o <output.tsv>\n'
  numloci = 0
  try:
    opts, args = getopt.getopt(argv,"hi:o:",["input=","output="])
  except getopt.GetoptError:
    print (strusage)
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print (strusage)
      sys.exit()
    elif opt in ("-i", "--input"):
      input = arg
    elif opt in ("-o", "--output"):
      output = arg
  if os.path.isfile(input):
    print ('input file is "', input, '"')
  else: 
    print ('input file is "', input, '" but does not exist')
    sys.exit(0)
  print ('output file is "', output, '"')
  mlst = mlst_calls(input)
  dist_mat, head_mat = mlst_distance(mlst)
  header = '\t'.join(head_mat)
  with open(output, "w") as output:
    output.write('\t' + header + '\n')
    i = 0
    for row in dist_mat:
      output.write(head_mat[i])
      i = i + 1
      for elem in row:
        output.write('\t' + str(elem))
      output.write('\n')
    

if __name__ == "__main__":
   main(sys.argv[1:])

