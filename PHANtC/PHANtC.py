#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
############################################################################
# Istituto Superiore di Sanita'
# European Union Reference Laboratory (EU-RL) for Escherichia coli, including Verotoxigenic E. coli (VTEC)
# Developer: Arnold Knijn arnold.knijn@iss.it
############################################################################
"""

import argparse
import configparser
import sys
import os
import fileinput
import getopt
import zlib
import textwrap
import shutil

TOOL_DIR = os.path.dirname(os.path.abspath(__file__))

def fasta_reader(filename):
  # read fasta file and return records
  from Bio.SeqIO.FastaIO import FastaIterator
  with open(filename) as handle:
    for record in FastaIterator(handle):
      yield record

def get_string_hash(string):
  return str(zlib.crc32(string.encode('utf8')))

def elaborate_fastafile(schemadirectory, file, outputdirectory):
  filename, file_extension = os.path.splitext(file)
  strhashnumberfile = filename + "_hash_number.txt"
  with open(os.path.join(outputdirectory, strhashnumberfile), "w") as hash_number:
    for entry in fasta_reader(os.path.join(schemadirectory, file)):
      # write allele hash to hash_number file
      hash_number.write(get_string_hash(str(entry.seq)))
      hash_number.write('\t')
      # write allele number to hash_number file
      hash_number.write(str(entry.id)[str(entry.id).rfind('_')+1:])
      hash_number.write('\n')

def get_hash_from_number(locus_hash_number_file, allele_number):
  if allele_number == "0" or allele_number in ['EXC', 'PLOT3', 'PLOT5', 'LOTSC', 'NIPH', 'NIPHEM', 'ALM', 'ASM', 'LNF']:
    return "0"
  else:
    if allele_number[0:4] == 'INF-':
      allele_number = allele_number[4:]
    # get the allele hash from the corresponding allele number
    with open(locus_hash_number_file) as locus_hashes_numbers:
      for allele_hash_number in locus_hashes_numbers:
        hash = 0
        hash, number = allele_hash_number.split('\t')
        if number.rstrip() == allele_number or number.rstrip() == '*' + allele_number:
          return hash

def get_number_from_hash(locus_hash_number_file, allele_hash):
  if allele_hash == "0" or allele_hash == "-":
    return "0"
  else:
    # get the allele number from the corresponding allele hash
    with open(locus_hash_number_file) as locus_hashes_numbers:
      for allele_hash_number in locus_hashes_numbers:
        hash, number = allele_hash_number.split('\t')
        if hash == allele_hash:
          return number.rstrip()

def elaborate_sample(samplefile, schemadirectory, outputdirectory):
  # a sample file was given, substitute the allele numbers with the allele sequence hashes
  samplefilename = os.path.basename(samplefile)
  with open(samplefile) as sample, open(os.path.join(outputdirectory, samplefilename), "w") as unhashed_sample:
    line_sample = sample.readline()
    header_sample = line_sample.split('\t')
    numcols = len(header_sample)
    unhashed_sample.write(line_sample)
    hashed_allele = [None] * numcols
    line_sample = sample.readline().split('\t')
    while len(line_sample[0]) > 0:
      unhashed_sample.write(line_sample[0].rstrip())
      for i in range(1, numcols):
        unhashed_sample.write('\t')
        unhashed_sample.write(str(get_number_from_hash(os.path.join(schemadirectory, header_sample[i].rstrip()[:-6] + "_hash_number.txt"), line_sample[i].rstrip()) or 'LNF'))
      unhashed_sample.write('\n')
      line_sample = sample.readline().split('\t')

def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--species', dest='species', help='sample species')
    parser.add_argument('--input', dest='input', help='samples')
    parser.add_argument('--hashprofiles', dest='hashprofiles', help='CRC32 hashes MLST profiles')
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(TOOL_DIR + '/../phantastic.conf')
    dataflowdir = config['fs']['dataflow_path']

    if args.species == "Escherichia coli":
        schemadirectory = dataflowdir + "/Chewie-NS/ecoli/ecoli_INNUENDO_wgMLST_ORIG"
    elif args.species == "Listeria monocytogenes":
        schemadirectory = dataflowdir + "/Chewie-NS/listeria/lmonocytogenes_Pasteur_cgMLST_ORIG"
    else:
        schemadirectory = dataflowdir + "/Chewie-NS"
    outputdirectory = "outputdirectory"
    profilefile = "cgMLST.tmp"
    # elaborate_fastafile for all fasta files in the directory
    if not os.path.isdir(outputdirectory):
      try:
        os.mkdir(outputdirectory)
      except OSError:
        print ("Creation of the directory %s failed" % path)
    for file in os.listdir(schemadirectory):
      if file.endswith(".fasta"):
        elaborate_fastafile(schemadirectory, file, outputdirectory)
    # create a profile file with assigned numbers instead of the hashes
    shutil.copyfile(args.input, profilefile)
    elaborate_sample(profilefile, outputdirectory, outputdirectory)
    shutil.copyfile(outputdirectory+"/"+profilefile, args.hashprofiles)

if __name__ == "__main__":
    __main__()
