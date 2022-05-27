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
import shutil
import subprocess
import json
import datetime
import mysql.connector
from mysql.connector import errorcode

TOOL_DIR = os.path.dirname(os.path.abspath(__file__))

def getIdFile(filename):
    splitFilename = filename.split("/")
    return splitFilename[5]

def getCoverage(inputfile):
    config = configparser.ConfigParser()
    config.read(TOOL_DIR + '/../phantastic.conf')
    dbhost = config['db']['host']
    dbdatabase = config['db']['database']
    dbuser = config['db']['user']
    dbpassword = config['db']['password']
    config = {
        'user': dbuser,
        'password': dbpassword,
        'host': dbhost,
        'database': dbdatabase
    }

    files_id = getIdFile(inputfile)
    sql = ("SELECT total_bases/genome_size from sequence_file_pair_files AS sfpf " +
          "inner join sample_sequencingobject AS sso on sfpf.pair_id = sso.sequencingobject_id " +
          "inner join  qc_entry  AS so on sso.sequencingObject_id = so.sequencingObject_id " +
          "inner join project_sample AS ps on sso.sample_id = ps.sample_id " +
          "inner join project AS p on ps.project_id = p.id " +
          "where genome_size>0 and sfpf.files_id = " + files_id)

    try:
        cnx = mysql.connector.connect(**config)
        cursor = cnx.cursor(buffered=True)
        cursor.execute(sql)
        row = cursor.fetchone()
        cursor.close()
        return row[0]
    except mysql.connector.Error as err:
        print(err)
    else:
        cnx.close()

def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('-1', '--input1', dest='input1', help='forward or single-end reads file in Sanger FASTQ format')
    parser.add_argument('--input_id', dest='input_id', help='sampleId')
    parser.add_argument('--genomeSize', dest='genomeSize', help='genome size')
    parser.add_argument('--contigs', dest='contigs', help='contigs')
    parser.add_argument('--json', dest='json', help='json')
    parser.add_argument('--quast', dest='quast', help='quast')
    args = parser.parse_args()

    report_data = {}
    report = open(args.json, 'w')
    # write JSON (json.dumps => TypeError: Decimal is not JSON serializable)
    #report_data["coverage"] = getCoverage(args.input1)
    #report.write(json.dumps(report_data))
    if args.input1.endswith(".fastq"):
      report.write("{\"coverage\": \"" + str(getCoverage(args.input1)) + "\"}")
    else:
      report.write("{\"coverage\": \"ND\"}")
    report.close()

    if args.input1.endswith(".fastq"):
      # FASTQ
      if float(getCoverage(args.input1)) < 100:
          # NO TRIMMING
          subprocess.call("ln -s " + args.input1 + " input_1.fq", shell=True)
      else:
          # TRIMMING
          subprocess.call("ln -s " + args.input1 + " fastq_in.fastqsanger", shell=True)
          os.system("ln -s " + os.popen("which trimmomatic.jar").read().strip() + " trimmomatic.jar")
          subprocess.call("java ${_JAVA_OPTIONS:--Xmx8G} -jar trimmomatic.jar SE -threads ${GALAXY_SLOTS:-6} -phred33 fastq_in.fastqsanger input_1.fq SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:55", shell=True)
      # ASSEMBLY
      strSequencer = ""
      # check if the file is ION Torrent
      with open('input_1.fq', 'r') as fq:
          if fq.readline().count(':') == 2:
              strSequencer = "--iontorrent"
      subprocess.call("perl " + TOOL_DIR + "/scripts/spades.pl spades_contigs spades_contig_stats spades_scaffolds spades_scaffold_stats spades_log NODE spades.py --disable-gzip-output --isolate -t ${GALAXY_SLOTS:-16} " + strSequencer + " -s input_1.fq", shell=True)
      subprocess.call("perl " + TOOL_DIR + "/scripts/filter_spades_repeats.pl -i spades_contigs -t spades_contig_stats -c 0.33 -r 1.75 -l 1000 -o output_with_repeats -u output_without_repeats -n repeat_sequences_only -e 5000 -f discarded_sequences -s summary", shell=True)
      shutil.move("output_without_repeats", args.contigs)
    else:
      # FASTA
      shutil.copy(args.input1, args.contigs)
      
    # QUAST
    genomeSizeBase = str(int(float(args.genomeSize) * 1000000))
    subprocess.call("quast --threads 4 -o outputdir --est-ref-size " + genomeSizeBase + " --min-contig 500 -l  '" + args.input_id + "' --contig-thresholds 0,1000 " + args.contigs, shell=True)
    shutil.move("outputdir/report.tsv", args.quast)

if __name__ == "__main__":
    __main__()
