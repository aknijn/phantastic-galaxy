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
from phantdb import IridaDb

TOOL_DIR = os.path.dirname(os.path.abspath(__file__))

def get_coverage(input_id):
    with IridaDb("Shiga toxin-producing Escherichia coli") as iridadb:
        return iridadb.get_singleend_coverage(input_id[2:])

def get_fastp(json_file):
    with open(json_file, "rb") as fastp_json:
        fastp_dict = json.load(fastp_json)
        read_mean_length = fastp_dict['summary']['after_filtering']['read1_mean_length']
        q30_rate = fastp_dict['summary']['after_filtering']['q30_rate']
        total_bases = fastp_dict['summary']['after_filtering']['total_bases']
    return (read_mean_length, q30_rate, total_bases)

def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('-1', '--input1', dest='input1', help='forward or single-end reads file in Sanger FASTQ format')
    parser.add_argument('--input_id', dest='input_id', help='sampleId')
    parser.add_argument('--genomeSize', dest='genomeSize', help='genome size')
    parser.add_argument('--contigs', dest='contigs', help='contigs')
    parser.add_argument('--json', dest='json', help='json')
    parser.add_argument('--quast', dest='quast', help='quast')
    args = parser.parse_args()

    # if filename ends with .dat a fastq.gz file was decomrpessed
    if args.input1.endswith(".fastq") or args.input1.endswith(".dat"):
        # FASTQ
        subprocess.call("ln -s " + args.input1 + " fastq_in.fastqsanger", shell=True)
        if float(get_coverage(args.input_id)) < 100:
            # NO TRIMMING
            subprocess.call("fastp --thread 4 -i fastq_in.fastqsanger -o input_1.fq -f 0 -t 0 -l 5 --cut_front_window_size 0 --cut_front_mean_quality 1 --cut_tail_window_size 0 --cut_tail_mean_quality 1", shell=True)
        else:
            # TRIMMING
            subprocess.call("fastp --thread 4 -i fastq_in.fastqsanger -o input_1.fq -f 3 -t 3 -l 55 --cut_front_window_size 5 --cut_front_mean_quality 20 --cut_tail_window_size 5 --cut_tail_mean_quality 20", shell=True)
        # ASSEMBLY
        str_sequencer = ""
        # check if the file is ION Torrent
        with open('input_1.fq', 'r') as fq:
            if fq.readline().count(':') == 2:
                str_sequencer = "--iontorrent"
        subprocess.call("perl " + TOOL_DIR + "/scripts/spades.pl spades_contigs spades_contig_stats spades_scaffolds spades_scaffold_stats spades_log NODE spades.py --disable-gzip-output --isolate -t ${GALAXY_SLOTS:-16} " + str_sequencer + " -s input_1.fq", shell=True)
        subprocess.call("perl " + TOOL_DIR + "/scripts/filter_spades_repeats.pl -i spades_contigs -t spades_contig_stats -c 0.33 -r 1.75 -l 1000 -o output_with_repeats -u output_without_repeats -n repeat_sequences_only -e 5000 -f discarded_sequences -s summary", shell=True)
        shutil.move("output_without_repeats", args.contigs)
    else:
        # FASTA
        shutil.copy(args.input1, args.contigs)
      
    # QUAST
    genome_size_base = str(int(float(args.genomeSize) * 1000000))
    subprocess.call("quast --threads 4 -o outputdir --est-ref-size " + genome_size_base + " --min-contig 500 -l  '" + args.input_id + "' --contig-thresholds 0,1000 " + args.contigs, shell=True)
    shutil.move("outputdir/report.tsv", args.quast)

    # JSON
    report_data = {}
    report = open(args.json, 'w')
    # write JSON (json.dumps => TypeError: Decimal is not JSON serializable)
    #report_data["coverage"] = getCoverage(args.input1)
    #report.write(json.dumps(report_data))
    if args.input1.endswith(".fastq") or args.input1.endswith(".dat"):
        report.write("{\"coverage\": \"" + str(get_coverage(args.input_id)) + "\",")
        (read_mean_length, q30_rate, total_bases) = get_fastp('fastp.json')
        report.write("\"read_mean_length\": \"" + str(read_mean_length) + "\",")
        report.write("\"q30_rate\": \"" + str(q30_rate) + "\",")
        report.write("\"total_bases\": \"" + str(total_bases) + "\"}")
    else:
        report.write("{\"coverage\": \"ND\"}")
    report.close()

if __name__ == "__main__":
    __main__()
