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
import shutil
import subprocess
import json
import datetime

TOOL_DIR = os.path.dirname(os.path.abspath(__file__))

def get_filetype(input_file):
    filetype = "error"
    with open(input_file) as infile:
        for line in infile:
            if line in ['\n', '\r\n']:
                continue
            elif line[0] == '@':
                filetype = "fastq"
                break
            elif line[0] == '>':
                filetype = "fasta"
                break
            else:
                filetype = "other"
                break
    return filetype

def get_coverage(input_fastq, genome_size):
    str_cmd = "cat " + input_fastq + " | paste - - - - | cut -f 2 | tr -d '\n' | wc -c"
    genomeSize = subprocess.run(str_cmd, shell=True, capture_output=True)
    coverage = "{:.2f}".format( float(genomeSize.stdout.decode("utf-8").strip()) / (float(genome_size)*1000000) )
    return coverage

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

    str_coverage = "ND"
    is_fastq = False
    if get_filetype(args.input1) == "fastq":
        is_fastq = True
    if is_fastq:
        # FASTQ
        os.symlink(args.input1, 'fastq_in.fastqsanger')
        str_coverage = str(get_coverage(args.input1, args.genomeSize))
        if float(str_coverage) < 100:
            # NO TRIMMING
            subprocess.run("fastp --thread 4 -i fastq_in.fastqsanger -o input_1.fq -f 0 -t 0 -l 5 --cut_front_window_size 0 --cut_front_mean_quality 1 --cut_tail_window_size 0 --cut_tail_mean_quality 1", shell=True)
        else:
            # TRIMMING
            subprocess.run("fastp --thread 4 -i fastq_in.fastqsanger -o input_1.fq -f 3 -t 3 -l 55 --cut_front_window_size 5 --cut_front_mean_quality 20 --cut_tail_window_size 5 --cut_tail_mean_quality 20", shell=True)
        # ASSEMBLY
        str_sequencer = ""
        # check if the file is ION Torrent
        with open('input_1.fq', 'r') as fq:
            if fq.readline().count(':') == 2:
                str_sequencer = "--iontorrent"
        subprocess.run("perl " + TOOL_DIR + "/scripts/spades.pl spades_contigs spades_contig_stats spades_scaffolds spades_scaffold_stats spades_log NODE spades.py --disable-gzip-output --isolate -t ${GALAXY_SLOTS:-16} " + str_sequencer + " -s input_1.fq", shell=True)
        subprocess.run("perl " + TOOL_DIR + "/scripts/filter_spades_repeats.pl -i spades_contigs -t spades_contig_stats -c 0.33 -r 1.75 -l 1000 -o output_with_repeats -u output_without_repeats -n repeat_sequences_only -e 5000 -f discarded_sequences -s summary", shell=True)
        shutil.move("output_without_repeats", args.contigs)
    else:
        # FASTA
        shutil.copy(args.input1, args.contigs)
      
    # QUAST
    genome_size_base = str(int(float(args.genomeSize) * 1000000))
    subprocess.run("quast --threads 4 -o outputdir --est-ref-size " + genome_size_base + " --min-contig 500 -l  '" + args.input_id + "' --contig-thresholds 0,1000 " + args.contigs, shell=True)
    shutil.move("outputdir/report.tsv", args.quast)

    # JSON
    report_data = {}
    report = open(args.json, 'w')
    # write JSON (json.dumps => TypeError: Decimal is not JSON serializable)
    #report_data["coverage"] = getCoverage(args.input1)
    #report.write(json.dumps(report_data))
    if is_fastq:
        report.write("{\"coverage\": \"" + str_coverage + "\",")
        (read_mean_length, q30_rate, total_bases) = get_fastp('fastp.json')
        report.write("\"read_mean_length\": \"" + str(read_mean_length) + "\",")
        report.write("\"q30_rate\": \"" + str(q30_rate) + "\",")
        report.write("\"total_bases\": \"" + str(total_bases) + "\"}")
    else:
        report.write("{\"coverage\": \"ND\"}")
    report.close()

if __name__ == "__main__":
    __main__()
