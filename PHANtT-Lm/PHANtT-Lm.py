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
import sys
import os
import shutil
import subprocess
import json
import datetime
from pathlib import Path

TOOL_DIR = os.path.dirname(os.path.abspath(__file__))

def openFileAsTable(filename):
    with open(filename) as table_in:
        table_data = [[str(col).rstrip() for col in row.split('\t')] for row in table_in]
    return table_data

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

def getAmplicons(sero_typing_tab):
    strAmplicon = "-"
    for i in [3,4,5,6,2]:
        if sero_typing_tab[1][i] == "FULL":
            strAmplicon += "," + sero_typing_tab[0][i]
    return strAmplicon.replace("-,", "")

def getCCandLineage(mlst_ST):
    strCC = "?"
    strLineage = "?"
    lst_types = openFileAsTable(TOOL_DIR + "/data/lmonocytogenes.txt")
    for lst_type in lst_types:
        if (mlst_ST == lst_type[0]):
            strCC = lst_type[8]
            strLineage = lst_type[9]
            break
    return strCC, strLineage

def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('-1', '--input1', dest='input1', help='forward or single-end reads file in Sanger FASTQ format')
    parser.add_argument('-2', '--input2', dest='input2', help='reverse reads file in Sanger FASTQ format')
    parser.add_argument('--fasta', dest='fasta', help='fasta input file')
    parser.add_argument('--input_id', dest='input_id', help='sample id')
    parser.add_argument('--region', dest='region', help='region')
    parser.add_argument('--year', dest='year', help='year')
    parser.add_argument('--output', dest='output', help='output report json file')
    parser.add_argument('--virulotypes', dest='virulotypes', help='strain virulotypes')
    parser.add_argument('--amrgenes', dest='amrgenes', help='strain AMR genes')
    parser.add_argument('--seqtype', dest='seqtype', help='MLST 7 genes')
    parser.add_argument('--samplereport', dest='samplereport', help='sample report')
    args = parser.parse_args()

    os.symlink(args.fasta, 'input.fasta')
    os.symlink(args.input1, 'input_1.fq')
    is_fastq = False
    if get_filetype(args.input1) == "fastq":
        is_fastq = True
    if args.input2:
        os.symlink(args.input2, 'input_2.fq')
    # AMRGENES (only if fastq are provided)
    if is_fastq:
        subprocess.run("abricate --db ncbi input.fasta > " + args.amrgenes, shell=True)
    else:
        Path(args.amrgenes).touch()
    # VIRULOTYPER (only if fastq are provided)
    if is_fastq:
        if args.input2:
            subprocess.run("perl " + TOOL_DIR + "/scripts/patho_typing.pl 'python " + TOOL_DIR + "/scripts/patho_typing.py -s Listeria monocytogenes -f input_1.fq input_2.fq -o output_dir -j 4 --minGeneCoverage 90 --minGeneIdentity 90 --minGeneDepth 15'", shell=True)
        else:
            subprocess.run("perl " + TOOL_DIR + "/scripts/patho_typing.pl 'python " + TOOL_DIR + "/scripts/patho_typing.py -s Listeria monocytogenesi -f input_1.fq -o output_dir -j 4 --minGeneCoverage 90 --minGeneIdentity 90 --minGeneDepth 15'", shell=True)
        subprocess.run("(head -n 1 pathotyper_rep_tot_tab && tail -n +2 pathotyper_rep_tot_tab | sort -k 2rn) > virulotyper", shell=True)
    else:
        Path("virulotyper").touch()
    # SEQUENCETYPER
    subprocess.run("mlst --legacy --scheme listeria_2 input.fasta | cut -f3,4,5,6,7,8,9,10 > " + args.seqtype, shell=True)
    sequence_typing = openFileAsTable(args.seqtype)
    # LISSEROTYPER
    subprocess.run("python LisSero.py input.fasta > output_tab", shell=True)
    sero_typing = openFileAsTable("output_tab")

    try:
        report_data = {}
        report = open(args.output, 'w')
        # write JSON
        report_data["information_name"] = args.input_id
        report_data["region"] = args.region
        report_data["year"] = args.year
        # write results
        if len(sero_typing) == 0:
            report_data["serotype_serogroup"] = "?"
        else:
            report_data["serotype_serogroup"] = sero_typing[1][1]
        if len(sero_typing) == 0:
            report_data["serotype_amplicons"] = "?"
        else:
            report_data["serotype_amplicons"] = getAmplicons(sero_typing)
        report_data["mlst_CC"] = "?"
        report_data["mlst_lineage"] = "?"
        if len(sequence_typing) < 2:
            report_data["mlst_ST"] = "ST?"
        elif sequence_typing[1][1] == "failed":
            report_data["mlst_ST"] = "ST?"
        elif sequence_typing[1][0] == "-":
            report_data["mlst_ST"] = "ST?"
        else:
            report_data["mlst_ST"] = "ST" + sequence_typing[1][0]
            report_data["mlst_CC"], report_data["mlst_lineage"] = getCCandLineage(sequence_typing[1][0])
        subprocess.run("cat virulotyper > " + args.virulotypes, shell=True)
        subprocess.run("sort virulotyper | awk '$2>50 && !seen[substr($1, 1, index($1, \"_\")-1)]++ { printf(\"%s%s\",sep,substr($1, 1, index($1, \"_\")-1));sep=\", \" }END{print \"\"}' > virulotyper_all", shell=True)
        
        if is_fastq:
            with open('virulotyper_all') as viruall:
                virulotypes_all = viruall.readline().strip() 
        else:
            virulotypes_all = "ND"
        report_data["virulotypes_all"] = virulotypes_all.replace(", ,", ",")

        sequence_qc_result = True
        if "?" in sequence_typing or "-" in sequence_typing:
            sequence_qc_result = False
        if sequence_qc_result:
            report_data['qc_status'] = "Passed"
            report_data['qc_messages'] = "Passed."
        else:
            report_data['qc_status'] = "Failed"
            report_data['qc_messages'] = "Accepted for outbreak investigation."
    finally:
        report.write(json.dumps(report_data))
        report.close()

if __name__ == "__main__":
    __main__()


