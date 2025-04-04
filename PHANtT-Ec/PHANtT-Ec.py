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
        subprocess.run("abricate --db resfinder input.fasta > " + args.amrgenes, shell=True)
    else:
        Path(args.amrgenes).touch()
    # VIRULOTYPER (only if fastq are provided)
    if is_fastq:
        if args.input2:
            subprocess.run("perl " + TOOL_DIR + "/scripts/patho_typing.pl 'python " + TOOL_DIR + "/scripts/patho_typing.py -s Escherichia coli -f input_1.fq input_2.fq -o output_dir -j 4 --minGeneCoverage 90 --minGeneIdentity 90 --minGeneDepth 15'", shell=True)
        else:
            subprocess.run("perl " + TOOL_DIR + "/scripts/patho_typing.pl 'python " + TOOL_DIR + "/scripts/patho_typing.py -s Escherichia coli -f input_1.fq -o output_dir -j 4 --minGeneCoverage 90 --minGeneIdentity 90 --minGeneDepth 15'", shell=True)
        subprocess.run("(head -n 1 pathotyper_rep_tot_tab && tail -n +2 pathotyper_rep_tot_tab | sort -k 2rn) > virulotyper", shell=True)
    else:
        Path("virulotyper").touch()
    # SHIGATOXIN TYPER (only if fastq are provided)
    if is_fastq:
        os.symlink(os.popen("which trimmomatic").read().strip()+".jar", 'trimmomatic.jar')
        if args.input2:
            # TRIMMING
            subprocess.run("java ${_JAVA_OPTIONS:--Xmx8G} -jar trimmomatic.jar PE -threads ${GALAXY_SLOTS:-6} -phred33 input_1.fq input_2.fq trimmed1 trimmed1unpaired trimmed2 trimmed2unpaired SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:36", shell=True)
            # ASSEMBLY
            subprocess.run("sh " + TOOL_DIR + "/scripts/stx_subtype_pe.sh " + TOOL_DIR + " trimmed1 trimmed2 input.fasta", shell=True)
        else:
            # TRIMMING
            subprocess.run("java ${_JAVA_OPTIONS:--Xmx8G} -jar trimmomatic.jar SE -threads ${GALAXY_SLOTS:-6} -phred33 input_1.fq trimmed1 SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:55", shell=True)
            # ASSEMBLY
            subprocess.run("sh " + TOOL_DIR + "/scripts/stx_subtype_se.sh " + TOOL_DIR + " trimmed1 input.fasta", shell=True)
        # SHIGATOXIN SEQUENCE SEARCH
        subprocess.run("sh " + TOOL_DIR + "/scripts/stx_subtype_fa.sh " + TOOL_DIR + " stx.fasta", shell=True)
    else:
        shutil.copy(args.input1, "input.fasta")
        Path("shigatoxin_fc").touch()
    # SEQUENCETYPER
    subprocess.run("mlst --legacy --scheme ecoli_achtman_4 input.fasta | cut -f3,4,5,6,7,8,9,10 > " + args.seqtype, shell=True)
    sequence_typing = openFileAsTable(args.seqtype)
    # SEROTYPER O&H
    if is_fastq:
        if args.input2:
            subprocess.run("sh " + TOOL_DIR + "/scripts/serotype.sh " + TOOL_DIR + " y input_1.fq input_2.fq input.fasta", shell=True)
        else:
            subprocess.run("sh " + TOOL_DIR + "/scripts/serotype.sh " + TOOL_DIR + " n input_1.fq xxx input.fasta", shell=True)
    else:
        subprocess.run("sh " + TOOL_DIR + "/scripts/serotype.sh " + TOOL_DIR + " 0 xxx xxx input.fasta", shell=True)
    # SEROTYPER O
    subprocess.run("awk -F '\t' '$4>800 { print $2 FS $3 FS $4 FS $16 }' serogroup_O | sort -nrk 2 -nrk 3 > serogroup_O_fc", shell=True)
    subprocess.run("awk -F , '!seen[$0]++' serogroup_O_fc > serogroup_O_fcd", shell=True)
    sero_typing_o = openFileAsTable("serogroup_O_fcd")
    # SEROTYPER H
    subprocess.run("awk -F '\t' '$4>800 { print $2 FS $3 FS $4 FS $16 }' serogroup_H | sort -nrk 2 -nrk 3 > serogroup_H_fc", shell=True)
    subprocess.run("awk -F , '!seen[$0]++' serogroup_H_fc > serogroup_H_fcd", shell=True)
    sero_typing_h = openFileAsTable("serogroup_H_fcd")
    try:
        report_data = {}
        report = open(args.output, 'w')
        # write JSON
        report_data["information_name"] = args.input_id
        report_data["region"] = args.region
        report_data["year"] = args.year
        # write results
        if len(sero_typing_o) == 0:
            report_data["serotype_o"] = "O?"
        else:
            report_data["serotype_o"] = sero_typing_o[0][0][sero_typing_o[0][0].rfind("O"):]
        if len(sero_typing_h) == 0:
            report_data["serotype_h"] = "H?"
        else:
            report_data["serotype_h"] = sero_typing_h[0][0][sero_typing_h[0][0].rfind("H"):]
        if len(sequence_typing) < 2:
            report_data["mlst_ST"] = "ST?"
        elif sequence_typing[1][1] == "failed":
            report_data["mlst_ST"] = "ST?"
        elif sequence_typing[1][0] == "-":
            report_data["mlst_ST"] = "ST?"
        else:
            report_data["mlst_ST"] = "ST" + sequence_typing[1][0]
        subprocess.run("cat virulotyper > " + args.virulotypes, shell=True)
        subprocess.run("sort virulotyper | awk '/eae_|stx1._|stx2._|ehxa_/ && $2>50 && !seen[substr($1, 1, index($1, \"_\")-1)]++ { printf(\"%s%s\",sep,substr($1, 1, index($1, \"_\")-1));sep=\", \" }END{print \"\"}' > virulotyper_rep", shell=True)
        subprocess.run("sort virulotyper | awk '$2>50 && !seen[substr($1, 1, index($1, \"_\")-1)]++ { printf(\"%s%s\",sep,substr($1, 1, index($1, \"_\")-1));sep=\", \" }END{print \"\"}' > virulotyper_all", shell=True)
        
        if is_fastq:
            with open('virulotyper_rep') as virurep:
                virulotype_eae = "-"
                virulotype_ehxa = "-"
                virulotype_stx1 = "-"
                virulotype_stx2 = "-"
                for line in virurep:
                    if "eae" in line:
                        virulotype_eae = "eae"
                    if "ehxa" in line:
                        virulotype_ehxa = "ehxa"
                    if "stx1" in line:
                        virulotype_stx1 = "stx1"
                    if "stx2" in line:
                        virulotype_stx2 = "stx2"
            with open('virulotyper_all') as viruall:
                virulotypes_all = viruall.readline().strip() 
        else:
            virulotype_eae = "ND"
            virulotype_ehxa = "ND"
            virulotype_stx1 = "ND"
            virulotype_stx2 = "ND"
            virulotypes_all = "ND"
        report_data["virulotype_eae"] = virulotype_eae
        report_data["virulotype_ehxa"] = virulotype_ehxa
        report_data["virulotype_stx1"] = virulotype_stx1
        report_data["virulotype_stx2"] = virulotype_stx2
        report_data["virulotypes_all"] = virulotypes_all.replace(", ,", ",")

        shigatoxin_typing = openFileAsTable("shigatoxin_fc")
        if len(shigatoxin_typing) == 0:
            if is_fastq:
                str_shigatoxin_subtype = "No subtype match found"
            else:
                str_shigatoxin_subtype = "ND"
        else:
            # get corresponding subtypes
            shigatoxin_subtypes = []
            shigatoxin_subtypes_raw = []
            shigatoxin_types = openFileAsTable(TOOL_DIR + "/data/stx_subtypes")
            for subtype in shigatoxin_typing:
                blast_pident_100 = float(subtype[1]) == 100
                if (blast_pident_100):
                    for item in shigatoxin_types:
                        if item[0] == subtype[0] and item[1] not in shigatoxin_subtypes_raw:
                            shigatoxin_subtypes.append(item[1])
                            shigatoxin_subtypes_raw.append(item[1])
           # partial matches
            for subtype in shigatoxin_typing:
                for item in shigatoxin_types:
                    if item[0] == subtype[0] and item[1] not in shigatoxin_subtypes_raw:
                        if item[1][0:4] == "stx1":
                            shigatoxin_subtypes.append(item[1] + "(*)")
                            shigatoxin_subtypes_raw.append(item[1])
                        if item[1][0:4] == "stx2":
                            shigatoxin_subtypes.append(item[1] + "(*)")
                            shigatoxin_subtypes_raw.append(item[1])
            shigatoxin_subtypes.sort()
            str_shigatoxin_subtype = " ".join(shigatoxin_subtypes)
        report_data["shigatoxin_subtype"] = str_shigatoxin_subtype
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


