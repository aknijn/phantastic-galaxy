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
import subprocess
import json

def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', dest='input', help='input file(s)')
    parser.add_argument('--species', dest='species', help='species')
    parser.add_argument('--contamination_json', dest='contamination', help='contamination_json')
    parser.add_argument('--contamination_txt', dest='contamination_txt', help='contamination_txt')
    args = parser.parse_args()

    strContamination = 'No'
    subprocess.run("kmerfinder.py  -i " + args.input + " -db /gfs/data-flow/KmerFinder/bacteria/bacteria.ATG -tax /gfs/data-flow/KmerFinder/bacteria/bacteria.tax -o output", shell=True))
    with open('output/results.txt', 'r') as table_in:
        table_data = [[str(col).rstrip() for col in row.split('\t')] for row in table_in]
    if len(table_data[0]) == 19:
        for row in table_data:
            if row[18] != args.species and row[18] != 'Species' and float(row[6]) > 10.0:
                strContamination = row[18]
    report_json = open(args.contamination_json, 'w')
    if strContamination != 'No':
        report_data = {}
        report_data['qc_status'] = "Failed"
        report_data['qc_messages'] = "Sample contaminated with " + strContamination
        report_json.write(json.dumps(report_data))
    report_json.close()
    report_txt = open(args.contamination_txt, 'w')
    report_txt.write(strContamination)
    report_txt.close()

if __name__ == "__main__":
    __main__()
