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
import subprocess
import json
import os

TOOL_DIR = os.path.dirname(os.path.abspath(__file__))

def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', dest='input', help='input file(s)')
    parser.add_argument('--species', dest='species', help='species')
    parser.add_argument('--contamination_json', dest='contamination_json', help='contamination_json')
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(TOOL_DIR + '/../phantastic.conf')
    dataflowdir = config['fs']['dataflow_path']

    strContamination = 'No'
    subprocess.run("kmerfinder.py  -i " + args.input + " -db " + dataflowdir + "/KmerFinder/bacteria/bacteria.ATG -tax " + dataflowdir + "/KmerFinder/bacteria/bacteria.tax -o output", shell=True)
    with open('output/results.txt', 'r') as table_in:
        table_data = [[str(col).rstrip() for col in row.split('\t')] for row in table_in]
    if len(table_data[0]) == 19:
        for row in table_data:
            # Query_Coverage>4.0 AND Template_length>100.0
            str_species = row[18]
            if str_species == "Shigella flexneri":
                str_species = "Escherichia coli"
            if str_species != args.species and str_species != 'Species' and float(row[5]) > 4.0 and float(row[4]) > 100.0:
                if strContamination != 'No':
                    if not (str_species in strContamination):
                        strContamination = strContamination + ", " + str_species
                else:
                    strContamination = str_species
    report_json = open(args.contamination_json, 'w')
    if strContamination != 'No':
        report_data = {}
        report_data['qc_status'] = "Failed"
        report_data['qc_messages'] = "Sample contaminated with " + strContamination
        report_json.write(json.dumps(report_data))
    report_json.close()

if __name__ == "__main__":
    __main__()
